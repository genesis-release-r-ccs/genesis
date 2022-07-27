#include <config.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <exception>
#include <mpi.h>
#include <bagel.h>
#include <json.h>

extern "C" {
  const char *qsimulate_interface(MPI_Fint*, const int*, const char*,
                              const int, const int*, const double*, const double*,
                              double*, double*, const double*, bool*);
}

static const char *qsimulate_interface_(MPI_Fint& fcomm, const int world_rank, const char* finput,
                              const int natoms, const int* const labels, const double* const coord, const double* const charges,
                              double* const force, double* const qmcharges, const double* const born_radii, bool& error);

const char *qsimulate_interface(MPI_Fint* fcomm, const int* world_rank, const char* input,
                              const int natoms, const int * labels, const double* coord, const double* charges,
                              double* force, double* qmcharges, const double* born_radii, bool* error) {

  return qsimulate_interface_(*fcomm, *world_rank, input,
                       natoms, labels, coord, charges, 
                       force, qmcharges, born_radii, *error);
}

using namespace std;

// quick look up table for elements
static const vector<string> elem = {"Q","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rr","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
    "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
    "Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo"};


shared_ptr<const void> saveref;
vector<char> output_buffer;

static const char *qsimulate_interface_(MPI_Fint& fcomm, const int world_rank, const char* finput,
                                 const int natoms, const int* const labels, const double* const coord, const double* const charges,
                                 double* const force, double* const qmcharges,
                                 const double* const born_radii, bool& error) {

  MPI_Comm comm = MPI_Comm_f2c(fcomm);

  int size, rank;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // assuming that the input will not change in the course of the calculation
  const string input(finput);
  stringstream ss(input);
  bagel::Json jinp;
  ss >> jinp;

  if (jinp.find("bagel") == jinp.end() || !jinp["bagel"].is_array() || jinp["bagel"].size() < 2)
    throw runtime_error("key \"bagel\" not found");

  if (jinp["bagel"][0].find("title") == jinp["bagel"][0].end() || jinp["bagel"][0].at("title") != "molecule")
    throw runtime_error("molecule block is not at the beginning of the input");

  // try to find "basis_link"
  auto basis_link = jinp["bagel"][0].find("basis_link");
  auto df_basis_link = jinp["bagel"][0].find("df_basis_link");
  const string blink = basis_link != jinp["bagel"][0].end() ? *basis_link : "";
  const string dlink = df_basis_link != jinp["bagel"][0].end() ? *df_basis_link : "";

  bagel::Json geometry = bagel::Json::array();
  int ncharges = 0;
  for (int i = 0; i != natoms; ++i) {
    if (labels[i] == 0) {
      ++ncharges;
      continue;
    }
    bagel::Json atom = bagel::Json({});
    atom["atom"] = elem.at(abs(labels[i]));
    atom["xyz"] = {coord[3*i], coord[3*i+1], coord[3*i+2]};

    // link atoms
    if (labels[i] < 0 and !blink.empty())
      atom["basis"] = blink;
    if (labels[i] < 0 and !dlink.empty())
      atom["df_basis"] = dlink;

    geometry.push_back(atom);
  }
  jinp["bagel"][0]["geometry"] = geometry;

  vector<double> extcharges(ncharges * 4);
  ncharges = 0;
  for (int i = 0; i != natoms; ++i) {
    if (labels[i] != 0)
      continue;
   extcharges[0 + 4 * ncharges] = coord[0 + 3 * i];
   extcharges[1 + 4 * ncharges] = coord[1 + 3 * i];
   extcharges[2 + 4 * ncharges] = coord[2 + 3 * i];
   extcharges[3 + 4 * ncharges] = charges[i];
    ++ncharges;
  }

  if (ncharges != 0) {
    jinp["bagel"][0]["charges"] = bagel::Json(move(extcharges));
  }

  // For GBSA
  if (born_radii) {
    vector<double> qm_born_radii(natoms - ncharges);
    vector<double> ext_born_radii(ncharges);
    size_t cnt_i = 0;
    size_t cnt_j = 0;
    for (int i = 0; i != natoms; ++i) {
      if (labels[i] == 0) {
        ext_born_radii[cnt_i++] = born_radii[i];
      } else {
        qm_born_radii[cnt_j++] = born_radii[i];
      }
    }
    assert(cnt_i == ncharges);

    for (auto& block : jinp["bagel"]) {
      if (block["title"] == "force") {
        for (auto& subblock : block["method"]) {
          if (subblock["title"] == "dftb" || subblock["title"] == "xtb") {
            subblock["born_rad"] = bagel::Json(move(qm_born_radii));
            subblock["external_born_rad"] = bagel::Json(move(ext_born_radii));
          }
        }
      }
      // For single point calculations
      if (block["title"] == "dftb" || block["title"] == "xtb") {
        block["born_rad"] = bagel::Json(move(qm_born_radii));
        block["external_born_rad"] = bagel::Json(move(ext_born_radii));
      }
    }
  }

  shared_ptr<const bagel::Json> out_c;
  error = false;
  try {
    tie(out_c, saveref) = bagel::run_bagel_from_json(make_shared<bagel::Json>(move(jinp)), &comm, world_rank, saveref);
  } catch (const exception& e) {
    cout << "exception was thrown in QSimulate: " << e.what() << endl;
    out_c = make_shared<const bagel::Json>();
    error = true;
  }
  shared_ptr<bagel::Json> out = const_pointer_cast<bagel::Json>(out_c);

  bool no_gradient = false;
  if (!error) {
    if (out->find("energy") == out->end()) {
      cout << "energy is not found in QSimulate" << endl;
      error = true;
    }
    if (out->find("gradient") == out->end()) {
      cout << "gradient is not found in QSimulate" << endl;
      // Not finding a gradient does not throw an error.
      // This is to allow single-shot calculations with genesis
      no_gradient = true;
    }
  }

  if (!error) {
    if (no_gradient) {
      double* cptr = qmcharges;
      for (auto& charge : out->at("charges"))
        *cptr++ = static_cast<double>(charge);
    } else {
      // retrieving nuclear force, and optionally qm charges
      double* ptr = force;
      double* cptr = qmcharges;
      for (auto& iatom : out->at("gradient")) {
        for (auto& xyz : iatom.at("xyz"))
          *ptr++ = xyz;
        *cptr++ = iatom.find("charge") != iatom.end() ? static_cast<double>(iatom.at("charge")) : 0.0;
      }
      out->erase("gradient");

      if (out->find("external_gradient") != out->end()) {
        const vector<double>& egrad = out->at("external_gradient");
        assert(egrad.size() == ncharges * 3);
        for (int i = 0; i != ncharges; ++i) {
          *ptr++ = egrad[0 + 3*i];
          *ptr++ = egrad[1 + 3*i];
          *ptr++ = egrad[2 + 3*i];
          *cptr++ = 0.0;
        }
        out->erase("external_gradient");
      }
    }
  }

  string output = out->dump();
  output_buffer.clear();
  std::copy(output.begin(), output.end() + 1, std::back_inserter(output_buffer));

  return output_buffer.data();

}

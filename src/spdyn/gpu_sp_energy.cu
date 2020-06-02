#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include "gpu_common.h"

#include "../config.h"
#include "cuda.h"

#include "nvToolsExt.h"

// #define DEBUG
// #define DO_VERIFY

#define S_ERROR  1e-3
#define D_ERROR  1e-6

#if defined(_SINGLE)
#define REAL  float
#define REALI float
#define ERROR  S_ERROR
#elif defined(_MIXED)
#define REAL  float
#define REALI double
#define ERROR  S_ERROR
#else
#define REAL  double
#define REALI double
#define ERROR  D_ERROR
#endif

typedef unsigned char  uchar;

//
#define coord_pbc(X,Y,Z)          _coord_pbc           [CALIDX3((X)-1,MaxAtom, (Y)-1,3, (Z)-1,ncel_max)]
#define force(X,Y,Z)              _force               [CALIDX3((X)-1,MaxAtom, (Y)-1,3, (Z)-1,ncel_max)]
#define coord(X,Y,Z)              _coord               [CALIDX3((X)-1,3, (Y)-1,MaxAtom, (Z)-1,ncel_max)]
#define trans1(X,Y,Z)             _trans1              [CALIDX3((X)-1,3, (Y)-1,MaxAtom, (Z)-1,ncel_max)]

#define cell_move(X,Y,Z)          _cell_move           [CALIDX3((X)-1,3, (Y)-1,ncel_max, (Z)-1,ncel_max)]
#define charge(Y,Z)               _charge              [CALIDX2((Y)-1,MaxAtom, (Z)-1,ncel_max)]
#define atmcls(Y,Z)               _atmcls              [CALIDX2((Y)-1,MaxAtom, (Z)-1,ncel_max)]
#define natom(Z)                  _natom               [(Z)-1]
#define nonb_lj12(Y,Z)            _nonb_lj12           [CALIDX2((Y)-1,num_atom_cls, (Z)-1,num_atom_cls)]
#define nonb_lj6(Y,Z)             _nonb_lj6            [CALIDX2((Y)-1,num_atom_cls, (Z)-1,num_atom_cls)]
#define nonb_lj6_factor(Z)        _nonb_lj6_factor     [(Z)-1]
#define table_ene(Z)              _table_ene           [(Z)-1]
#define table_grad(Z)             _table_grad          [(Z)-1]
#define univ_cell_pairlist1(Y,Z)  _univ_cell_pairlist1 [CALIDX2((Y)-1,2, (Z)-1,univ_maxcell)]
#define univ_mask2(Y,Z)           _univ_mask2          [CALIDX2((Y)-1,univ_mask2_size, (Z)-1,univ_ncell_near)]
#define univ_ix_list(Y,Z)         _univ_ix_list        [CALIDX2((Y)-1,MaxAtom, (Z)-1,ncel_max)]
#define univ_iy_list(Y,Z)         _univ_iy_list        [CALIDX2((Y)-1,MaxAtom, (Z)-1,ncel_max)]
#define univ_ix_natom(Z)          _univ_ix_natom       [(Z)-1]
#define univ_iy_natom(Z)          _univ_iy_natom       [(Z)-1]
#define univ_ij_sort_list(Z)      _univ_ij_sort_list   [(Z)-1]
#define virial_check(Y,Z)         _virial_check        [CALIDX2((Y)-1,ncel_max, (Z)-1,ncel_max)]

#define ene_virial(Z)      _ene_virial[(Z)-1]
#define ene_viri_mid(Y,Z)  _ene_viri_mid[CALIDX2((Y)-1,5, (Z)-1,univ_maxcell)]
#define force_local(Z)     _force_local[(Z)-1]
#define virial(Z)          _virial[(Z)-1]
#define force_ix_local(Z)  _force_ix_local[(Z)-1]

/* for debug */
#define TMP_univ_ix_list(Y,Z)      tmp_univ_ix_list        [CALIDX2((Y)-1,MaxAtom, (Z)-1,ncel_max)]
#define TMP_univ_iy_list(Y,Z)      tmp_univ_iy_list        [CALIDX2((Y)-1,MaxAtom, (Z)-1,ncel_max)]
#define TMP_univ_ix_natom(Z)       tmp_univ_ix_natom       [(Z)-1]
#define TMP_univ_iy_natom(Z)       tmp_univ_iy_natom       [(Z)-1]

/* */
#define WARP_RSUM_12(X)    { X+=__shfl_xor(X,1,32); X+=__shfl_xor(X,2,32); }
#define WARP_RSUM_345(X)   { X+=__shfl_xor(X,4,32); X+=__shfl_xor(X,8,32); X+=__shfl_xor(X,16,32); }
#define WARP_RSUM_12345(X) { X+=__shfl_xor(X,1,32); X+=__shfl_xor(X,2,32); X+=__shfl_xor(X,4,32); X+=__shfl_xor(X,8,32); X+=__shfl_xor(X,16,32); }


static REAL    *dev_coord_pbc = NULL;
static REAL    *dev_force = NULL;
static REAL    *dev_force_recip = NULL;
static double  *dev_ene_virial = NULL;
static REALI   *dev_coord = NULL;
static REAL    *dev_trans1 = NULL;
static char    *dev_cell_move = NULL;
static REAL    *dev_charge = NULL;
static int     *dev_atmcls = NULL;
static int     *dev_natom = NULL;
static REAL    *dev_nonb_lj12 = NULL;
static REAL    *dev_nonb_lj6 = NULL;
static REAL    *dev_nonb_lj6_factor = NULL;
static REAL    *dev_table_ene = NULL;
static REAL    *dev_table_grad = NULL;
static int     *dev_univ_cell_pairlist1 = NULL;
static int     *dev_univ_ix_natom = NULL;
static uchar   *dev_univ_ix_list = NULL;
static int     *dev_univ_iy_natom = NULL;
static uchar   *dev_univ_iy_list = NULL;
static int     *dev_univ_ij_sort_list = NULL;
static char    *dev_virial_check = NULL;
static double  *dev_ene_viri_mid = NULL;
static char    *dev_univ_mask2 = NULL;
static char    *dev_pack_univ_mask2 = NULL;

static size_t  size_coord_pbc = 0;
static size_t  size_force = 0;
static size_t  size_ene_virial = 0;
static size_t  size_coord = 0;
static size_t  size_trans1 = 0;
static size_t  size_cell_move = 0;
static size_t  size_charge = 0;
static size_t  size_atmcls = 0;
static size_t  size_natom = 0;
static size_t  size_nonb_lj12 = 0;
static size_t  size_nonb_lj6 = 0;
static size_t  size_nonb_lj6_factor = 0;
static size_t  size_table_ene = 0;
static size_t  size_table_grad = 0;
static size_t  size_univ_cell_pairlist1 = 0;
static size_t  size_univ_ix_natom = 0;
static size_t  size_univ_ix_list = 0;
static size_t  size_univ_iy_natom = 0;
static size_t  size_univ_iy_list = 0;
static size_t  size_univ_ij_sort_list = 0;
static size_t  size_virial_check = 0;
static size_t  size_ene_viri_mid = 0;
static size_t  size_univ_mask2 = 0;
static size_t  max_size_univ_mask2 = 0;

/* for debug */
static REAL    *tmp_force = NULL;
static uchar   *tmp_univ_ix_list = NULL;
static uchar   *tmp_univ_iy_list = NULL;
static int     *tmp_univ_ix_natom = NULL;
static int     *tmp_univ_iy_natom = NULL;

/* for recip */
#define qdf_real(Z)        _qdf_real[(Z)-1]
#define qdf(X,Y,Z)         _qdf_real[CALIDX3((X)-1,nlocalx, (Y)-1,nlocaly, (Z)-1,nlocalz)]
#define vi(X,Y,Z)          _vi[CALIDX3((X)-1,3, (Y)-1,MaxAtom, (Z)-1,ncell)]
#define bsc(V,X,Y,Z)       _bsc[CALIDX4((V)-1,4, (X)-1,3, (Y)-1,MaxAtom, (Z)-1,ncell)]
#define bscd(V,X,Y,Z)      _bscd[CALIDX4((V)-1,4, (X)-1,3, (Y)-1,MaxAtom, (Z)-1,ncell)]
#define r_scale(Z)         _r_scale[(Z)-1]
#define ngrid(Z)           _ngrid[(Z)-1]

#define TMP_force(X,Y,Z)   tmp_force[CALIDX3((X)-1,MaxAtom, (Y)-1,3, (Z)-1,ncell)]
#define TMP_qdf_real(Z)    tmp_qdf_real[(Z)-1]
#define TMP_qdf(X,Y,Z)     tmp_qdf_real[CALIDX3((X)-1,nlocalx, (Y)-1,nlocaly, (Z)-1,nlocalz)]
#define TMP_vi(X,Y,Z)      tmp_vi[CALIDX3((X)-1,3, (Y)-1,MaxAtom, (Z)-1,ncell)]
#define TMP_bsc(V,X,Y,Z)   tmp_bsc[CALIDX4((V)-1,4, (X)-1,3, (Y)-1,MaxAtom, (Z)-1,ncell)]
#define TMP_bscd(V,X,Y,Z)  tmp_bscd[CALIDX4((V)-1,4, (X)-1,3, (Y)-1,MaxAtom, (Z)-1,ncell)]

static size_t  size_qdf_real = 0;
static size_t  size_vi = 0;
static size_t  size_bsc = 0;
static size_t  size_bscd = 0;
static size_t  size_r_scale = 0;
static size_t  size_ngrid = 0;

static REAL    *dev_qdf_real = NULL;
static int     *dev_vi = NULL;
static REAL    *dev_bsc = NULL;
static REAL    *dev_bscd = NULL;
static REAL    *dev_r_scale = NULL;
static int     *dev_ngrid = NULL;

static REAL    *tmp_qdf_real = NULL;
static int     *tmp_vi = NULL;
static REAL    *tmp_bsc = NULL;
static REAL    *tmp_bscd = NULL;

/* */
extern __shared__ REAL  smem[];

#define NUM_STREAM  4
#define NUM_EVENT   4
static cudaStream_t  stream[NUM_STREAM];
static cudaEvent_t   event[NUM_EVENT];

#define NUM_P_STREAM  1
static cudaStream_t  p_stream[NUM_STREAM];
static struct cudaDeviceProp  prop;

/* flags */
static int flag_build_pairlist_on_GPU = 0;

/* */
static int *dummy_buf = NULL;

// #define SMEM_SIZE (48*1024)
#define SMEM_SIZE (prop.sharedMemPerBlock)

/*
 *
 */
void gpu_init();

/*
 *
 */
void show_size()
{
#ifdef DEBUG
    printf( "size_coord_pbc          : %9d\n", size_coord_pbc           );
    printf( "size_force              : %9d\n", size_force               );
    printf( "size_ene_virial         : %9d\n", size_ene_virial          );
    printf( "size_coord              : %9d\n", size_coord               );
    printf( "size_trans1             : %9d\n", size_trans1              );
    printf( "size_cell_move          : %9d\n", size_cell_move           );
    printf( "size_charge             : %9d\n", size_charge              );
    printf( "size_atmcls             : %9d\n", size_atmcls              );
    printf( "size_natom              : %9d\n", size_natom               );
    printf( "size_nonb_lj12          : %9d\n", size_nonb_lj12           );
    printf( "size_nonb_lj6           : %9d\n", size_nonb_lj6            );
    printf( "size_nonb_lj6_factor    : %9d\n", size_nonb_lj6_factor     );
    printf( "size_table_ene          : %9d\n", size_table_ene           );
    printf( "size_table_grad         : %9d\n", size_table_grad          );
    printf( "size_univ_cell_pairlist1: %9d\n", size_univ_cell_pairlist1 );
    printf( "size_univ_ix_natom      : %9d\n", size_univ_ix_natom       );
    printf( "size_univ_ix_list       : %9d\n", size_univ_ix_list        );
    printf( "size_univ_iy_natom      : %9d\n", size_univ_iy_natom       );
    printf( "size_univ_iy_list       : %9d\n", size_univ_iy_list        );
    printf( "size_univ_ij_sort_list  : %9d\n", size_univ_ij_sort_list   );
    printf( "size_virial_check       : %9d\n", size_virial_check        );
    printf( "size_ene_viri_mid       : %9d\n", size_ene_viri_mid        );
    printf( "size_univ_mask2         : %9d\n", size_univ_mask2          );
    printf( "max_size_univ_mask2     : %9d\n", max_size_univ_mask2      );
#endif
}

/*
 * do cudaMalloc() when necessary
 */
cudaError_t cudaMalloc_WN( void** devPtr, size_t size )
{
    if ( *devPtr != NULL ) return cudaSuccess;
    return cudaMalloc( devPtr, size );
}

/*
 *
 */
void gpu_init_buffer(
    const REAL   *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    const REAL   *_force,               // ( 1:MaxAtom, 1:3, 1:ncel_max, 1:nthread )
    const double *_ene_virial,          // ( 2 )
    const REALI  *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL   *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char   *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL   *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int    *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int    *_natom,               // ( 1:ncel_max )
    const REAL   *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL   *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL   *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL   *_table_ene,           // ( 1:6*cutoff_int )
    const REAL   *_table_grad,          // ( 1:6*cutoff_int )
    const int    *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char   *_univ_mask2,
    const int    *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar  *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int    *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar  *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const char   *_virial_check,        // ( 1:ncel_max, 1:ncel_max )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_ncell_near,
    int  univ_mask2_size
    )
{
    size_coord_pbc           = sizeof(REAL)   * 3 * MaxAtom * ncel_max;
    size_force               = sizeof(REAL)   * 3 * MaxAtom * ncel_max;
    size_ene_virial          = sizeof(double) * 5;
    size_coord               = sizeof(REALI) * 3 * MaxAtom * ncel_max;
    size_trans1              = sizeof(REAL)   * 3 * MaxAtom * ncel_max;
    size_cell_move           = sizeof(char)   * 3 * ncel_max * ncel_max;
    size_charge              = sizeof(REAL)   * MaxAtom * ncel_max;
    size_atmcls              = sizeof(int )   * MaxAtom * ncel_max;
    size_natom               = sizeof(int )   * ncel_max;
    size_nonb_lj12           = sizeof(REAL)   * num_atom_cls * num_atom_cls;
    size_nonb_lj6            = sizeof(REAL)   * num_atom_cls * num_atom_cls;
    size_nonb_lj6_factor     = sizeof(REAL)   * num_atom_cls;
    size_table_ene           = sizeof(REAL)   * 6*cutoff_int;
    size_table_grad          = sizeof(REAL)   * 6*cutoff_int;
    size_univ_cell_pairlist1 = sizeof(int )   * 2 * univ_maxcell;
    size_univ_ix_list        = sizeof(uchar)  * MaxAtom * univ_maxcell1;
    size_univ_iy_list        = sizeof(uchar)  * MaxAtom * univ_maxcell1;
    size_univ_ix_natom       = sizeof(int)    * univ_maxcell1;
    size_univ_iy_natom       = sizeof(int)    * univ_maxcell1;
    size_univ_ij_sort_list   = sizeof(int )   * univ_maxcell1;
    size_ene_viri_mid        = sizeof(double) * 5 * univ_maxcell;
    size_virial_check        = sizeof(char)   * ncel_max * ncel_max;
    // size_univ_mask2          = sizeof(char)   * univ_mask2_size * univ_ncell_near;
    // max_size_univ_mask2      = size_univ_mask2;

    show_size();

    CUDA_CALL( cudaMalloc_WN( (void**) &dev_coord_pbc          , size_coord_pbc           ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_force              , size_force               ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_force_recip        , size_force               ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_ene_virial         , size_ene_virial          ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_coord              , size_coord               ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_trans1             , size_trans1              ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_cell_move          , size_cell_move           ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_charge             , size_charge              ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_atmcls             , size_atmcls              ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_natom              , size_natom               ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_nonb_lj12          , size_nonb_lj12           ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_nonb_lj6           , size_nonb_lj6            ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_nonb_lj6_factor    , size_nonb_lj6_factor     ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_table_ene          , size_table_ene           ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_table_grad         , size_table_grad          ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_cell_pairlist1, size_univ_cell_pairlist1 ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_ix_natom      , size_univ_ix_natom       ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_ix_list       , size_univ_ix_list        ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_iy_natom      , size_univ_iy_natom       ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_iy_list       , size_univ_iy_list        ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_ij_sort_list  , size_univ_ij_sort_list   ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_ene_viri_mid       , size_ene_viri_mid        ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_virial_check       , size_virial_check        ) );
    // CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_mask2         , size_univ_mask2          ) );

    // CUDA_CALL( cudaDeviceSetSharedMemConfig( cudaSharedMemBankSizeFourByte ) );
    CUDA_CALL( cudaDeviceSetSharedMemConfig( cudaSharedMemBankSizeEightByte ) );
#ifdef DEBUG
    cudaSharedMemConfig pConfig;
    CUDA_CALL( cudaDeviceGetSharedMemConfig ( &pConfig ) );
    printf( "SharedMemConfig: %d\n", pConfig );
#endif
}

/*
 * JJ : send data from host to device before energy/force calculation
 */
void gpu_memcpy_h2d_energy(
    const REAL   *_coord_pbc,
    const REAL   *_force,
    const double  *_ene_virial,
    const REALI  *_coord,
    const REAL   *_trans1,
    const char   *_cell_move,
    const REAL   *_charge,
    const int    *_atmcls,
    const int    *_natom,
    const REAL   *_nonb_lj12,
    const REAL   *_nonb_lj6,
    const REAL   *_nonb_lj6_factor,
    const REAL   *_table_ene,
    const REAL   *_table_grad,
    const int    *_univ_cell_pairlist1,
    const char   *_univ_mask2,
    const int    *_univ_ix_natom,
    const uchar  *_univ_ix_list,
    const int    *_univ_iy_natom,
    const uchar  *_univ_iy_list,
    const int    *_univ_ij_sort_list,
    const char   *_virial_check,
    int  MaxAtom,
    int  univ_maxcell,
    int  univ_ncell_near,
    int  univ_mask2_size,
    const int   update,
    const int   check_virial,
    int   first
    )
{
    CUDA_CALL( cudaMemsetAsync( dev_force,  0, size_force, stream[0] ) );
    CUDA_CALL( cudaMemsetAsync( dev_ene_virial, 0, size_ene_virial, stream[0] ) );

    CUDA_CALL( cudaMemcpyAsync( dev_coord, _coord, size_coord, cudaMemcpyHostToDevice, stream[0] ) );
    CUDA_CALL( cudaEventRecord( event[0], stream[0] ) ); /* dev_coord */

    if ( first ) {
        CUDA_CALL( cudaMemcpy( dev_charge             , _charge             , size_charge             , cudaMemcpyHostToDevice ) );
        CUDA_CALL( cudaMemcpy( dev_atmcls             , _atmcls             , size_atmcls             , cudaMemcpyHostToDevice ) );
        CUDA_CALL( cudaMemcpy( dev_nonb_lj12          , _nonb_lj12          , size_nonb_lj12          , cudaMemcpyHostToDevice ) );
        CUDA_CALL( cudaMemcpy( dev_nonb_lj6           , _nonb_lj6           , size_nonb_lj6           , cudaMemcpyHostToDevice ) );
        CUDA_CALL( cudaMemcpy( dev_nonb_lj6_factor    , _nonb_lj6_factor    , size_nonb_lj6_factor    , cudaMemcpyHostToDevice ) );
        CUDA_CALL( cudaMemcpy( dev_table_ene          , _table_ene          , size_table_ene          , cudaMemcpyHostToDevice ) );
        CUDA_CALL( cudaMemcpy( dev_table_grad         , _table_grad         , size_table_grad         , cudaMemcpyHostToDevice ) );
        CUDA_CALL( cudaMemcpy( dev_natom              , _natom              , size_natom              , cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaEventRecord( event[0], stream[0] ) ); /* dev_natom */

	if ( flag_build_pairlist_on_GPU == 0 ) {
	    CUDA_CALL( cudaMemcpy( dev_trans1             , _trans1             , size_trans1             , cudaMemcpyHostToDevice ) );
	    CUDA_CALL( cudaMemcpy( dev_cell_move          , _cell_move          , size_cell_move          , cudaMemcpyHostToDevice ) );

	    CUDA_CALL( cudaMemcpy( dev_univ_cell_pairlist1, _univ_cell_pairlist1, size_univ_cell_pairlist1, cudaMemcpyHostToDevice ) );
	    CUDA_CALL( cudaMemcpy( dev_univ_ix_list       , _univ_ix_list       , size_univ_ix_list       , cudaMemcpyHostToDevice ) );
	    CUDA_CALL( cudaMemcpy( dev_univ_iy_list       , _univ_iy_list       , size_univ_iy_list       , cudaMemcpyHostToDevice ) );
	    CUDA_CALL( cudaMemcpy( dev_univ_ix_natom      , _univ_ix_natom      , size_univ_ix_natom      , cudaMemcpyHostToDevice ) );
	    CUDA_CALL( cudaMemcpy( dev_univ_iy_natom      , _univ_iy_natom      , size_univ_iy_natom      , cudaMemcpyHostToDevice ) );
	}
        CUDA_CALL( cudaMemcpy( dev_univ_ij_sort_list  , _univ_ij_sort_list  , size_univ_ij_sort_list  , cudaMemcpyHostToDevice ) );
        CUDA_CALL( cudaMemcpy( dev_virial_check       , _virial_check       , size_virial_check       , cudaMemcpyHostToDevice ) );
        // CUDA_CALL( cudaMemcpy( dev_univ_mask2         , _univ_mask2         , size_univ_mask2         , cudaMemcpyHostToDevice ) );

        return;
    }

    if ( update ) {
        CUDA_CALL( cudaMemcpy( dev_charge             , _charge             , size_charge             , cudaMemcpyHostToDevice ) );
        CUDA_CALL( cudaMemcpy( dev_atmcls             , _atmcls             , size_atmcls             , cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy( dev_natom              , _natom              , size_natom              , cudaMemcpyHostToDevice ) );

	if ( flag_build_pairlist_on_GPU == 0 ) {
	    CUDA_CALL( cudaMemcpy( dev_trans1             , _trans1             , size_trans1             , cudaMemcpyHostToDevice ) );

	    // CUDA_CALL( cudaMemcpy( dev_univ_cell_pairlist1, _univ_cell_pairlist1, size_univ_cell_pairlist1, cudaMemcpyHostToDevice ) );
	    CUDA_CALL( cudaMemcpy( dev_univ_ix_list       , _univ_ix_list       , size_univ_ix_list       , cudaMemcpyHostToDevice ) );
	    CUDA_CALL( cudaMemcpy( dev_univ_iy_list       , _univ_iy_list       , size_univ_iy_list       , cudaMemcpyHostToDevice ) );
	    CUDA_CALL( cudaMemcpy( dev_univ_ix_natom      , _univ_ix_natom      , size_univ_ix_natom      , cudaMemcpyHostToDevice ) );
	    CUDA_CALL( cudaMemcpy( dev_univ_iy_natom      , _univ_iy_natom      , size_univ_iy_natom      , cudaMemcpyHostToDevice ) );
	}
        CUDA_CALL( cudaMemcpy( dev_univ_ij_sort_list  , _univ_ij_sort_list  , size_univ_ij_sort_list  , cudaMemcpyHostToDevice ) );

        // size_univ_mask2 = sizeof(char) * univ_mask2_size * univ_ncell_near ;
        // if (max_size_univ_mask2 < size_univ_mask2){
        //    max_size_univ_mask2 = size_univ_mask2;
        //    CUDA_CALL( cudaFree( dev_univ_mask2 ) );
	//    dev_univ_mask2 = NULL;
        //    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_mask2, max_size_univ_mask2 ) );
        //    show_size();
        // }
        // CUDA_CALL( cudaMemcpy( dev_univ_mask2         , _univ_mask2         , size_univ_mask2         , cudaMemcpyHostToDevice ) );

        return;
    }
    else if ( check_virial) {
  	    CUDA_CALL( cudaMemcpy( dev_trans1             , _trans1             , size_trans1             , cudaMemcpyHostToDevice ) );
        return;
    }
}

/*
 * JJ : send data from host to device
 */
void gpu_memcpy_h2d_force(
    const REAL   *_coord_pbc,
    const REAL   *_force,
    const double *_ene_virial,
    const REALI  *_coord,
    const REAL   *_trans1,
    const char   *_cell_move,
    const REAL   *_charge,
    const int    *_atmcls,
    const int    *_natom,
    const REAL   *_nonb_lj12,
    const REAL   *_nonb_lj6,
    const REAL   *_nonb_lj6_factor,
    const REAL   *_table_grad,
    const int    *_univ_cell_pairlist1,
    const char   *_univ_mask2,
    const int    *_univ_ix_natom,
    const uchar  *_univ_ix_list,
    const int    *_univ_iy_natom,
    const uchar  *_univ_iy_list,
    const int    *_univ_ij_sort_list,
    const char   *_virial_check,
    int  MaxAtom,
    int  univ_maxcell,
    int  univ_ncell_near,
    int  univ_mask2_size,
    const int   update,
    int   check_virial,
    cudaStream_t  st
    )
{
    CUDA_CALL( cudaMemsetAsync( dev_force,      0, size_force,      st ) );
    CUDA_CALL( cudaMemsetAsync( dev_ene_virial, 0, size_ene_virial, st ) );

    if ( flag_build_pairlist_on_GPU == 0 || update == 0 ) {
	CUDA_CALL( cudaMemcpyAsync( dev_coord, _coord, size_coord, cudaMemcpyHostToDevice, st ) );
	CUDA_CALL( cudaEventRecord( event[0], st ) ); /* dev_coord */
    }

    if ( update ) {
	CUDA_CALL( cudaMemcpyAsync( dev_charge, _charge, size_charge, cudaMemcpyHostToDevice, st ) );
	CUDA_CALL( cudaMemcpyAsync( dev_atmcls, _atmcls, size_atmcls, cudaMemcpyHostToDevice, st ) );
        CUDA_CALL( cudaMemcpyAsync( dev_natom , _natom , size_natom , cudaMemcpyHostToDevice, st ) );
	CUDA_CALL( cudaEventRecord( event[0], st ) ); /* dev_natom */

	if ( flag_build_pairlist_on_GPU == 0 ) {
	    CUDA_CALL( cudaMemcpyAsync( dev_trans1       , _trans1       , size_trans1       , cudaMemcpyHostToDevice, st ) );
	    CUDA_CALL( cudaMemcpyAsync( dev_univ_ix_list , _univ_ix_list , size_univ_ix_list , cudaMemcpyHostToDevice, st ) );
	    CUDA_CALL( cudaMemcpyAsync( dev_univ_iy_list , _univ_iy_list , size_univ_iy_list , cudaMemcpyHostToDevice, st ) );
	    CUDA_CALL( cudaMemcpyAsync( dev_univ_ix_natom, _univ_ix_natom, size_univ_ix_natom, cudaMemcpyHostToDevice, st ) );
	    CUDA_CALL( cudaMemcpyAsync( dev_univ_iy_natom, _univ_iy_natom, size_univ_iy_natom, cudaMemcpyHostToDevice, st ) );
	}
	CUDA_CALL( cudaMemcpyAsync( dev_univ_ij_sort_list, _univ_ij_sort_list, size_univ_ij_sort_list, cudaMemcpyHostToDevice, st ) );

        // size_univ_mask2 = sizeof(char) * univ_mask2_size * univ_ncell_near ;
        // if (max_size_univ_mask2 < size_univ_mask2){
	//     max_size_univ_mask2 = size_univ_mask2;
	//     CUDA_CALL( cudaFree( dev_univ_mask2 ) );
	//     dev_univ_mask2 = NULL;
	//     CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_mask2, max_size_univ_mask2 ) );
	//     show_size();
        // }
        // CUDA_CALL( cudaMemcpy( dev_univ_mask2         , _univ_mask2         , size_univ_mask2         , cudaMemcpyHostToDevice ) );

	return;
    }
    else if ( check_virial) {
            CUDA_CALL( cudaMemcpyAsync( dev_trans1   , _trans1   , size_trans1   , cudaMemcpyHostToDevice, st ) );
        return;
    }

}

/*
 *
 */
void gpu_init_buffer_pairlist(
    const REAL   *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    const REALI  *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL   *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char   *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const int    *_natom,               // ( 1:ncel_max )
    const int    *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const uchar  *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const uchar  *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int    *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const int    *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    int  MaxAtom,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  univ_maxcell,
    int  univ_maxcell1
    )
{
    size_coord_pbc           = sizeof(REAL)   * 3 * MaxAtom * ncel_max;
    size_coord               = sizeof(REALI)  * 3 * MaxAtom * ncel_max;
    size_trans1              = sizeof(REAL)   * 3 * MaxAtom * ncel_max;
    size_cell_move           = sizeof(char)   * 3 * ncel_max * ncel_max;
    size_natom               = sizeof(int )   * ncel_max;
    size_univ_cell_pairlist1 = sizeof(int )   * 2 * univ_maxcell;
    size_univ_ix_list        = sizeof(uchar)  * MaxAtom * univ_maxcell1;
    size_univ_iy_list        = sizeof(uchar)  * MaxAtom * univ_maxcell1;
    size_univ_ix_natom       = sizeof(int)    * univ_maxcell1;
    size_univ_iy_natom       = sizeof(int)    * univ_maxcell1;

    show_size();

    CUDA_CALL( cudaMalloc_WN( (void**) &dev_coord_pbc          , size_coord_pbc           ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_coord              , size_coord               ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_trans1             , size_trans1              ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_cell_move          , size_cell_move           ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_natom              , size_natom               ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_cell_pairlist1, size_univ_cell_pairlist1 ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_ix_list       , size_univ_ix_list        ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_iy_list       , size_univ_iy_list        ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_ix_natom      , size_univ_ix_natom       ) );
    CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_iy_natom      , size_univ_iy_natom       ) );

    // CUDA_CALL( cudaDeviceSetSharedMemConfig( cudaSharedMemBankSizeFourByte ) );
    CUDA_CALL( cudaDeviceSetSharedMemConfig( cudaSharedMemBankSizeEightByte ) );
#ifdef DEBUG
    cudaSharedMemConfig pConfig;
    CUDA_CALL( cudaDeviceGetSharedMemConfig ( &pConfig ) );
    printf( "SharedMemConfig: %d\n", pConfig );
#endif

    /* debug */
    tmp_univ_ix_list = (uchar*)malloc( size_univ_ix_list );
    tmp_univ_iy_list = (uchar*)malloc( size_univ_iy_list );
    tmp_univ_ix_natom = (int*)malloc( size_univ_ix_natom );
    tmp_univ_iy_natom = (int*)malloc( size_univ_iy_natom );
}

/*
 *
 */
void gpu_memcpy_h2d_pairlist(
    REAL         *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    const REALI  *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL   *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char   *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const int    *_natom,               // ( 1:ncel_max )
    const int    *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    uchar        *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    uchar        *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    int          *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    int          *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    int  first,
    cudaStream_t  st
    )
{
    CUDA_CALL( cudaMemsetAsync( dev_univ_ix_list, 0, size_univ_ix_list, st ) );
    CUDA_CALL( cudaMemsetAsync( dev_univ_iy_list, 0, size_univ_iy_list, st ) );

    CUDA_CALL( cudaMemcpyAsync( dev_coord, _coord, size_coord, cudaMemcpyHostToDevice, st ) );

    CUDA_CALL( cudaMemcpyAsync( dev_trans1, _trans1, size_trans1, cudaMemcpyHostToDevice, st ) );
    CUDA_CALL( cudaMemcpyAsync( dev_natom, _natom, size_natom, cudaMemcpyHostToDevice, st ) );

    if ( first ) {
	CUDA_CALL( cudaMemcpyAsync( dev_cell_move, _cell_move, size_cell_move,
				    cudaMemcpyHostToDevice, st ) );
	CUDA_CALL( cudaMemcpyAsync( dev_univ_cell_pairlist1, _univ_cell_pairlist1, size_univ_cell_pairlist1,
				    cudaMemcpyHostToDevice, st ) );
    }
}

/*
 *
 */
void gpu_memcpy_d2h_pairlist(
    uchar         *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    uchar         *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    int           *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    int           *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    cudaStream_t  st,
    int  dummy
    )
{
    /* if you use pairlist on CPU, the following two d2h memcpy are necessary */
    // CUDA_CALL( cudaMemcpy( _univ_ix_list, dev_univ_ix_list, size_univ_ix_list, cudaMemcpyDeviceToHost ) );
    // CUDA_CALL( cudaMemcpy( _univ_iy_list, dev_univ_iy_list, size_univ_iy_list, cudaMemcpyDeviceToHost ) );

    CUDA_CALL( cudaMemcpyAsync( _univ_ix_natom, dev_univ_ix_natom, size_univ_ix_natom,
				cudaMemcpyDeviceToHost, st ) );
    CUDA_CALL( cudaMemcpyAsync( _univ_iy_natom, dev_univ_iy_natom, size_univ_iy_natom,
				cudaMemcpyDeviceToHost, st ) );
}

/*
 * JJ : return force, energy, and virial valudes from the device to host
 */
void gpu_memcpy_d2h_energy(
    REAL  *_coord_pbc,
    REAL  *_force,
    double *_ene_virial
    )
{
    CUDA_CALL( cudaMemcpy( _coord_pbc,  dev_coord_pbc,  size_coord_pbc,  cudaMemcpyDeviceToHost ) );
    CUDA_CALL( cudaMemcpy( _force,  dev_force,  size_force,  cudaMemcpyDeviceToHost ) );
    CUDA_CALL( cudaMemcpy( _ene_virial, dev_ene_virial, size_ene_virial, cudaMemcpyDeviceToHost ) );
}

/*
 *
 */
void gpu_memcpy_d2h_force(
    REAL  *_coord_pbc,
    REAL  *_force,
    double *_ene_virial,
    int  check_virial,
    cudaStream_t  st
    )
{
    CUDA_CALL( cudaMemcpyAsync( _force, dev_force, size_force, cudaMemcpyDeviceToHost, st ) );

    if ( check_virial != 0 ) {
       CUDA_CALL( cudaMemcpyAsync( _coord_pbc,  dev_coord_pbc,  size_coord_pbc,  cudaMemcpyDeviceToHost, st ) );
       CUDA_CALL( cudaMemcpyAsync( _ene_virial, dev_ene_virial, size_ene_virial, cudaMemcpyDeviceToHost, st ) );
    }
}

/*
 *
 */
#if (CUDA_VERSION < 6600)
#if __CUDA_ARCH__ < 600
#warning CUDA_ARCH < 600
__device__ double atomicAdd( double* address, double val )
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif
#else

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
#warning CUDA_ARCH < 600
__device__ double atomicAdd( double* address, double val )
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif
#endif

/*
 *
 */
#if (CUDA_VERSION < 6050)
__device__ __inline__ double __shfl_xor( double var, int mask, int width )
{
    int  lo, hi;
    asm volatile( "mov.b64 {%0,%1}, %2;" : "=r"(lo), "=r"(hi) : "d"(var) );
    lo = __shfl_xor( lo, mask, width );
    hi = __shfl_xor( hi, mask, width );
    asm volatile( "mov.b64 %0, {%1,%2};" : "=d"(var) : "r"(lo), "r"(hi) );
    return var;
}
#endif

/*
 *
 */
__global__ void kern_compute_force_nonbond_table_linear_univ__coord_pbc(
    REAL        *_coord_pbc,    // ( 1:MaxAtom, 1:3, 1:ncel_max )
    const REALI *_coord,        // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,       // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const int   *_natom,        // ( 1:ncel_max )
    /* */
    const int  MaxAtom,
    const int  ncel_local,
    const int  ncel_bound,
    const int  ncel_max )
{
    const int  id_block = blockIdx.x;
    const int  id_thread = threadIdx.x;
    const int  num_thread = blockDim.x;

    int i = id_block + 1;
    if ( i > ncel_local + ncel_bound )
	return;

    int k;
    for ( k = id_thread; k < natom(i) * 3; k += num_thread ) {
	int n  = (k % 3) + 1;
	int ix = (k / 3) + 1;
	coord_pbc(ix,n,i) = coord(n,ix,i) + trans1(n,ix,i);
    }
}

/*
 *
 */
#if defined(_MIXED) || defined(_SINGLE)
//__launch_bounds__(128,12)  // for SP
#else
//__launch_bounds__(128,8)  // for DP
#endif
__global__ void kern_compute_energy_nonbond_table_linear_univ__energyforce_intra_cell(
    REAL        * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:ncel_max)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _table_ene,           // ( 1:6*cutoff_int )
    const REAL  * _table_grad,          // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_start,
    int  index_end,
    REAL  density,
    REAL  cutoff2
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];
    double val_energy_elec = 0.0;
    double val_energy_evdw = 0.0;

    int  index = index_start + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_end ) return;
    int  univ_ij = index;

    const int  i = univ_ij;
    const int  j = univ_ij;

    int  ix_natom = natom(i);
    int  iy_natom = natom(j);
    if ( ix_natom * iy_natom <= 0 ) return;

    for ( int ix = id_thread_xx + 1; ix <= ix_natom; ix += num_thread_xx ) {

        REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i));
        REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i));
        REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i));
        REAL  iqtmp   = __ldg(&charge(ix,i));
        int   iatmcls = __ldg(&atmcls(ix,i));

        force_local(1) = 0.0;
        force_local(2) = 0.0;
        force_local(3) = 0.0;
        REAL   elec_temp = 0.0;
        REAL   evdw_temp = 0.0;

        for ( int iy = id_thread_xy + 1; iy <= iy_natom; iy += num_thread_xy ) {

            REAL  grad_coef = 0.0;
            REAL  dij1 = 0.0;
            REAL  dij2 = 0.0;
            REAL  dij3 = 0.0;
            REAL  rij2;

            int idx = iy + (ix-1)*univ_natom_max;
            if (univ_mask2(idx,univ_ij)) {

                dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

                if ( rij2 < cutoff2) {

                    rij2 = cutoff2 * density / rij2;

                    REAL  jqtmp   = __ldg(&charge(iy,j));
                    int   jatmcls = __ldg(&atmcls(iy,j));
                    REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
                    REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));

                    int   L  = int(rij2);
                    REAL  R  = rij2 - L;
                    int   L1 = 3*L - 2;

                    REAL  tg0  = __ldg(&table_ene(L1  ));
                    REAL  tg1  = __ldg(&table_ene(L1+1));
                    REAL  tg2  = __ldg(&table_ene(L1+2));
                    REAL  tg3  = __ldg(&table_ene(L1+3));
                    REAL  tg4  = __ldg(&table_ene(L1+4));
                    REAL  tg5  = __ldg(&table_ene(L1+5));
                    REAL  term_lj12 = tg0 + R*(tg3-tg0);
                    REAL  term_lj6  = tg1 + R*(tg4-tg1);
                    REAL  term_elec = tg2 + R*(tg5-tg2);

                    elec_temp += iqtmp*jqtmp*term_elec;
                    evdw_temp += term_lj12*lj12 - term_lj6*lj6;

                    tg0  = __ldg(&table_grad(L1  ));
                    tg1  = __ldg(&table_grad(L1+1));
                    tg2  = __ldg(&table_grad(L1+2));
                    tg3  = __ldg(&table_grad(L1+3));
                    tg4  = __ldg(&table_grad(L1+4));
                    tg5  = __ldg(&table_grad(L1+5));
                    term_lj12 = tg0 + R*(tg3-tg0);
                    term_lj6  = tg1 + R*(tg4-tg1);
                    term_elec = tg2 + R*(tg5-tg2);

                    grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;
                }
            }

            REAL  work1 = grad_coef*dij1;
            REAL  work2 = grad_coef*dij2;
            REAL  work3 = grad_coef*dij3;

            force_local(1) -= work1;
            force_local(2) -= work2;
            force_local(3) -= work3;
        }
        val_energy_elec += elec_temp/2.0;
        val_energy_evdw += evdw_temp/2.0;

        // update energy/force(:,:,i)
        WARP_RSUM_12( force_local(1) );
        WARP_RSUM_12( force_local(2) );
        WARP_RSUM_12( force_local(3) );
        if ( id_thread_xy == 0 ) {
            if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
            if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
            if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
        }
    }

    WARP_RSUM_12345 ( val_energy_elec );
    WARP_RSUM_12345 ( val_energy_evdw );
    if (id_thread_x < 5) {
        int n = id_thread_x + 1;
        if ( n == 1 ) ene_viri_mid(n,index) = val_energy_elec;
        if ( n == 2 ) ene_viri_mid(n,index) = val_energy_evdw;
        if ( n == 3 ) ene_viri_mid(n,index) = 0.0;
        if ( n == 4 ) ene_viri_mid(n,index) = 0.0;
        if ( n == 5 ) ene_viri_mid(n,index) = 0.0;
    }
}


__global__ void kern_compute_energy_nonbond_table_ljpme_univ__energyforce_intra_cell(
    REAL        * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:ncel_max)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  * _table_ene,           // ( 1:6*cutoff_int )
    const REAL  * _table_grad,          // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_start,
    int  index_end,
    REAL  density,
    REAL  cutoff2
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];
    double val_energy_elec = 0.0;
    double val_energy_evdw = 0.0;

    int  index = index_start + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_end ) return;
    int  univ_ij = index;

    const int  i = univ_ij;
    const int  j = univ_ij;

    int  ix_natom = natom(i);
    int  iy_natom = natom(j);
    if ( ix_natom * iy_natom <= 0 ) return;

    REAL inv_cutoff2 = 1.0 / cutoff2;
    REAL inv_cutoff6 = inv_cutoff2 * inv_cutoff2 * inv_cutoff2;

    for ( int ix = id_thread_xx + 1; ix <= ix_natom; ix += num_thread_xx ) {

        REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i));
        REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i));
        REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i));
        REAL  iqtmp   = __ldg(&charge(ix,i));
        int   iatmcls = __ldg(&atmcls(ix,i));
        REAL  lj6_i   = __ldg(&nonb_lj6_factor(iatmcls));

        force_local(1) = 0.0;
        force_local(2) = 0.0;
        force_local(3) = 0.0;
        REAL   elec_temp = 0.0;
        REAL   evdw_temp = 0.0;

        for ( int iy = id_thread_xy + 1; iy <= iy_natom; iy += num_thread_xy ) {

            REAL  grad_coef = 0.0;
            REAL  dij1 = 0.0;
            REAL  dij2 = 0.0;
            REAL  dij3 = 0.0;
            REAL  rij2;

            int idx = iy + (ix-1)*univ_natom_max;
            if (univ_mask2(idx,univ_ij)) {

                dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

                if ( rij2 < cutoff2) {

                    REAL  inv_r2  = 1.0 / rij2;
                    REAL  inv_r6  = inv_r2 * inv_r2 * inv_r2;
                    REAL  term_lj12 = inv_r6 * inv_r6;
                    rij2 = cutoff2 * density * inv_r2;

                    REAL  jqtmp   = __ldg(&charge(iy,j));
                    int   jatmcls = __ldg(&atmcls(iy,j));
                    REAL  lj6_j   = __ldg(&nonb_lj6_factor(jatmcls));
                    REAL  lj6_ij  = lj6_i * lj6_j;
                    REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
                    REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls)) - lj6_ij;

                    int   L  = int(rij2);
                    REAL  R  = rij2 - L;
                    int   L1 = 2*L - 1;

                    REAL  tg0  = __ldg(&table_ene(L1  ));
                    REAL  tg1  = __ldg(&table_ene(L1+1));
                    REAL  tg2  = __ldg(&table_ene(L1+2));
                    REAL  tg3  = __ldg(&table_ene(L1+3));
                    REAL  term_lj6  = tg0 + R*(tg2-tg0);
                    REAL  term_elec = tg1 + R*(tg3-tg1);
                    REAL  term_temp = inv_r6 - inv_cutoff6;

                    elec_temp += iqtmp*jqtmp*term_elec;
                    evdw_temp += term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij;

                    tg0  = __ldg(&table_grad(L1  ));
                    tg1  = __ldg(&table_grad(L1+1));
                    tg2  = __ldg(&table_grad(L1+2));
                    tg3  = __ldg(&table_grad(L1+3));
                    term_lj6  = tg0 + R*(tg2-tg0);
                    term_elec = tg1 + R*(tg3-tg1);
                    term_lj12 = -12.0 * term_lj12 * inv_r2;
                    term_temp = - 6.0 * inv_r6 * inv_r2;

                    grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij + iqtmp*jqtmp*term_elec;
                }
            }

            REAL  work1 = grad_coef*dij1;
            REAL  work2 = grad_coef*dij2;
            REAL  work3 = grad_coef*dij3;

            force_local(1) -= work1;
            force_local(2) -= work2;
            force_local(3) -= work3;
        }
        val_energy_elec += elec_temp/2.0;
        val_energy_evdw += evdw_temp/2.0;

        // update energy/force(:,:,i)
        WARP_RSUM_12( force_local(1) );
        WARP_RSUM_12( force_local(2) );
        WARP_RSUM_12( force_local(3) );
        if ( id_thread_xy == 0 ) {
            if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
            if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
            if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
        }
    }

    WARP_RSUM_12345 ( val_energy_elec );
    WARP_RSUM_12345 ( val_energy_evdw );
    if (id_thread_x < 5) {
        int n = id_thread_x + 1;
        if ( n == 1 ) ene_viri_mid(n,index) = val_energy_elec;
        if ( n == 2 ) ene_viri_mid(n,index) = val_energy_evdw;
        if ( n == 3 ) ene_viri_mid(n,index) = 0.0;
        if ( n == 4 ) ene_viri_mid(n,index) = 0.0;
        if ( n == 5 ) ene_viri_mid(n,index) = 0.0;
    }
}
/*
 *
 */
#if defined(_MIXED) || defined(_SINGLE)
//__launch_bounds__(128,12)  // for SP
#else
//__launch_bounds__(128,8)  // for DP
#endif
__global__ void kern_compute_energy_nonbond_notable_univ__energyforce_intra_cell(
    REAL        * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:ncel_max)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _table_ene,           // ( 1:6*cutoff_int )
    const REAL  * _table_grad,          // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_start,
    int  index_end,
    REAL  density,
    REAL  cutoff2
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];
    double val_energy_elec = 0.0;
    double val_energy_evdw = 0.0;

    int  index = index_start + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_end ) return;
    int  univ_ij = index;

    const int  i = univ_ij;
    const int  j = univ_ij;

    int  ix_natom = natom(i);
    int  iy_natom = natom(j);
    if ( ix_natom * iy_natom <= 0 ) return;

    for ( int ix = id_thread_xx + 1; ix <= ix_natom; ix += num_thread_xx ) {

        REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i));
        REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i));
        REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i));
        REAL  iqtmp   = __ldg(&charge(ix,i));
        int   iatmcls = __ldg(&atmcls(ix,i));

        force_local(1) = 0.0;
        force_local(2) = 0.0;
        force_local(3) = 0.0;
        REAL   elec_temp = 0.0;
        REAL   evdw_temp = 0.0;

        for ( int iy = id_thread_xy + 1; iy <= iy_natom; iy += num_thread_xy ) {

            REAL  grad_coef = 0.0;
            REAL  dij1 = 0.0;
            REAL  dij2 = 0.0;
            REAL  dij3 = 0.0;
            REAL  rij2;

            int idx = iy + (ix-1)*univ_natom_max;
            if (univ_mask2(idx,univ_ij)) {

                dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

                if ( rij2 < cutoff2) {

                    REAL rij2_inv = 1.0 / rij2;
                    rij2 = cutoff2 * density * rij2_inv;

                    REAL  jqtmp   = __ldg(&charge(iy,j));
                    int   jatmcls = __ldg(&atmcls(iy,j));
                    REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
                    REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));

                    int   L  = int(rij2);
                    REAL  R  = rij2 - L;

                    REAL  tg0  = __ldg(&table_ene(L  ));
                    REAL  tg1  = __ldg(&table_ene(L+1));
                    REAL  term_elec = tg0 + R*(tg1-tg0);
                    REAL  term_lj6  = rij2_inv * rij2_inv * rij2_inv;
                    REAL  term_lj12 = term_lj6 * term_lj6;

                    elec_temp += iqtmp*jqtmp*term_elec;
                    evdw_temp += term_lj12*lj12 - term_lj6*lj6;

                    tg0  = __ldg(&table_grad(L  ));
                    tg1  = __ldg(&table_grad(L+1));
                    term_elec = tg0 + R*(tg1-tg0);
                    term_lj12 = -12.0 * term_lj12 * rij2_inv;
                    term_lj6  = - 6.0 * term_lj6  * rij2_inv;

                    grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;
                }
            }

            REAL  work1 = grad_coef*dij1;
            REAL  work2 = grad_coef*dij2;
            REAL  work3 = grad_coef*dij3;

            force_local(1) -= work1;
            force_local(2) -= work2;
            force_local(3) -= work3;
        }
        val_energy_elec += elec_temp/2.0;
        val_energy_evdw += evdw_temp/2.0;

        // update energy/force(:,:,i)
        WARP_RSUM_12( force_local(1) );
        WARP_RSUM_12( force_local(2) );
        WARP_RSUM_12( force_local(3) );
        if ( id_thread_xy == 0 ) {
            if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
            if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
            if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
        }
    }

    WARP_RSUM_12345 ( val_energy_elec );
    WARP_RSUM_12345 ( val_energy_evdw );
    if (id_thread_x < 5) {
        int n = id_thread_x + 1;
        if ( n == 1 ) ene_viri_mid(n,index) = val_energy_elec;
        if ( n == 2 ) ene_viri_mid(n,index) = val_energy_evdw;
        if ( n == 3 ) ene_viri_mid(n,index) = 0.0;
        if ( n == 4 ) ene_viri_mid(n,index) = 0.0;
        if ( n == 5 ) ene_viri_mid(n,index) = 0.0;
    }
}

/*
 *
 */
#if defined(_MIXED) || defined(_SINGLE)
//__launch_bounds__(128,12)  // for SP
#else
//__launch_bounds__(128,9)  // for DP
#endif
__global__ void kern_compute_force_nonbond_table_linear_univ__force_intra_cell(
    const REAL  * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:ncel_max)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _table_grad,          // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near )
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_s,
    int  index_e,
    int  check_virial,
    REAL  density,
    REAL  cutoff2,
    REAL  pairlistdist2
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];

    int  index = index_s + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_e ) return;
    int  univ_ij = index;

    const int  i = univ_ij;
    const int  j = univ_ij;

    int  ix_natom = natom(i);
    int  iy_natom = natom(j);
    if ( ix_natom * iy_natom <= 0 ) return;

    for ( int ix = id_thread_xx + 1; ix <= ix_natom; ix += num_thread_xx ) {

	REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i));
	REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i));
	REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i));
	REAL  iqtmp   = __ldg(&charge(ix,i));
	int   iatmcls = __ldg(&atmcls(ix,i));

	force_local(1) = 0.0;
	force_local(2) = 0.0;
	force_local(3) = 0.0;

	for ( int iy = id_thread_xy + 1; iy <= iy_natom; iy += num_thread_xy ) {

            REAL  grad_coef = 0.0;
            REAL  dij1 = 0.0;
            REAL  dij2 = 0.0;
            REAL  dij3 = 0.0;
            REAL  rij2;

            int idx = iy + (ix-1)*univ_natom_max;

            if (univ_mask2(idx,univ_ij)) {

                dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

                if ( rij2 < cutoff2 && rij2 > 0.0) {

                    rij2 = cutoff2 * density / rij2;

                    REAL  jqtmp   = __ldg(&charge(iy,j));
                    int   jatmcls = __ldg(&atmcls(iy,j));
                    REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
                    REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));

                    int   L  = int(rij2);
                    REAL  R  = rij2 - L;
                    int   L1 = 3*L - 2;
                    REAL  tg0  = __ldg(&table_grad(L1  ));
                    REAL  tg1  = __ldg(&table_grad(L1+1));
                    REAL  tg2  = __ldg(&table_grad(L1+2));
                    REAL  tg3  = __ldg(&table_grad(L1+3));
                    REAL  tg4  = __ldg(&table_grad(L1+4));
                    REAL  tg5  = __ldg(&table_grad(L1+5));
                    REAL  term_lj12 = tg0 + R*(tg3-tg0);
                    REAL  term_lj6  = tg1 + R*(tg4-tg1);
                    REAL  term_elec = tg2 + R*(tg5-tg2);

                    grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;
                }

	        REAL  work1 = grad_coef*dij1;
	        REAL  work2 = grad_coef*dij2;
        	REAL  work3 = grad_coef*dij3;

    	        force_local(1) -= work1;
	        force_local(2) -= work2;
	        force_local(3) -= work3;
    	    }
        }

	// update force(:,:,i)
	WARP_RSUM_12( force_local(1) );
	WARP_RSUM_12( force_local(2) );
	WARP_RSUM_12( force_local(3) );
  	if ( id_thread_xy == 0 ) {
           if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
           if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
           if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
  	}
    }

    if (check_virial != 0 && id_thread_x < 3) {
       int n = id_thread_x + 1;
       if (n == 1) ene_viri_mid(n,index) = 0.0;
       if (n == 2) ene_viri_mid(n,index) = 0.0;
       if (n == 3) ene_viri_mid(n,index) = 0.0;
    }
}


__global__ void kern_compute_force_nonbond_table_ljpme_univ__force_intra_cell(
    const REAL  * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:ncel_max)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  * _table_grad,          // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near )
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_s,
    int  index_e,
    int  check_virial,
    REAL  density,
    REAL  cutoff2,
    REAL  pairlistdist2
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];

    int  index = index_s + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_e ) return;
    int  univ_ij = index;

    const int  i = univ_ij;
    const int  j = univ_ij;

    int  ix_natom = natom(i);
    int  iy_natom = natom(j);
    if ( ix_natom * iy_natom <= 0 ) return;

    for ( int ix = id_thread_xx + 1; ix <= ix_natom; ix += num_thread_xx ) {

	REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i));
	REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i));
	REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i));
	REAL  iqtmp   = __ldg(&charge(ix,i));
	int   iatmcls = __ldg(&atmcls(ix,i));
        REAL  lj6_i   = __ldg(&nonb_lj6_factor(iatmcls));

	force_local(1) = 0.0;
	force_local(2) = 0.0;
	force_local(3) = 0.0;

	for ( int iy = id_thread_xy + 1; iy <= iy_natom; iy += num_thread_xy ) {

            REAL  grad_coef = 0.0;
            REAL  dij1 = 0.0;
            REAL  dij2 = 0.0;
            REAL  dij3 = 0.0;
            REAL  rij2;

            int idx = iy + (ix-1)*univ_natom_max;

            if (univ_mask2(idx,univ_ij)) {

                dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

                if ( rij2 < cutoff2 && rij2 > 0.0) {

                    REAL inv_r2 = 1.0 / rij2;
                    REAL inv_r6 = inv_r2 * inv_r2 * inv_r2;
                    REAL inv_r12 = inv_r6 * inv_r6;
                    rij2 = cutoff2 * density * inv_r2;

                    REAL  jqtmp   = __ldg(&charge(iy,j));
                    int   jatmcls = __ldg(&atmcls(iy,j));
                    REAL  lj6_j   = __ldg(&nonb_lj6_factor(jatmcls));
                    REAL  lj6_ij  = lj6_i * lj6_j;
                    REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
                    REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls)) - lj6_ij;

                    int   L  = int(rij2);
                    REAL  R  = rij2 - L;
                    int   L1 = 2*L - 1;
                    REAL  tg0  = __ldg(&table_grad(L1  ));
                    REAL  tg1  = __ldg(&table_grad(L1+1));
                    REAL  tg2  = __ldg(&table_grad(L1+2));
                    REAL  tg3  = __ldg(&table_grad(L1+3));
                    REAL  term_lj6  = tg0 + R*(tg2-tg0);
                    REAL  term_elec = tg1 + R*(tg3-tg1);
                    REAL  term_lj12 = -12.0 * inv_r12 * inv_r2;
                    REAL  term_temp = - 6.0 * inv_r6  * inv_r2;
                     
                    grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij + iqtmp*jqtmp*term_elec;
                }

	        REAL  work1 = grad_coef*dij1;
	        REAL  work2 = grad_coef*dij2;
        	REAL  work3 = grad_coef*dij3;

    	        force_local(1) -= work1;
	        force_local(2) -= work2;
	        force_local(3) -= work3;
    	    }
        }

	// update force(:,:,i)
	WARP_RSUM_12( force_local(1) );
	WARP_RSUM_12( force_local(2) );
	WARP_RSUM_12( force_local(3) );
  	if ( id_thread_xy == 0 ) {
           if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
           if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
           if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
  	}
    }

    if (check_virial != 0 && id_thread_x < 3) {
       int n = id_thread_x + 1;
       if (n == 1) ene_viri_mid(n,index) = 0.0;
       if (n == 2) ene_viri_mid(n,index) = 0.0;
       if (n == 3) ene_viri_mid(n,index) = 0.0;
    }
}
/*
 *
 */
#if defined(_MIXED) || defined(_SINGLE)
//__launch_bounds__(128,12)  // for SP
#else
//__launch_bounds__(128,9)  // for DP
#endif
__global__ void kern_compute_force_nonbond_notable_univ__force_intra_cell(
    const REAL  * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:ncel_max)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _table_grad,          // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near )
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_s,
    int  index_e,
    int  check_virial,
    REAL  density,
    REAL  cutoff2,
    REAL  pairlistdist2
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];

    int  index = index_s + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_e ) return;
    int  univ_ij = index;

    const int  i = univ_ij;
    const int  j = univ_ij;

    int  ix_natom = natom(i);
    int  iy_natom = natom(j);
    if ( ix_natom * iy_natom <= 0 ) return;

    for ( int ix = id_thread_xx + 1; ix <= ix_natom; ix += num_thread_xx ) {

	REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i));
	REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i));
	REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i));
	REAL  iqtmp   = __ldg(&charge(ix,i));
	int   iatmcls = __ldg(&atmcls(ix,i));

	force_local(1) = 0.0;
	force_local(2) = 0.0;
	force_local(3) = 0.0;

	for ( int iy = id_thread_xy + 1; iy <= iy_natom; iy += num_thread_xy ) {

            REAL  grad_coef = 0.0;
            REAL  dij1 = 0.0;
            REAL  dij2 = 0.0;
            REAL  dij3 = 0.0;
            REAL  rij2;

            int idx = iy + (ix-1)*univ_natom_max;

            if (univ_mask2(idx,univ_ij)) {

                dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

                if ( rij2 < cutoff2) {

                    REAL rij2_inv = 1.0 / rij2;

                    rij2 = cutoff2 * density * rij2_inv;

                    REAL  jqtmp   = __ldg(&charge(iy,j));
                    int   jatmcls = __ldg(&atmcls(iy,j));
                    REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
                    REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));

                    int   L  = int(rij2);
                    REAL  R  = rij2 - L;
                    REAL  tg0  = __ldg(&table_grad(L  ));
                    REAL  tg1  = __ldg(&table_grad(L+1));
                    REAL  term_elec = tg0 + R*(tg1-tg0);
                    REAL  term_lj6  = rij2_inv * rij2_inv * rij2_inv;
                    REAL  term_lj12 = term_lj6 * term_lj6;
                    term_lj12 = -12.0 * term_lj12 * rij2_inv;
                    term_lj6  = - 6.0 * term_lj6  * rij2_inv;

                    grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;
                }

	        REAL  work1 = grad_coef*dij1;
	        REAL  work2 = grad_coef*dij2;
        	REAL  work3 = grad_coef*dij3;

    	        force_local(1) -= work1;
	        force_local(2) -= work2;
	        force_local(3) -= work3;
    	    }
        }

	// update force(:,:,i)
	WARP_RSUM_12( force_local(1) );
	WARP_RSUM_12( force_local(2) );
	WARP_RSUM_12( force_local(3) );
  	if ( id_thread_xy == 0 ) {
           if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
           if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
           if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
  	}
    }

    if (check_virial != 0 && id_thread_x < 3) {
       int n = id_thread_x + 1;
       if (n == 1) ene_viri_mid(n,index) = 0.0;
       if (n == 2) ene_viri_mid(n,index) = 0.0;
       if (n == 3) ene_viri_mid(n,index) = 0.0;
    }
}


/* */

#if defined(_MIXED) || defined(_SINGLE)
#define NUM_CTA__ENERGYFORCE_INTER_CELL 12
#else
#define NUM_CTA__ENERGYFORCE_INTER_CELL  8
#endif
/* */
//__launch_bounds__(128,NUM_CTA__ENERGYFORCE_INTER_CELL)
__global__ void kern_compute_energy_nonbond_table_linear_univ__energyforce_inter_cell(
    const REAL  * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:univ_maxcell)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _table_ene,           // ( 1:6*cutoff_int )
    const REAL  * _table_grad,          // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  * _virial_check,        // ( 1:ncel_max, 1:ncel_max )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_ncell_near,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_start,
    int  index_end,
    int  max_iy_natom,
    REAL  density,
    REAL  cutoff2,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

    const int  num_thread_x = blockDim.x;
    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];

    /* shared memory */
    REAL  *_force_iy_smem = & smem[id_thread_y * 3*max_iy_natom]; // 3 * iy_natom
#define force_iy_smem(X,Y) _force_iy_smem[CALIDX2((X)-1,3, (Y)-1,iy_natom)]

    int  index = index_start + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_end ) return;
    int  univ_ij = univ_ij_sort_list( index );

    int  ix_natom = univ_ix_natom(univ_ij);
    int  iy_natom = univ_iy_natom(univ_ij);
    if ( ix_natom * iy_natom <= 0 ) return;

    const int  i = univ_cell_pairlist1(1,univ_ij);
    const int  j = univ_cell_pairlist1(2,univ_ij);
    int  k;

#define sumval(Z) _sumval[(Z)-1]
    double _sumval[5];
    sumval(1) = 0.0;  // elec
    sumval(2) = 0.0;  // evdw
    sumval(3) = 0.0;  // virial(1)
    sumval(4) = 0.0;  // virial(2)
    sumval(5) = 0.0;  // virial(3)
    const int  check_virial = virial_check(j,i);

    if (univ_ij > univ_ncell_near){

        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
            int  iiy_e = iiy_s + max_iy_natom - 1;
            if ( iiy_e > iy_natom ) iiy_e = iy_natom;

            // initialize force_iy at shared memory
            for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
                int   n   = (k % 3) + 1;
                int   iiy = (k / 3) + iiy_s;
                force_iy_smem(n,iiy-iiy_s+1) = 0.0;
            }

            for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
                int  ix = __ldg(&univ_ix_list(iix,univ_ij));

                REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
                REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
                REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
                REAL  iqtmp   = __ldg(&charge(ix,i));
                int   iatmcls = __ldg(&atmcls(ix,i));

                force_local(1) = 0.0;
                force_local(2) = 0.0;
                force_local(3) = 0.0;
                REAL  elec_temp = 0.0;
                REAL  evdw_temp = 0.0;

                for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
                    int  iy = __ldg(&univ_iy_list(iiy,univ_ij));

                    REAL  grad_coef = 0.0;
                    REAL  dij1 = 0.0;
                    REAL  dij2 = 0.0;
                    REAL  dij3 = 0.0;
                    REAL  rij2;

                    dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                    dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                    dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                    rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;
                    grad_coef = 0.0;

                    if ( rij2 < cutoff2 ) {

                        rij2 = cutoff2 * density / rij2;

                        REAL  jqtmp   = __ldg(&charge(iy,j));
                        int   jatmcls = __ldg(&atmcls(iy,j));
                        REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
                        REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));

                        int   L  = int(rij2);
                        REAL  R  = rij2 - L;
                        int   L1 = 3*L - 2;

                        REAL  tg0  = __ldg(&table_grad(L1  ));
                        REAL  tg1  = __ldg(&table_grad(L1+1));
                        REAL  tg2  = __ldg(&table_grad(L1+2));
                        REAL  tg3  = __ldg(&table_grad(L1+3));
                        REAL  tg4  = __ldg(&table_grad(L1+4));
                        REAL  tg5  = __ldg(&table_grad(L1+5));
                        REAL  term_lj12 = tg0 + R*(tg3-tg0);
                        REAL  term_lj6  = tg1 + R*(tg4-tg1);
                        REAL  term_elec = tg2 + R*(tg5-tg2);

                        grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;

                        /* energy */
                        tg0  = __ldg(&table_ene(L1  ));
                        tg1  = __ldg(&table_ene(L1+1));
                        tg2  = __ldg(&table_ene(L1+2));
                        tg3  = __ldg(&table_ene(L1+3));
                        tg4  = __ldg(&table_ene(L1+4));
                        tg5  = __ldg(&table_ene(L1+5));
                        term_lj12 = tg0 + R*(tg3-tg0);
                        term_lj6  = tg1 + R*(tg4-tg1);
                        term_elec = tg2 + R*(tg5-tg2);
                        evdw_temp += term_lj12*lj12 - term_lj6*lj6;
                        elec_temp += iqtmp*jqtmp*term_elec;
                    }

                    REAL  work1 = grad_coef*dij1;
                    REAL  work2 = grad_coef*dij2;
                    REAL  work3 = grad_coef*dij3;

                    force_local(1) -= work1;
                    force_local(2) -= work2;
                    force_local(3) -= work3;

                    // update force_iy(:,iiy) at smem
                    WARP_RSUM_345( work1 );
                    WARP_RSUM_345( work2 );
                    WARP_RSUM_345( work3 );
                    if ( id_thread_xx == 0 ) {
                        force_iy_smem(1,iiy-iiy_s+1) += work1;
                        force_iy_smem(2,iiy-iiy_s+1) += work2;
                        force_iy_smem(3,iiy-iiy_s+1) += work3;
                    }
                }
                // energy and virial
                sumval(1) += elec_temp;
                sumval(2) += evdw_temp;
                if (check_virial != 0) {
                    sumval(3) += force_local(1);  // virial(1)
                    sumval(4) += force_local(2);  // virial(2)
                    sumval(5) += force_local(3);  // virial(3)
                }

                // update force(:,:,i)
                WARP_RSUM_12( force_local(1) );
                WARP_RSUM_12( force_local(2) );
                WARP_RSUM_12( force_local(3) );
                if ( id_thread_xy == 0 ) {
                    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
                    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
                    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
                }
            }

            // __syncthreads();

            // update force(:,:,j)
            for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
                int   n   = (k % 3) + 1;
                int   iiy = (k / 3) + iiy_s;
                int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
                REAL  val = force_iy_smem(n,iiy-iiy_s+1);
                if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
            }
        }
    }

    else{

        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
  	    int  iiy_e = iiy_s + max_iy_natom - 1;
    	    if ( iiy_e > iy_natom ) iiy_e = iy_natom;

	// initialize force_iy at shared memory
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int   n   = (k % 3) + 1;
	        int   iiy = (k / 3) + iiy_s;
	        force_iy_smem(n,iiy-iiy_s+1) = 0.0;
	    }

	    for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
	        int  ix = __ldg(&univ_ix_list(iix,univ_ij));

    	        REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
	        REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
	        REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
	        REAL  iqtmp   = __ldg(&charge(ix,i));
	        int   iatmcls = __ldg(&atmcls(ix,i));

	        force_local(1) = 0.0;
	        force_local(2) = 0.0;
	        force_local(3) = 0.0;
	        REAL  elec_temp = 0.0;
	        REAL  evdw_temp = 0.0;

	        for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
		    int  iy = __ldg(&univ_iy_list(iiy,univ_ij));

                    REAL  grad_coef = 0.0;
                    REAL  dij1 = 0.0;
                    REAL  dij2 = 0.0;
                    REAL  dij3 = 0.0;
                    REAL  rij2;

                    int idx = iy + (ix-1)*univ_natom_max;
                    if (univ_mask2(idx,univ_ij)) {

                        dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                        dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                        dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;
                        grad_coef = 0.0;

                        if ( rij2 < cutoff2 ) {

                            rij2 = cutoff2 * density / rij2;

    		            REAL  jqtmp   = __ldg(&charge(iy,j));
		            int   jatmcls = __ldg(&atmcls(iy,j));
		            REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
		            REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));

		            int   L  = int(rij2);
		            REAL  R  = rij2 - L;
		            int   L1 = 3*L - 2;

		            REAL  tg0  = __ldg(&table_grad(L1  ));
		            REAL  tg1  = __ldg(&table_grad(L1+1));
		            REAL  tg2  = __ldg(&table_grad(L1+2));
		            REAL  tg3  = __ldg(&table_grad(L1+3));
		            REAL  tg4  = __ldg(&table_grad(L1+4));
		            REAL  tg5  = __ldg(&table_grad(L1+5));
		            REAL  term_lj12 = tg0 + R*(tg3-tg0);
		            REAL  term_lj6  = tg1 + R*(tg4-tg1);
		            REAL  term_elec = tg2 + R*(tg5-tg2);

		            grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;

		            /* energy */
		            tg0  = __ldg(&table_ene(L1  ));
		            tg1  = __ldg(&table_ene(L1+1));
		            tg2  = __ldg(&table_ene(L1+2));
		            tg3  = __ldg(&table_ene(L1+3));
		            tg4  = __ldg(&table_ene(L1+4));
		            tg5  = __ldg(&table_ene(L1+5));
		            term_lj12 = tg0 + R*(tg3-tg0);
		            term_lj6  = tg1 + R*(tg4-tg1);
		            term_elec = tg2 + R*(tg5-tg2);
                            evdw_temp += term_lj12*lj12 - term_lj6*lj6;
                            elec_temp += iqtmp*jqtmp*term_elec;
		        }
                    }

		    REAL  work1 = grad_coef*dij1;
		    REAL  work2 = grad_coef*dij2;
		    REAL  work3 = grad_coef*dij3;

		    force_local(1) -= work1;
		    force_local(2) -= work2;
		    force_local(3) -= work3;

		    WARP_RSUM_345( work1 );
		    WARP_RSUM_345( work2 );
		    WARP_RSUM_345( work3 );
		    if ( id_thread_xx == 0 ) {
		        force_iy_smem(1,iiy-iiy_s+1) += work1;
		        force_iy_smem(2,iiy-iiy_s+1) += work2;
		        force_iy_smem(3,iiy-iiy_s+1) += work3;
		    }
	        }

	        // energy and virial
	        sumval(1) += elec_temp;
	        sumval(2) += evdw_temp;
	        if (check_virial != 0) {
		    sumval(3) += force_local(1);  // virial(1)
		    sumval(4) += force_local(2);  // virial(2)
		    sumval(5) += force_local(3);  // virial(3)
	        }

	        // update force(:,:,i)
	        WARP_RSUM_12( force_local(1) );
	        WARP_RSUM_12( force_local(2) );
	        WARP_RSUM_12( force_local(3) );
	        if ( id_thread_xy == 0 ) {
		    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
		    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
		    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
	        }
	    }

	    // __syncthreads();

	    // update force(:,:,j)
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int   n   = (k % 3) + 1;
	        int   iiy = (k / 3) + iiy_s;
	        int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
	        REAL  val = force_iy_smem(n,iiy-iiy_s+1);
	        if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
	    }
        }
    }

    // update energy and virial
    sumval(3) *= __ldg(&cell_move(1,j,i))*system_x;  // virial(1)
    sumval(4) *= __ldg(&cell_move(2,j,i))*system_y;  // virial(2)
    sumval(5) *= __ldg(&cell_move(3,j,i))*system_z;  // virial(3)

    WARP_RSUM_12345( sumval(1) );  // elec
    WARP_RSUM_12345( sumval(2) );  // evdw
    WARP_RSUM_12345( sumval(3) );  // virial(1)
    WARP_RSUM_12345( sumval(4) );  // virial(2)
    WARP_RSUM_12345( sumval(5) );  // virial(3)
    if (id_thread_x < 5) {
        int n = id_thread_x + 1;
	ene_viri_mid(n,index) = sumval(n);
    }
}


__global__ void kern_compute_energy_nonbond_table_ljpme_univ__energyforce_inter_cell(
    const REAL  * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:univ_maxcell)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6_factor,     // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _table_ene,           // ( 1:6*cutoff_int )
    const REAL  * _table_grad,          // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  * _virial_check,        // ( 1:ncel_max, 1:ncel_max )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_ncell_near,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_start,
    int  index_end,
    int  max_iy_natom,
    REAL  density,
    REAL  cutoff2,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

    const int  num_thread_x = blockDim.x;
    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];

    /* shared memory */
    REAL  *_force_iy_smem = & smem[id_thread_y * 3*max_iy_natom]; // 3 * iy_natom
#define force_iy_smem(X,Y) _force_iy_smem[CALIDX2((X)-1,3, (Y)-1,iy_natom)]

    int  index = index_start + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_end ) return;
    int  univ_ij = univ_ij_sort_list( index );

    int  ix_natom = univ_ix_natom(univ_ij);
    int  iy_natom = univ_iy_natom(univ_ij);
    if ( ix_natom * iy_natom <= 0 ) return;

    const int  i = univ_cell_pairlist1(1,univ_ij);
    const int  j = univ_cell_pairlist1(2,univ_ij);
    int  k;

    REAL inv_cutoff2 = 1.0 / cutoff2;
    REAL inv_cutoff6 = inv_cutoff2 * inv_cutoff2 * inv_cutoff2;

#define sumval(Z) _sumval[(Z)-1]
    double _sumval[5];
    sumval(1) = 0.0;  // elec
    sumval(2) = 0.0;  // evdw
    sumval(3) = 0.0;  // virial(1)
    sumval(4) = 0.0;  // virial(2)
    sumval(5) = 0.0;  // virial(3)
    const int  check_virial = virial_check(j,i);

    if (univ_ij > univ_ncell_near){

        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
            int  iiy_e = iiy_s + max_iy_natom - 1;
            if ( iiy_e > iy_natom ) iiy_e = iy_natom;

            // initialize force_iy at shared memory
            for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
                int   n   = (k % 3) + 1;
                int   iiy = (k / 3) + iiy_s;
                force_iy_smem(n,iiy-iiy_s+1) = 0.0;
            }

            for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
                int  ix = __ldg(&univ_ix_list(iix,univ_ij));

                REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
                REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
                REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
                REAL  iqtmp   = __ldg(&charge(ix,i));
                int   iatmcls = __ldg(&atmcls(ix,i));
                REAL  lj6_i   = __ldg(&nonb_lj6_factor(iatmcls));

                force_local(1) = 0.0;
                force_local(2) = 0.0;
                force_local(3) = 0.0;
                REAL  elec_temp = 0.0;
                REAL  evdw_temp = 0.0;

                for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
                    int  iy = __ldg(&univ_iy_list(iiy,univ_ij));

                    REAL  grad_coef = 0.0;
                    REAL  dij1 = 0.0;
                    REAL  dij2 = 0.0;
                    REAL  dij3 = 0.0;
                    REAL  rij2;

                    dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                    dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                    dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                    rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;
                    grad_coef = 0.0;

                    if ( rij2 < cutoff2 ) {

                        REAL inv_r2 = 1.0 / rij2;
                        REAL inv_r6 = inv_r2 * inv_r2 * inv_r2;
                        REAL term_lj12 = inv_r6 * inv_r6;
                        rij2 = cutoff2 * density * inv_r2;

                        REAL  jqtmp   = __ldg(&charge(iy,j));
                        int   jatmcls = __ldg(&atmcls(iy,j));
                        REAL  lj6_j   = __ldg(&nonb_lj6_factor(jatmcls));
                        REAL  lj6_ij  = lj6_i * lj6_j;
                        REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
                        REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls)) - lj6_ij;

                        int   L  = int(rij2);
                        REAL  R  = rij2 - L;
                        int   L1 = 2*L - 1;

                        /* energy */
                        REAL tg0  = __ldg(&table_ene(L1  ));
                        REAL tg1  = __ldg(&table_ene(L1+1));
                        REAL tg2  = __ldg(&table_ene(L1+2));
                        REAL tg3  = __ldg(&table_ene(L1+3));
                        REAL term_lj6  = tg0 + R*(tg2-tg0);
                        REAL term_elec = tg1 + R*(tg3-tg1);
                        REAL term_temp = inv_r6 - inv_cutoff6;

                        evdw_temp += term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij;
                        elec_temp += iqtmp*jqtmp*term_elec;

                        tg0  = __ldg(&table_grad(L1  ));
                        tg1  = __ldg(&table_grad(L1+1));
                        tg2  = __ldg(&table_grad(L1+2));
                        tg3  = __ldg(&table_grad(L1+3));
                        term_lj6  = tg0 + R*(tg2-tg0);
                        term_elec = tg1 + R*(tg3-tg1);
                        term_lj12 = -12.0 * term_lj12 * inv_r2;
                        term_temp = - 6.0 * inv_r6 * inv_r2;
                        
                        grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij + iqtmp*jqtmp*term_elec;

                    }

                    REAL  work1 = grad_coef*dij1;
                    REAL  work2 = grad_coef*dij2;
                    REAL  work3 = grad_coef*dij3;

                    force_local(1) -= work1;
                    force_local(2) -= work2;
                    force_local(3) -= work3;

                    // update force_iy(:,iiy) at smem
                    WARP_RSUM_345( work1 );
                    WARP_RSUM_345( work2 );
                    WARP_RSUM_345( work3 );
                    if ( id_thread_xx == 0 ) {
                        force_iy_smem(1,iiy-iiy_s+1) += work1;
                        force_iy_smem(2,iiy-iiy_s+1) += work2;
                        force_iy_smem(3,iiy-iiy_s+1) += work3;
                    }
                }
                // energy and virial
                sumval(1) += elec_temp;
                sumval(2) += evdw_temp;
                if (check_virial != 0) {
                    sumval(3) += force_local(1);  // virial(1)
                    sumval(4) += force_local(2);  // virial(2)
                    sumval(5) += force_local(3);  // virial(3)
                }

                // update force(:,:,i)
                WARP_RSUM_12( force_local(1) );
                WARP_RSUM_12( force_local(2) );
                WARP_RSUM_12( force_local(3) );
                if ( id_thread_xy == 0 ) {
                    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
                    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
                    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
                }
            }

            // __syncthreads();

            // update force(:,:,j)
            for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
                int   n   = (k % 3) + 1;
                int   iiy = (k / 3) + iiy_s;
                int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
                REAL  val = force_iy_smem(n,iiy-iiy_s+1);
                if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
            }
        }
    }

    else{

        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
  	    int  iiy_e = iiy_s + max_iy_natom - 1;
    	    if ( iiy_e > iy_natom ) iiy_e = iy_natom;

	// initialize force_iy at shared memory
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int   n   = (k % 3) + 1;
	        int   iiy = (k / 3) + iiy_s;
	        force_iy_smem(n,iiy-iiy_s+1) = 0.0;
	    }

	    for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
	        int  ix = __ldg(&univ_ix_list(iix,univ_ij));

    	        REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
	        REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
	        REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
	        REAL  iqtmp   = __ldg(&charge(ix,i));
	        int   iatmcls = __ldg(&atmcls(ix,i));
                REAL  lj6_i   = __ldg(&nonb_lj6_factor(iatmcls));

	        force_local(1) = 0.0;
	        force_local(2) = 0.0;
	        force_local(3) = 0.0;
	        REAL  elec_temp = 0.0;
	        REAL  evdw_temp = 0.0;

	        for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
		    int  iy = __ldg(&univ_iy_list(iiy,univ_ij));

                    REAL  grad_coef = 0.0;
                    REAL  dij1 = 0.0;
                    REAL  dij2 = 0.0;
                    REAL  dij3 = 0.0;
                    REAL  rij2;

                    int idx = iy + (ix-1)*univ_natom_max;
                    if (univ_mask2(idx,univ_ij)) {

                        dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                        dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                        dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;
                        grad_coef = 0.0;

                        if ( rij2 < cutoff2 ) {

                            REAL  inv_r2 = 1.0 / rij2;
                            REAL  inv_r6 = inv_r2 * inv_r2 * inv_r2;
                            REAL  term_lj12 = inv_r6 * inv_r6;
                            rij2 = cutoff2 * density * inv_r2;

    		            REAL  jqtmp   = __ldg(&charge(iy,j));
		            int   jatmcls = __ldg(&atmcls(iy,j));
                            REAL  lj6_j   = __ldg(&nonb_lj6_factor(jatmcls));
                            REAL  lj6_ij  = lj6_i * lj6_j;
		            REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
		            REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls)) - lj6_ij;

		            int   L  = int(rij2);
		            REAL  R  = rij2 - L;
		            int   L1 = 2*L - 1;

                            /* energy */
                            REAL  tg0  = __ldg(&table_ene(L1  ));
                            REAL  tg1  = __ldg(&table_ene(L1+1));
                            REAL  tg2  = __ldg(&table_ene(L1+2));
                            REAL  tg3  = __ldg(&table_ene(L1+3));
                            REAL  term_lj6  = tg0 + R*(tg2-tg0);
                            REAL  term_elec = tg1 + R*(tg3-tg1);
                            REAL  term_temp = inv_r6 - inv_cutoff6;
                            evdw_temp += term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij;
                            elec_temp += iqtmp*jqtmp*term_elec;

		            tg0  = __ldg(&table_grad(L1  ));
		            tg1  = __ldg(&table_grad(L1+1));
		            tg2  = __ldg(&table_grad(L1+2));
		            tg3  = __ldg(&table_grad(L1+3));
		            term_lj6  = tg0 + R*(tg2-tg0);
		            term_elec = tg1 + R*(tg3-tg1);
		            term_lj12 = -12.0 * term_lj12 * inv_r2;
                            term_temp = - 6.0 * inv_r6 * inv_r2;

		            grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij + iqtmp*jqtmp*term_elec;
		        }
                    }

		    REAL  work1 = grad_coef*dij1;
		    REAL  work2 = grad_coef*dij2;
		    REAL  work3 = grad_coef*dij3;

		    force_local(1) -= work1;
		    force_local(2) -= work2;
		    force_local(3) -= work3;

		    WARP_RSUM_345( work1 );
		    WARP_RSUM_345( work2 );
		    WARP_RSUM_345( work3 );
		    if ( id_thread_xx == 0 ) {
		        force_iy_smem(1,iiy-iiy_s+1) += work1;
		        force_iy_smem(2,iiy-iiy_s+1) += work2;
		        force_iy_smem(3,iiy-iiy_s+1) += work3;
		    }
	        }

	        // energy and virial
	        sumval(1) += elec_temp;
	        sumval(2) += evdw_temp;
	        if (check_virial != 0) {
		    sumval(3) += force_local(1);  // virial(1)
		    sumval(4) += force_local(2);  // virial(2)
		    sumval(5) += force_local(3);  // virial(3)
	        }

	        // update force(:,:,i)
	        WARP_RSUM_12( force_local(1) );
	        WARP_RSUM_12( force_local(2) );
	        WARP_RSUM_12( force_local(3) );
	        if ( id_thread_xy == 0 ) {
		    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
		    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
		    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
	        }
	    }

	    // __syncthreads();

	    // update force(:,:,j)
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int   n   = (k % 3) + 1;
	        int   iiy = (k / 3) + iiy_s;
	        int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
	        REAL  val = force_iy_smem(n,iiy-iiy_s+1);
	        if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
	    }
        }
    }

    // update energy and virial
    sumval(3) *= __ldg(&cell_move(1,j,i))*system_x;  // virial(1)
    sumval(4) *= __ldg(&cell_move(2,j,i))*system_y;  // virial(2)
    sumval(5) *= __ldg(&cell_move(3,j,i))*system_z;  // virial(3)

    WARP_RSUM_12345( sumval(1) );  // elec
    WARP_RSUM_12345( sumval(2) );  // evdw
    WARP_RSUM_12345( sumval(3) );  // virial(1)
    WARP_RSUM_12345( sumval(4) );  // virial(2)
    WARP_RSUM_12345( sumval(5) );  // virial(3)
    if (id_thread_x < 5) {
        int n = id_thread_x + 1;
	ene_viri_mid(n,index) = sumval(n);
    }
}
/* */

#if defined(_MIXED) || defined(_SINGLE)
#define NUM_CTA__ENERGYFORCE_INTER_CELL 12
#else
#define NUM_CTA__ENERGYFORCE_INTER_CELL  8
#endif
/* */
//__launch_bounds__(128,NUM_CTA__ENERGYFORCE_INTER_CELL)
__global__ void kern_compute_energy_nonbond_notable_univ__energyforce_inter_cell(
    const REAL  * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:univ_maxcell)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _table_ene,           // ( 1:6*cutoff_int )
    const REAL  * _table_grad,          // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  * _virial_check,        // ( 1:ncel_max, 1:ncel_max )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_ncell_near,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_start,
    int  index_end,
    int  max_iy_natom,
    REAL  density,
    REAL  cutoff2,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

    const int  num_thread_x = blockDim.x;
    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];

    /* shared memory */
    REAL  *_force_iy_smem = & smem[id_thread_y * 3*max_iy_natom]; // 3 * iy_natom
#define force_iy_smem(X,Y) _force_iy_smem[CALIDX2((X)-1,3, (Y)-1,iy_natom)]

    int  index = index_start + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_end ) return;
    int  univ_ij = univ_ij_sort_list( index );

    int  ix_natom = univ_ix_natom(univ_ij);
    int  iy_natom = univ_iy_natom(univ_ij);
    if ( ix_natom * iy_natom <= 0 ) return;

    const int  i = univ_cell_pairlist1(1,univ_ij);
    const int  j = univ_cell_pairlist1(2,univ_ij);
    int  k;

#define sumval(Z) _sumval[(Z)-1]
    double _sumval[5];
    sumval(1) = 0.0;  // elec
    sumval(2) = 0.0;  // evdw
    sumval(3) = 0.0;  // virial(1)
    sumval(4) = 0.0;  // virial(2)
    sumval(5) = 0.0;  // virial(3)
    const int  check_virial = virial_check(j,i);

    if (univ_ij > univ_ncell_near){

        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
            int  iiy_e = iiy_s + max_iy_natom - 1;
            if ( iiy_e > iy_natom ) iiy_e = iy_natom;

            // initialize force_iy at shared memory
            for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
                int   n   = (k % 3) + 1;
                int   iiy = (k / 3) + iiy_s;
                force_iy_smem(n,iiy-iiy_s+1) = 0.0;
            }

            for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
                int  ix = __ldg(&univ_ix_list(iix,univ_ij));

                REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
                REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
                REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
                REAL  iqtmp   = __ldg(&charge(ix,i));
                int   iatmcls = __ldg(&atmcls(ix,i));

                force_local(1) = 0.0;
                force_local(2) = 0.0;
                force_local(3) = 0.0;
                REAL  elec_temp = 0.0;
                REAL  evdw_temp = 0.0;

                for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
                    int  iy = __ldg(&univ_iy_list(iiy,univ_ij));

                    REAL  grad_coef = 0.0;
                    REAL  dij1 = 0.0;
                    REAL  dij2 = 0.0;
                    REAL  dij3 = 0.0;
                    REAL  rij2;

                    dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                    dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                    dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                    rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;
                    grad_coef = 0.0;

                    if ( rij2 < cutoff2 ) {

                        REAL rij2_inv = 1.0 / rij2;
                        rij2 = cutoff2 * density * rij2_inv;

                        REAL  jqtmp   = __ldg(&charge(iy,j));
                        int   jatmcls = __ldg(&atmcls(iy,j));
                        REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
                        REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));

                        int   L  = int(rij2);
                        REAL  R  = rij2 - L;

                        /* energy */
                        REAL tg0  = __ldg(&table_ene(L  ));
                        REAL tg1  = __ldg(&table_ene(L+1));
                        REAL term_elec = tg0 + R*(tg1-tg0);
                        REAL term_lj6  = rij2_inv * rij2_inv * rij2_inv;
                        REAL term_lj12 = term_lj6 * term_lj6;
                        evdw_temp += term_lj12*lj12 - term_lj6*lj6;
                        elec_temp += iqtmp*jqtmp*term_elec;
 
                        tg0  = __ldg(&table_grad(L  ));
                        tg1  = __ldg(&table_grad(L+1));
                        term_elec = tg0 + R*(tg1-tg0);
                        term_lj12 = -12.0 * term_lj12 * rij2_inv;
                        term_lj6  = - 6.0 * term_lj6  * rij2_inv;

                        grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;
                    }

                    REAL  work1 = grad_coef*dij1;
                    REAL  work2 = grad_coef*dij2;
                    REAL  work3 = grad_coef*dij3;

                    force_local(1) -= work1;
                    force_local(2) -= work2;
                    force_local(3) -= work3;

                    // update force_iy(:,iiy) at smem
                    WARP_RSUM_345( work1 );
                    WARP_RSUM_345( work2 );
                    WARP_RSUM_345( work3 );
                    if ( id_thread_xx == 0 ) {
                        force_iy_smem(1,iiy-iiy_s+1) += work1;
                        force_iy_smem(2,iiy-iiy_s+1) += work2;
                        force_iy_smem(3,iiy-iiy_s+1) += work3;
                    }
                }
                // energy and virial
                sumval(1) += elec_temp;
                sumval(2) += evdw_temp;
                if (check_virial != 0) {
                    sumval(3) += force_local(1);  // virial(1)
                    sumval(4) += force_local(2);  // virial(2)
                    sumval(5) += force_local(3);  // virial(3)
                }

                // update force(:,:,i)
                WARP_RSUM_12( force_local(1) );
                WARP_RSUM_12( force_local(2) );
                WARP_RSUM_12( force_local(3) );
                if ( id_thread_xy == 0 ) {
                    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
                    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
                    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
                }
            }

            // __syncthreads();

            // update force(:,:,j)
            for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
                int   n   = (k % 3) + 1;
                int   iiy = (k / 3) + iiy_s;
                int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
                REAL  val = force_iy_smem(n,iiy-iiy_s+1);
                if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
            }
        }
    }

    else{

        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
  	    int  iiy_e = iiy_s + max_iy_natom - 1;
    	    if ( iiy_e > iy_natom ) iiy_e = iy_natom;

	// initialize force_iy at shared memory
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int   n   = (k % 3) + 1;
	        int   iiy = (k / 3) + iiy_s;
	        force_iy_smem(n,iiy-iiy_s+1) = 0.0;
	    }

	    for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
	        int  ix = __ldg(&univ_ix_list(iix,univ_ij));

    	        REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
	        REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
	        REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
	        REAL  iqtmp   = __ldg(&charge(ix,i));
	        int   iatmcls = __ldg(&atmcls(ix,i));

	        force_local(1) = 0.0;
	        force_local(2) = 0.0;
	        force_local(3) = 0.0;
	        REAL  elec_temp = 0.0;
	        REAL  evdw_temp = 0.0;

	        for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
		    int  iy = __ldg(&univ_iy_list(iiy,univ_ij));

                    REAL  grad_coef = 0.0;
                    REAL  dij1 = 0.0;
                    REAL  dij2 = 0.0;
                    REAL  dij3 = 0.0;
                    REAL  rij2;

                    int idx = iy + (ix-1)*univ_natom_max;
                    if (univ_mask2(idx,univ_ij)) {

                        dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
                        dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
                        dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
                        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;
                        grad_coef = 0.0;

                        if ( rij2 < cutoff2 ) {

                            REAL rij2_inv = 1.0 / rij2;
                            rij2 = cutoff2 * density * rij2_inv;

                            REAL  jqtmp   = __ldg(&charge(iy,j));
                            int   jatmcls = __ldg(&atmcls(iy,j));
                            REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
                            REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));

                            int   L  = int(rij2);
                            REAL  R  = rij2 - L;

                            /* energy */
                            REAL tg0  = __ldg(&table_ene(L  ));
                            REAL tg1  = __ldg(&table_ene(L+1));
                            REAL term_elec = tg0 + R*(tg1-tg0);
                            REAL term_lj6  = rij2_inv * rij2_inv * rij2_inv;
                            REAL term_lj12 = term_lj6 * term_lj6;
                            evdw_temp += term_lj12*lj12 - term_lj6*lj6;
                            elec_temp += iqtmp*jqtmp*term_elec;

                            tg0  = __ldg(&table_grad(L  ));
                            tg1  = __ldg(&table_grad(L+1));
                            term_elec = tg0 + R*(tg1-tg0);
                            term_lj12 = -12.0 * term_lj12 * rij2_inv;
                            term_lj6  = - 6.0 * term_lj6  * rij2_inv;
      
                            grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;
		        }
                    }

		    REAL  work1 = grad_coef*dij1;
		    REAL  work2 = grad_coef*dij2;
		    REAL  work3 = grad_coef*dij3;

		    force_local(1) -= work1;
		    force_local(2) -= work2;
		    force_local(3) -= work3;

		    WARP_RSUM_345( work1 );
		    WARP_RSUM_345( work2 );
		    WARP_RSUM_345( work3 );
		    if ( id_thread_xx == 0 ) {
		        force_iy_smem(1,iiy-iiy_s+1) += work1;
		        force_iy_smem(2,iiy-iiy_s+1) += work2;
		        force_iy_smem(3,iiy-iiy_s+1) += work3;
		    }
	        }

	        // energy and virial
	        sumval(1) += elec_temp;
	        sumval(2) += evdw_temp;
	        if (check_virial != 0) {
		    sumval(3) += force_local(1);  // virial(1)
		    sumval(4) += force_local(2);  // virial(2)
		    sumval(5) += force_local(3);  // virial(3)
	        }

                // update force(:,:,i)
                WARP_RSUM_12( force_local(1) );
                WARP_RSUM_12( force_local(2) );
                WARP_RSUM_12( force_local(3) );
                if ( id_thread_xy == 0 ) {
                    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
                    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
                    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
                }
            }

            // __syncthreads();

            // update force(:,:,j)
            for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
                int   n   = (k % 3) + 1;
                int   iiy = (k / 3) + iiy_s;
                int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
                REAL  val = force_iy_smem(n,iiy-iiy_s+1);
                if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
            }
        }
    }

    // update energy and virial
    sumval(3) *= __ldg(&cell_move(1,j,i))*system_x;  // virial(1)
    sumval(4) *= __ldg(&cell_move(2,j,i))*system_y;  // virial(2)
    sumval(5) *= __ldg(&cell_move(3,j,i))*system_z;  // virial(3)

    WARP_RSUM_12345( sumval(1) );  // elec
    WARP_RSUM_12345( sumval(2) );  // evdw
    WARP_RSUM_12345( sumval(3) );  // virial(1)
    WARP_RSUM_12345( sumval(4) );  // virial(2)
    WARP_RSUM_12345( sumval(5) );  // virial(3)
    if (id_thread_x < 5) {
        int n = id_thread_x + 1;
        ene_viri_mid(n,index) = sumval(n);
    }
}



__global__ void kern_compute_energy_nonbond_table_linear_univ_energy_sum(
    double       *_ene_virial,
    const double *_ene_viri_mid,
    int          ncel_local,
    int          ncel_max,
    int          univ_maxcell,
    int          univ_ncell_nonzero
    )
{
    const int  num_thread_x = blockDim.x;
    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;

    const int num_thread   = ( blockDim.x * blockDim.y );
    const int id_thread    = ( threadIdx.x + blockDim.x * threadIdx.y );

    __shared__ double _ene_virial_smem[5*32];
#define ene_virial_smem(Y,Z)  _ene_virial_smem[CALIDX2((Y)-1,5, (Z)-1,32)]

    double energy_elec = 0.0;
    double energy_evdw = 0.0;
    double _virial[3];
    virial(1) = 0.0;
    virial(2) = 0.0;
    virial(3) = 0.0;

    int ij;
    for ( ij = id_thread+1 ; ij <= univ_ncell_nonzero ; ij += num_thread ) {
       energy_elec += ene_viri_mid(1,ij);
       energy_evdw += ene_viri_mid(2,ij);
       virial(1)   += ene_viri_mid(3,ij);
       virial(2)   += ene_viri_mid(4,ij);
       virial(3)   += ene_viri_mid(5,ij);
    }

    int width;
    int mask;
    width = num_thread_x;
    for ( mask = 1 ; mask < width ; mask *=2 ) {
       energy_elec += __shfl_xor(energy_elec, mask, width);
       energy_evdw += __shfl_xor(energy_evdw, mask, width);
       virial(1)   += __shfl_xor(virial(1),   mask, width);
       virial(2)   += __shfl_xor(virial(2),   mask, width);
       virial(3)   += __shfl_xor(virial(3),   mask, width);
    }

    if ( id_thread_x == 0 ) {
       ene_virial_smem(1,id_thread_y+1) = energy_elec;
       ene_virial_smem(2,id_thread_y+1) = energy_evdw;
       ene_virial_smem(3,id_thread_y+1) = virial(1);
       ene_virial_smem(4,id_thread_y+1) = virial(2);
       ene_virial_smem(5,id_thread_y+1) = virial(3);
    }

    __syncthreads();

    if ( id_thread_y == 0 ) {
       energy_elec = ene_virial_smem(1,id_thread_x+1);
       energy_evdw = ene_virial_smem(2,id_thread_x+1);
       virial(1)   = ene_virial_smem(3,id_thread_x+1);
       virial(2)   = ene_virial_smem(4,id_thread_x+1);
       virial(3)   = ene_virial_smem(5,id_thread_x+1);

       width = num_thread_y;
       for ( mask = 1 ; mask < width ; mask *= 2) {
           energy_elec += __shfl_xor (energy_elec, mask, width);
           energy_evdw += __shfl_xor (energy_evdw, mask, width);
           virial(1)   += __shfl_xor (virial(1),   mask, width);
           virial(2)   += __shfl_xor (virial(2),   mask, width);
           virial(3)   += __shfl_xor (virial(3),   mask, width);
       }

       if (id_thread_x == 0 ) {
           ene_virial(1) += energy_elec;
           ene_virial(2) += energy_evdw;
           ene_virial(3) += virial(1);
           ene_virial(4) += virial(2);
           ene_virial(5) += virial(3);
       }
    }
}


/* */

#if defined(_MIXED) || defined(_SINGLE)
#define NUM_CTA__FORCE_INTER_CELL 12
#else
#define NUM_CTA__FORCE_INTER_CELL  9
#endif
/* */
//__launch_bounds__(128,NUM_CTA__FORCE_INTER_CELL)
__global__ void kern_compute_force_nonbond_table_linear_univ__force_inter_cell(
    const REAL  * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:univ_maxcell)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _table_grad,         // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  * _virial_check,        // ( 1:ncel_max, 1:ncel_max )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_ncell_near,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_s,
    int  index_e,
    int  max_iy_natom,
    int  check_virial,
    REAL  density,
    REAL  cutoff2,
    REAL  pairlistdist2,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

#define num_thread_x   32
    // const int  num_thread_x = blockDim.x;
    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];

    /* shared memory */
    REAL  *_warp_smem = & smem[id_thread_y * (max_iy_natom*6)]; // iy_natom*6

#define force_iy_smem(X,Y)  _warp_smem[CALIDX2((X)-1,3, (Y)-1,max_iy_natom)]
#define coord_pbc_smem(X,Y) _warp_smem[CALIDX2((X)-1,3, (Y)-1,max_iy_natom) + (max_iy_natom*3)]

    int  index = index_s + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_e ) return;
    int  univ_ij = univ_ij_sort_list( index );

    int  ix_natom = univ_ix_natom(univ_ij);
    int  iy_natom = univ_iy_natom(univ_ij);
    if ( ix_natom * iy_natom <= 0 ) return;

    const int  i = univ_cell_pairlist1(1,univ_ij);
    const int  j = univ_cell_pairlist1(2,univ_ij);
    int  k;

#define sumval(Z) _sumval[(Z)-1]
    double _sumval[3];
    sumval(1) = 0.0;  // virial(1)
    sumval(2) = 0.0;  // virial(2)
    sumval(3) = 0.0;  // virial(3)
    char  check_virial_ij = 0;
    if ( (check_virial != 0) && (virial_check(j,i) != 0) ) {
	check_virial_ij = 1;
    }

    if (univ_ij > univ_ncell_near) {
	/* far */
        //
        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
    	    int  iiy_e = iiy_s + max_iy_natom - 1;
	    if ( iiy_e > iy_natom ) iiy_e = iy_natom;

	    // initialize force_iy at shared memory
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int  n   = (k % 3) + 1;
	        int  iiy = (k / 3) + iiy_s;
	        force_iy_smem(n, iiy-iiy_s+1) = 0.0;
		int  iy = __ldg(&univ_iy_list(iiy,univ_ij)); /* */
		coord_pbc_smem(n, iiy-iiy_s+1) = __ldg(&coord_pbc(iy, n, j));
	    }

	    for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
	        int  ix = __ldg(&univ_ix_list(iix,univ_ij));

	        REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
	        REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
	        REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
	        REAL  iqtmp   = __ldg(&charge(ix,i));
	        int   iatmcls = __ldg(&atmcls(ix,i));

	        force_local(1) = 0.0;
	        force_local(2) = 0.0;
	        force_local(3) = 0.0;

	        for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
		    int  iy = __ldg(&univ_iy_list(iiy,univ_ij)); /* */

                    REAL  grad_coef = 0.0;
                    REAL  dij1 = 0.0;
                    REAL  dij2 = 0.0;
                    REAL  dij3 = 0.0;
                    REAL  rij2;

                    // dij1 = rtmp1 - __ldg(&coord_pbc(1,iy,j));
                    // dij2 = rtmp2 - __ldg(&coord_pbc(2,iy,j));
                    // dij3 = rtmp3 - __ldg(&coord_pbc(3,iy,j));
                    dij1 = rtmp1 - coord_pbc_smem(1,iiy-iiy_s+1);
                    dij2 = rtmp2 - coord_pbc_smem(2,iiy-iiy_s+1);
                    dij3 = rtmp3 - coord_pbc_smem(3,iiy-iiy_s+1);
                    rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

		    if ( rij2 < cutoff2 ) {

		        rij2 = cutoff2 * density / rij2;

		        REAL  jqtmp   = __ldg(&charge(iy,j));
		        int   jatmcls = __ldg(&atmcls(iy,j));
		        REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
		        REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));

		        int   L  = int(rij2);
		        REAL  R  = rij2 - L;

		        int   L1 = 3*L - 2;
		        REAL  tg0  = __ldg(&table_grad(L1  ));
		        REAL  tg1  = __ldg(&table_grad(L1+1));
		        REAL  tg2  = __ldg(&table_grad(L1+2));
		        REAL  tg3  = __ldg(&table_grad(L1+3));
		        REAL  tg4  = __ldg(&table_grad(L1+4));
		        REAL  tg5  = __ldg(&table_grad(L1+5));
		        REAL  term_lj12 = tg0 + R*(tg3-tg0);
		        REAL  term_lj6  = tg1 + R*(tg4-tg1);
		        REAL  term_elec = tg2 + R*(tg5-tg2);

		        grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;
		    }

		    REAL  work1 = grad_coef*dij1;
		    REAL  work2 = grad_coef*dij2;
		    REAL  work3 = grad_coef*dij3;

		    force_local(1) -= work1;
		    force_local(2) -= work2;
		    force_local(3) -= work3;

		    // update force_iy(:,iiy) at smem
		    WARP_RSUM_345( work1 );
		    WARP_RSUM_345( work2 );
		    WARP_RSUM_345( work3 );
		    if ( id_thread_xx == 0 ) {
		        force_iy_smem(1,iiy-iiy_s+1) += work1;
		        force_iy_smem(2,iiy-iiy_s+1) += work2;
		        force_iy_smem(3,iiy-iiy_s+1) += work3;
		    }
	        }

	        // update virial
                if ( check_virial_ij != 0 ) {
                    sumval(1) += force_local(1);
                    sumval(2) += force_local(2);
                    sumval(3) += force_local(3);
                }

	        // update force(:,:,i)
	        WARP_RSUM_12( force_local(1) );
	        WARP_RSUM_12( force_local(2) );
	        WARP_RSUM_12( force_local(3) );
	        if ( id_thread_xy == 0 ) {
		    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
		    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
		    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
	        }
	    }

	    // __syncthreads();

	    // update force(:,:,j)
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int   n   = (k % 3) + 1;
	        int   iiy = (k / 3) + iiy_s;
	        int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
	        REAL  val = force_iy_smem(n,iiy-iiy_s+1);
	        if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
	    }

        }
    }
    else {
	/* near */
        //
        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
    	    int  iiy_e = iiy_s + max_iy_natom - 1;
	    if ( iiy_e > iy_natom ) iiy_e = iy_natom;

	    // initialize force_iy at shared memory
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int  n   = (k % 3) + 1;
	        int  iiy = (k / 3) + iiy_s;
	        force_iy_smem(n, iiy-iiy_s+1) = 0.0;
		int  iy = __ldg(&univ_iy_list(iiy,univ_ij)); /* */
		coord_pbc_smem(n, iiy-iiy_s+1) = __ldg(&coord_pbc(iy, n, j));
	    }

	    for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
		int  ix = __ldg(&univ_ix_list(iix,univ_ij));

		REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
		REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
		REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
		REAL  iqtmp   = __ldg(&charge(ix,i));
		int   iatmcls = __ldg(&atmcls(ix,i));

		force_local(1) = 0.0;
		force_local(2) = 0.0;
		force_local(3) = 0.0;

		int   idx0 = (ix-1) * univ_natom_max;

		for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
		    int  iy = __ldg(&univ_iy_list(iiy,univ_ij));

		    REAL  grad_coef = 0.0;
		    REAL  dij1 = 0.0;
		    REAL  dij2 = 0.0;
		    REAL  dij3 = 0.0;
		    REAL  rij2;

		    // int idx = iy + (ix-1)*univ_natom_max;
		    int  idx = iy + idx0;
		    if (univ_mask2(idx,univ_ij)) {

			// dij1 = rtmp1 - __ldg(&coord_pbc(1,iy,j));
			// dij2 = rtmp2 - __ldg(&coord_pbc(2,iy,j));
			// dij3 = rtmp3 - __ldg(&coord_pbc(3,iy,j));
			dij1 = rtmp1 - coord_pbc_smem(1,iiy-iiy_s+1);
			dij2 = rtmp2 - coord_pbc_smem(2,iiy-iiy_s+1);
			dij3 = rtmp3 - coord_pbc_smem(3,iiy-iiy_s+1);
			rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

			if ( rij2 < cutoff2 ) {

			    rij2 = cutoff2 * density / rij2;

			    REAL  jqtmp   = __ldg(&charge(iy,j));
			    int   jatmcls = __ldg(&atmcls(iy,j));
			    REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
			    REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));

			    int   L  = int(rij2);
			    REAL  R  = rij2 - L;

			    int   L1 = 3*L - 2;
			    REAL  tg0  = __ldg(&table_grad(L1  ));
			    REAL  tg1  = __ldg(&table_grad(L1+1));
			    REAL  tg2  = __ldg(&table_grad(L1+2));
			    REAL  tg3  = __ldg(&table_grad(L1+3));
			    REAL  tg4  = __ldg(&table_grad(L1+4));
			    REAL  tg5  = __ldg(&table_grad(L1+5));
			    REAL  term_lj12 = tg0 + R*(tg3-tg0);
			    REAL  term_lj6  = tg1 + R*(tg4-tg1);
			    REAL  term_elec = tg2 + R*(tg5-tg2);

			    grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;
			}
		    }
		    REAL  work1 = grad_coef*dij1;
		    REAL  work2 = grad_coef*dij2;
		    REAL  work3 = grad_coef*dij3;

		    force_local(1) -= work1;
		    force_local(2) -= work2;
		    force_local(3) -= work3;

		    // update force_iy(:,iiy) at smem
		    WARP_RSUM_345( work1 );
		    WARP_RSUM_345( work2 );
		    WARP_RSUM_345( work3 );
		    if ( id_thread_xx == 0 ) {
			force_iy_smem(1,iiy-iiy_s+1) += work1;
			force_iy_smem(2,iiy-iiy_s+1) += work2;
			force_iy_smem(3,iiy-iiy_s+1) += work3;
		    }
		}

		// update virial
                if ( check_virial_ij != 0 ) {
                    sumval(1) += force_local(1);
                    sumval(2) += force_local(2);
                    sumval(3) += force_local(3);
                }

		// update force(:,:,i)
		WARP_RSUM_12( force_local(1) );
		WARP_RSUM_12( force_local(2) );
		WARP_RSUM_12( force_local(3) );
		if ( id_thread_xy == 0 ) {
		    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
		    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
		    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
		}
	    }

	    // __syncthreads();

	    // update force(:,:,j)
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int   n   = (k % 3) + 1;
	        int   iiy = (k / 3) + iiy_s;
	        int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
	        REAL  val = force_iy_smem(n,iiy-iiy_s+1);
	        if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
	    }

        }
    }


    // update virial
    if (check_virial != 0) {
        WARP_RSUM_12345( sumval(1) );  // virial(1)
        WARP_RSUM_12345( sumval(2) );  // virial(2)
        WARP_RSUM_12345( sumval(3) );  // virial(3)
        if (id_thread_x < 3) {
            int n = id_thread_x + 1;
            if (n == 1) sumval(n) *= __ldg(&cell_move(n,j,i))*system_x;
            if (n == 2) sumval(n) *= __ldg(&cell_move(n,j,i))*system_y;
            if (n == 3) sumval(n) *= __ldg(&cell_move(n,j,i))*system_z;
            ene_viri_mid(n,index) = sumval(n);
        }
    }
}


__global__ void kern_compute_force_nonbond_table_ljpme_univ__force_inter_cell(
    const REAL  * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:univ_maxcell)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  * _table_grad,         // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  * _virial_check,        // ( 1:ncel_max, 1:ncel_max )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_ncell_near,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_s,
    int  index_e,
    int  max_iy_natom,
    int  check_virial,
    REAL  density,
    REAL  cutoff2,
    REAL  pairlistdist2,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

#define num_thread_x   32
    // const int  num_thread_x = blockDim.x;
    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];

    /* shared memory */
    REAL  *_warp_smem = & smem[id_thread_y * (max_iy_natom*6)]; // iy_natom*6

#define force_iy_smem(X,Y)  _warp_smem[CALIDX2((X)-1,3, (Y)-1,max_iy_natom)]
#define coord_pbc_smem(X,Y) _warp_smem[CALIDX2((X)-1,3, (Y)-1,max_iy_natom) + (max_iy_natom*3)]

    int  index = index_s + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_e ) return;
    int  univ_ij = univ_ij_sort_list( index );

    int  ix_natom = univ_ix_natom(univ_ij);
    int  iy_natom = univ_iy_natom(univ_ij);
    if ( ix_natom * iy_natom <= 0 ) return;

    const int  i = univ_cell_pairlist1(1,univ_ij);
    const int  j = univ_cell_pairlist1(2,univ_ij);
    int  k;

#define sumval(Z) _sumval[(Z)-1]
    double _sumval[3];
    sumval(1) = 0.0;  // virial(1)
    sumval(2) = 0.0;  // virial(2)
    sumval(3) = 0.0;  // virial(3)
    char  check_virial_ij = 0;
    if ( (check_virial != 0) && (virial_check(j,i) != 0) ) {
	check_virial_ij = 1;
    }

    if (univ_ij > univ_ncell_near) {
	/* far */
        //
        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
    	    int  iiy_e = iiy_s + max_iy_natom - 1;
	    if ( iiy_e > iy_natom ) iiy_e = iy_natom;

	    // initialize force_iy at shared memory
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int  n   = (k % 3) + 1;
	        int  iiy = (k / 3) + iiy_s;
	        force_iy_smem(n, iiy-iiy_s+1) = 0.0;
		int  iy = __ldg(&univ_iy_list(iiy,univ_ij)); /* */
		coord_pbc_smem(n, iiy-iiy_s+1) = __ldg(&coord_pbc(iy, n, j));
	    }

	    for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
	        int  ix = __ldg(&univ_ix_list(iix,univ_ij));

	        REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
	        REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
	        REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
	        REAL  iqtmp   = __ldg(&charge(ix,i));
	        int   iatmcls = __ldg(&atmcls(ix,i));
                REAL  lj6_i   = __ldg(&nonb_lj6_factor(iatmcls));

	        force_local(1) = 0.0;
	        force_local(2) = 0.0;
	        force_local(3) = 0.0;

	        for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
		    int  iy = __ldg(&univ_iy_list(iiy,univ_ij)); /* */

                    REAL  grad_coef = 0.0;
                    REAL  dij1 = 0.0;
                    REAL  dij2 = 0.0;
                    REAL  dij3 = 0.0;
                    REAL  rij2;

                    dij1 = rtmp1 - coord_pbc_smem(1,iiy-iiy_s+1);
                    dij2 = rtmp2 - coord_pbc_smem(2,iiy-iiy_s+1);
                    dij3 = rtmp3 - coord_pbc_smem(3,iiy-iiy_s+1);
                    rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

		    if ( rij2 < cutoff2 ) {

                        REAL  inv_r2  = 1.0 / rij2;
                        REAL  inv_r6  = inv_r2 * inv_r2 * inv_r2;
                        REAL  inv_r12 = inv_r6 * inv_r6;
		        rij2 = cutoff2 * density * inv_r2;

		        REAL  jqtmp   = __ldg(&charge(iy,j));
		        int   jatmcls = __ldg(&atmcls(iy,j));
                        REAL  lj6_j   = __ldg(&nonb_lj6_factor(jatmcls));
                        REAL  lj6_ij  = lj6_i * lj6_j;
		        REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
		        REAL  lj6  = __ldg(&nonb_lj6 (jatmcls,iatmcls)) - lj6_ij;

		        int   L  = int(rij2);
		        REAL  R  = rij2 - L;
		        int   L1 = 2*L - 1;

		        REAL  tg0  = __ldg(&table_grad(L1  ));
		        REAL  tg1  = __ldg(&table_grad(L1+1));
		        REAL  tg2  = __ldg(&table_grad(L1+2));
		        REAL  tg3  = __ldg(&table_grad(L1+3));
		        REAL  term_lj6  = tg0 + R*(tg2-tg0);
		        REAL  term_elec = tg1 + R*(tg3-tg1);
                        REAL  term_lj12 = -12.0 * inv_r12 * inv_r2;
                        REAL  term_temp = - 6.0 * inv_r6  * inv_r2;

		        grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij + iqtmp*jqtmp*term_elec;
		    }

		    REAL  work1 = grad_coef*dij1;
		    REAL  work2 = grad_coef*dij2;
		    REAL  work3 = grad_coef*dij3;

		    force_local(1) -= work1;
		    force_local(2) -= work2;
		    force_local(3) -= work3;

		    // update force_iy(:,iiy) at smem
		    WARP_RSUM_345( work1 );
		    WARP_RSUM_345( work2 );
		    WARP_RSUM_345( work3 );
		    if ( id_thread_xx == 0 ) {
		        force_iy_smem(1,iiy-iiy_s+1) += work1;
		        force_iy_smem(2,iiy-iiy_s+1) += work2;
		        force_iy_smem(3,iiy-iiy_s+1) += work3;
		    }
	        }

	        // update virial
                if ( check_virial_ij != 0 ) {
                    sumval(1) += force_local(1);
                    sumval(2) += force_local(2);
                    sumval(3) += force_local(3);
                }

	        // update force(:,:,i)
	        WARP_RSUM_12( force_local(1) );
	        WARP_RSUM_12( force_local(2) );
	        WARP_RSUM_12( force_local(3) );
	        if ( id_thread_xy == 0 ) {
		    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
		    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
		    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
	        }
	    }

	    // __syncthreads();

	    // update force(:,:,j)
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int   n   = (k % 3) + 1;
	        int   iiy = (k / 3) + iiy_s;
	        int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
	        REAL  val = force_iy_smem(n,iiy-iiy_s+1);
	        if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
	    }

        }
    }
    else {
	/* near */
        //
        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
    	    int  iiy_e = iiy_s + max_iy_natom - 1;
	    if ( iiy_e > iy_natom ) iiy_e = iy_natom;

	    // initialize force_iy at shared memory
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int  n   = (k % 3) + 1;
	        int  iiy = (k / 3) + iiy_s;
	        force_iy_smem(n, iiy-iiy_s+1) = 0.0;
		int  iy = __ldg(&univ_iy_list(iiy,univ_ij)); /* */
		coord_pbc_smem(n, iiy-iiy_s+1) = __ldg(&coord_pbc(iy, n, j));
	    }

	    for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
		int  ix = __ldg(&univ_ix_list(iix,univ_ij));

		REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
		REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
		REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
		REAL  iqtmp   = __ldg(&charge(ix,i));
		int   iatmcls = __ldg(&atmcls(ix,i));
                REAL  lj6_i   = __ldg(&nonb_lj6_factor(iatmcls));

		force_local(1) = 0.0;
		force_local(2) = 0.0;
		force_local(3) = 0.0;

		int   idx0 = (ix-1) * univ_natom_max;

		for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
		    int  iy = __ldg(&univ_iy_list(iiy,univ_ij));

		    REAL  grad_coef = 0.0;
		    REAL  dij1 = 0.0;
		    REAL  dij2 = 0.0;
		    REAL  dij3 = 0.0;
		    REAL  rij2;

		    // int idx = iy + (ix-1)*univ_natom_max;
		    int  idx = iy + idx0;
		    if (univ_mask2(idx,univ_ij)) {

			dij1 = rtmp1 - coord_pbc_smem(1,iiy-iiy_s+1);
			dij2 = rtmp2 - coord_pbc_smem(2,iiy-iiy_s+1);
			dij3 = rtmp3 - coord_pbc_smem(3,iiy-iiy_s+1);
			rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

			if ( rij2 < cutoff2 ) {

                            REAL  inv_r2  = 1.0 / rij2;
                            REAL  inv_r6  = inv_r2 * inv_r2 * inv_r2;
                            REAL  inv_r12 = inv_r6 * inv_r6;
			    rij2 = cutoff2 * density * inv_r2;

			    REAL  jqtmp   = __ldg(&charge(iy,j));
			    int   jatmcls = __ldg(&atmcls(iy,j));
                            REAL  lj6_j   = __ldg(&nonb_lj6_factor(jatmcls));
                            REAL  lj6_ij  = lj6_i * lj6_j;
			    REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
			    REAL  lj6  = __ldg(&nonb_lj6 (jatmcls,iatmcls)) - lj6_ij;

			    int   L  = int(rij2);
			    REAL  R  = rij2 - L;

			    int   L1 = 2*L - 1;
			    REAL  tg0  = __ldg(&table_grad(L1  ));
			    REAL  tg1  = __ldg(&table_grad(L1+1));
			    REAL  tg2  = __ldg(&table_grad(L1+2));
			    REAL  tg3  = __ldg(&table_grad(L1+3));
			    REAL  term_lj6  = tg0 + R*(tg2-tg0);
			    REAL  term_elec = tg1 + R*(tg3-tg1);
			    REAL  term_lj12 = -12.0 * inv_r12 * inv_r2;
                            REAL  term_temp = - 6.0 * inv_r6  * inv_r2;

			    grad_coef = term_lj12*lj12 - term_temp*lj6 - term_lj6*lj6_ij + iqtmp*jqtmp*term_elec;
			}
		    }
		    REAL  work1 = grad_coef*dij1;
		    REAL  work2 = grad_coef*dij2;
		    REAL  work3 = grad_coef*dij3;

		    force_local(1) -= work1;
		    force_local(2) -= work2;
		    force_local(3) -= work3;

		    // update force_iy(:,iiy) at smem
		    WARP_RSUM_345( work1 );
		    WARP_RSUM_345( work2 );
		    WARP_RSUM_345( work3 );
		    if ( id_thread_xx == 0 ) {
			force_iy_smem(1,iiy-iiy_s+1) += work1;
			force_iy_smem(2,iiy-iiy_s+1) += work2;
			force_iy_smem(3,iiy-iiy_s+1) += work3;
		    }
		}

		// update virial
                if ( check_virial_ij != 0 ) {
                    sumval(1) += force_local(1);
                    sumval(2) += force_local(2);
                    sumval(3) += force_local(3);
                }

		// update force(:,:,i)
		WARP_RSUM_12( force_local(1) );
		WARP_RSUM_12( force_local(2) );
		WARP_RSUM_12( force_local(3) );
		if ( id_thread_xy == 0 ) {
		    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
		    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
		    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
		}
	    }

	    // __syncthreads();

	    // update force(:,:,j)
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int   n   = (k % 3) + 1;
	        int   iiy = (k / 3) + iiy_s;
	        int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
	        REAL  val = force_iy_smem(n,iiy-iiy_s+1);
	        if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
	    }

        }
    }


    // update virial
    if (check_virial != 0) {
        WARP_RSUM_12345( sumval(1) );  // virial(1)
        WARP_RSUM_12345( sumval(2) );  // virial(2)
        WARP_RSUM_12345( sumval(3) );  // virial(3)
        if (id_thread_x < 3) {
            int n = id_thread_x + 1;
            if (n == 1) sumval(n) *= __ldg(&cell_move(n,j,i))*system_x;
            if (n == 2) sumval(n) *= __ldg(&cell_move(n,j,i))*system_y;
            if (n == 3) sumval(n) *= __ldg(&cell_move(n,j,i))*system_z;
            ene_viri_mid(n,index) = sumval(n);
        }
    }
}


#if defined(_MIXED) || defined(_SINGLE)
#define NUM_CTA__FORCE_INTER_CELL 12
#else
#define NUM_CTA__FORCE_INTER_CELL  9
#endif
/* */
__launch_bounds__(128,NUM_CTA__FORCE_INTER_CELL)
__global__ void kern_compute_force_nonbond_notable_univ__force_inter_cell(
    const REAL  * _coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        * _force,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    double      * _ene_virial,          // ( 5 )
    double      * _ene_viri_mid,        // ( 1:5, 1:univ_maxcell)
    const char  * _cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  * _charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   * _natom,               // ( 1:ncel_max )
    const REAL  * _nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  * _table_grad,         // ( 1:6*cutoff_int )
    const int   * _univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  * _univ_mask2,
    const int   * _univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar * _univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   * _univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  * _virial_check,        // ( 1:ncel_max, 1:ncel_max )
    int  MaxAtom,
    int  MaxAtomCls,
    int  num_atom_cls,
    int  ncel_local,
    int  ncel_bound,
    int  ncel_max,
    int  cutoff_int,
    int  univ_maxcell,
    int  univ_maxcell1,
    int  univ_ncell_near,
    int  univ_mask2_size,
    int  univ_natom_max,
    int  index_s,
    int  index_e,
    int  max_iy_natom,
    int  check_virial,
    REAL  density,
    REAL  cutoff2,
    REAL  pairlistdist2,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

#define num_thread_x   32
    // const int  num_thread_x = blockDim.x;
    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;
    REAL  _force_local[3];

    /* shared memory */
    REAL  *_warp_smem = & smem[id_thread_y * (max_iy_natom*6)]; // iy_natom*6

#define force_iy_smem(X,Y)  _warp_smem[CALIDX2((X)-1,3, (Y)-1,max_iy_natom)]
#define coord_pbc_smem(X,Y) _warp_smem[CALIDX2((X)-1,3, (Y)-1,max_iy_natom) + (max_iy_natom*3)]

    int  index = index_s + id_thread_y + (num_thread_y * blockIdx.x);
    if ( index > index_e ) return;
    int  univ_ij = univ_ij_sort_list( index );

    int  ix_natom = univ_ix_natom(univ_ij);
    int  iy_natom = univ_iy_natom(univ_ij);
    if ( ix_natom * iy_natom <= 0 ) return;

    const int  i = univ_cell_pairlist1(1,univ_ij);
    const int  j = univ_cell_pairlist1(2,univ_ij);
    int  k;

#define sumval(Z) _sumval[(Z)-1]
    double _sumval[3];
    sumval(1) = 0.0;  // virial(1)
    sumval(2) = 0.0;  // virial(2)
    sumval(3) = 0.0;  // virial(3)
    char  check_virial_ij = 0;
    if ( (check_virial != 0) && (virial_check(j,i) != 0) ) {
	check_virial_ij = 1;
    }

    if (univ_ij > univ_ncell_near) {
	/* far */
        //
        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
    	    int  iiy_e = iiy_s + max_iy_natom - 1;
	    if ( iiy_e > iy_natom ) iiy_e = iy_natom;

	    // initialize force_iy at shared memory
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int  n   = (k % 3) + 1;
	        int  iiy = (k / 3) + iiy_s;
	        force_iy_smem(n, iiy-iiy_s+1) = 0.0;
		int  iy = __ldg(&univ_iy_list(iiy,univ_ij)); /* */
		coord_pbc_smem(n, iiy-iiy_s+1) = __ldg(&coord_pbc(iy, n, j));
	    }

	    for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
	        int  ix = __ldg(&univ_ix_list(iix,univ_ij));

	        REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
	        REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
	        REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
	        REAL  iqtmp   = __ldg(&charge(ix,i));
	        int   iatmcls = __ldg(&atmcls(ix,i));

	        force_local(1) = 0.0;
	        force_local(2) = 0.0;
	        force_local(3) = 0.0;

	        for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
		    int  iy = __ldg(&univ_iy_list(iiy,univ_ij)); /* */

                    REAL  grad_coef = 0.0;
                    REAL  dij1 = 0.0;
                    REAL  dij2 = 0.0;
                    REAL  dij3 = 0.0;
                    REAL  rij2;

                    // dij1 = rtmp1 - __ldg(&coord_pbc(1,iy,j));
                    // dij2 = rtmp2 - __ldg(&coord_pbc(2,iy,j));
                    // dij3 = rtmp3 - __ldg(&coord_pbc(3,iy,j));
                    dij1 = rtmp1 - coord_pbc_smem(1,iiy-iiy_s+1);
                    dij2 = rtmp2 - coord_pbc_smem(2,iiy-iiy_s+1);
                    dij3 = rtmp3 - coord_pbc_smem(3,iiy-iiy_s+1);
                    rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

		    if ( rij2 < cutoff2 ) {

                        REAL rij2_inv = 1.0 / rij2;
		        rij2 = cutoff2 * density * rij2_inv;

		        REAL  jqtmp   = __ldg(&charge(iy,j));
		        int   jatmcls = __ldg(&atmcls(iy,j));
		        REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
		        REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));

		        int   L  = int(rij2);
		        REAL  R  = rij2 - L;

		        REAL  tg0  = __ldg(&table_grad(L  ));
		        REAL  tg1  = __ldg(&table_grad(L+1));
		        REAL  term_elec = tg0 + R*(tg1-tg0);
                        REAL  term_lj6  = rij2_inv * rij2_inv * rij2_inv;
                        REAL  term_lj12 = term_lj6 * term_lj6;
                        term_lj12 = -12.0 * term_lj12 * rij2_inv;
                        term_lj6  = - 6.0 * term_lj6  * rij2_inv;

		        grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;
		    }

		    REAL  work1 = grad_coef*dij1;
		    REAL  work2 = grad_coef*dij2;
		    REAL  work3 = grad_coef*dij3;

		    force_local(1) -= work1;
		    force_local(2) -= work2;
		    force_local(3) -= work3;

		    // update force_iy(:,iiy) at smem
		    WARP_RSUM_345( work1 );
		    WARP_RSUM_345( work2 );
		    WARP_RSUM_345( work3 );
		    if ( id_thread_xx == 0 ) {
		        force_iy_smem(1,iiy-iiy_s+1) += work1;
		        force_iy_smem(2,iiy-iiy_s+1) += work2;
		        force_iy_smem(3,iiy-iiy_s+1) += work3;
		    }
	        }

	        // update virial
                if ( check_virial_ij != 0 ) {
                    sumval(1) += force_local(1);
                    sumval(2) += force_local(2);
                    sumval(3) += force_local(3);
                }

	        // update force(:,:,i)
	        WARP_RSUM_12( force_local(1) );
	        WARP_RSUM_12( force_local(2) );
	        WARP_RSUM_12( force_local(3) );
	        if ( id_thread_xy == 0 ) {
		    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
		    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
		    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
	        }
	    }

	    // __syncthreads();

	    // update force(:,:,j)
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int   n   = (k % 3) + 1;
	        int   iiy = (k / 3) + iiy_s;
	        int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
	        REAL  val = force_iy_smem(n,iiy-iiy_s+1);
	        if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
	    }

        }
    }
    else {
	/* near */
        //
        for ( int iiy_s = 1; iiy_s <= iy_natom; iiy_s += max_iy_natom ) { // blocking
    	    int  iiy_e = iiy_s + max_iy_natom - 1;
	    if ( iiy_e > iy_natom ) iiy_e = iy_natom;

	    // initialize force_iy at shared memory
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int  n   = (k % 3) + 1;
	        int  iiy = (k / 3) + iiy_s;
	        force_iy_smem(n, iiy-iiy_s+1) = 0.0;
		int  iy = __ldg(&univ_iy_list(iiy,univ_ij)); /* */
		coord_pbc_smem(n, iiy-iiy_s+1) = __ldg(&coord_pbc(iy, n, j));
	    }

	    for ( int iix = id_thread_xx + 1; iix <= ix_natom; iix += num_thread_xx ) {
		int  ix = __ldg(&univ_ix_list(iix,univ_ij));

		REAL  rtmp1   = __ldg(&coord_pbc(ix,1,i)) + __ldg(&cell_move(1,j,i))*system_x;
		REAL  rtmp2   = __ldg(&coord_pbc(ix,2,i)) + __ldg(&cell_move(2,j,i))*system_y;
		REAL  rtmp3   = __ldg(&coord_pbc(ix,3,i)) + __ldg(&cell_move(3,j,i))*system_z;
		REAL  iqtmp   = __ldg(&charge(ix,i));
		int   iatmcls = __ldg(&atmcls(ix,i));

		force_local(1) = 0.0;
		force_local(2) = 0.0;
		force_local(3) = 0.0;

		int   idx0 = (ix-1) * univ_natom_max;

		for ( int iiy = id_thread_xy + iiy_s; iiy <= iiy_e; iiy += num_thread_xy ) {
		    int  iy = __ldg(&univ_iy_list(iiy,univ_ij));

		    REAL  grad_coef = 0.0;
		    REAL  dij1 = 0.0;
		    REAL  dij2 = 0.0;
		    REAL  dij3 = 0.0;
		    REAL  rij2;

		    // int idx = iy + (ix-1)*univ_natom_max;
		    int  idx = iy + idx0;
		    if (univ_mask2(idx,univ_ij)) {

			// dij1 = rtmp1 - __ldg(&coord_pbc(1,iy,j));
			// dij2 = rtmp2 - __ldg(&coord_pbc(2,iy,j));
			// dij3 = rtmp3 - __ldg(&coord_pbc(3,iy,j));
			dij1 = rtmp1 - coord_pbc_smem(1,iiy-iiy_s+1);
			dij2 = rtmp2 - coord_pbc_smem(2,iiy-iiy_s+1);
			dij3 = rtmp3 - coord_pbc_smem(3,iiy-iiy_s+1);
			rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;

			if ( rij2 < cutoff2 ) {

                            REAL rij2_inv = 1.0 / rij2;
                            rij2 = cutoff2 * density * rij2_inv;
  
                            REAL  jqtmp   = __ldg(&charge(iy,j));
                            int   jatmcls = __ldg(&atmcls(iy,j));
                            REAL  lj12 = __ldg(&nonb_lj12(jatmcls,iatmcls));
                            REAL  lj6  = __ldg(&nonb_lj6( jatmcls,iatmcls));
  
                            int   L  = int(rij2);
                            REAL  R  = rij2 - L;
  
                            REAL  tg0  = __ldg(&table_grad(L  ));
                            REAL  tg1  = __ldg(&table_grad(L+1));
                            REAL  term_elec = tg0 + R*(tg1-tg0);
                            REAL  term_lj6  = rij2_inv * rij2_inv * rij2_inv;
                            REAL  term_lj12 = term_lj6 * term_lj6;
                            term_lj12 = -12.0 * term_lj12 * rij2_inv;
                            term_lj6  = - 6.0 * term_lj6  * rij2_inv;

			    grad_coef = term_lj12*lj12 - term_lj6*lj6 + iqtmp*jqtmp*term_elec;
			}
		    }
		    REAL  work1 = grad_coef*dij1;
		    REAL  work2 = grad_coef*dij2;
		    REAL  work3 = grad_coef*dij3;

		    force_local(1) -= work1;
		    force_local(2) -= work2;
		    force_local(3) -= work3;

		    // update force_iy(:,iiy) at smem
		    WARP_RSUM_345( work1 );
		    WARP_RSUM_345( work2 );
		    WARP_RSUM_345( work3 );
		    if ( id_thread_xx == 0 ) {
			force_iy_smem(1,iiy-iiy_s+1) += work1;
			force_iy_smem(2,iiy-iiy_s+1) += work2;
			force_iy_smem(3,iiy-iiy_s+1) += work3;
		    }
		}

		// update virial
                if ( check_virial_ij != 0 ) {
                    sumval(1) += force_local(1);
                    sumval(2) += force_local(2);
                    sumval(3) += force_local(3);
                }

		// update force(:,:,i)
		WARP_RSUM_12( force_local(1) );
		WARP_RSUM_12( force_local(2) );
		WARP_RSUM_12( force_local(3) );
		if ( id_thread_xy == 0 ) {
		    if ( force_local(1) != 0.0 ) atomicAdd( &(force(ix,1,i)), force_local(1) );
		    if ( force_local(2) != 0.0 ) atomicAdd( &(force(ix,2,i)), force_local(2) );
		    if ( force_local(3) != 0.0 ) atomicAdd( &(force(ix,3,i)), force_local(3) );
		}
	    }

	    // __syncthreads();

	    // update force(:,:,j)
	    for ( k = id_thread_x; k < 3 * (iiy_e - iiy_s + 1); k += num_thread_x ) {
	        int   n   = (k % 3) + 1;
	        int   iiy = (k / 3) + iiy_s;
	        int   iy  = __ldg(&univ_iy_list(iiy,univ_ij));
	        REAL  val = force_iy_smem(n,iiy-iiy_s+1);
	        if ( val != 0.0 ) atomicAdd( &(force(iy,n,j)), val );
	    }

        }
    }

    // update virial
    if (check_virial != 0) {
        WARP_RSUM_12345( sumval(1) );  // virial(1)
        WARP_RSUM_12345( sumval(2) );  // virial(2)
        WARP_RSUM_12345( sumval(3) );  // virial(3)
        if (id_thread_x < 3) {
            int n = id_thread_x + 1;
            if (n == 1) sumval(n) *= __ldg(&cell_move(n,j,i))*system_x;
            if (n == 2) sumval(n) *= __ldg(&cell_move(n,j,i))*system_y;
            if (n == 3) sumval(n) *= __ldg(&cell_move(n,j,i))*system_z;
            ene_viri_mid(n,index) = sumval(n);
        }
    }
}

__global__ void kern_compute_force_nonbond_table_linear_univ_sum(
    double       *_ene_virial,
    const double *_ene_viri_mid,
    int          ncel_local,
    int          ncel_max,
    int          univ_maxcell,
    int          univ_gpu_start,
    int          univ_ncell_nonzero
    )
{
#undef  num_thread_xx
#undef  num_thread_xy
#undef  id_thread_xx
#undef  id_thread_xy
#define num_thread_xx   8
#define num_thread_xy   4
#define id_thread_xx  (id_thread_x / num_thread_xy)
#define id_thread_xy  (id_thread_x & (num_thread_xy-1))

#define num_thread_x   32
    // const int  num_thread_x = blockDim.x;
    const int  num_thread_y = blockDim.y;
    const int  id_thread_x = threadIdx.x;
    const int  id_thread_y = threadIdx.y;

    const int num_thread   = ( blockDim.x * blockDim.y );
    const int id_thread    = ( threadIdx.x + blockDim.x * threadIdx.y );

    __shared__ double _virial_smem[3*32];
#define virial_smem(Y,Z)  _virial_smem[CALIDX2((Y)-1,3, (Z)-1,32)]

    double _virial[3];
    virial(1) = 0.0;
    virial(2) = 0.0;
    virial(3) = 0.0;

    int ij;
    for ( ij = id_thread+univ_gpu_start+1 ; ij <= univ_ncell_nonzero ; ij += num_thread ) {
       virial(1)   += ene_viri_mid(1,ij);
       virial(2)   += ene_viri_mid(2,ij);
       virial(3)   += ene_viri_mid(3,ij);
    }

    int width;
    int mask;
    width = num_thread_x;
    for ( mask = 1 ; mask < width ; mask *=2 ) {
       virial(1)   += __shfl_xor(virial(1),   mask, width);
       virial(2)   += __shfl_xor(virial(2),   mask, width);
       virial(3)   += __shfl_xor(virial(3),   mask, width);
    }

    if ( id_thread_x == 0 ) {
       virial_smem(1,id_thread_y+1) = virial(1);
       virial_smem(2,id_thread_y+1) = virial(2);
       virial_smem(3,id_thread_y+1) = virial(3);
    }

    __syncthreads();

    if ( id_thread_y == 0 ) {
       virial(1)   = virial_smem(1,id_thread_x+1);
       virial(2)   = virial_smem(2,id_thread_x+1);
       virial(3)   = virial_smem(3,id_thread_x+1);

       width = num_thread_y;
       for ( mask = 1 ; mask < width ; mask *= 2) {
           virial(1)   += __shfl_xor (virial(1),   mask, width);
           virial(2)   += __shfl_xor (virial(2),   mask, width);
           virial(3)   += __shfl_xor (virial(3),   mask, width);
       }

       if (id_thread_x == 0 ) {
           ene_virial(1) += virial(1);
           ene_virial(2) += virial(2);
           ene_virial(3) += virial(3);
       }
    }
}

/*
 *
 */
__device__ int ceil_pow2( int num )
{
    int ans = num;
    ans -= 1;
    ans |= (ans >>  1);
    ans |= (ans >>  2);
    ans |= (ans >>  4);
    ans |= (ans >>  8);
    ans |= (ans >> 16);
    ans += 1;
    return ans;
}

/*
 *
 */
#if defined(_MIXED) || defined(_SINGLE)
#define NUM_CTA__BUILD_PAIRLIST 16
#else
#define NUM_CTA__BUILD_PAIRLIST 12
#endif
/* */
template <int NUM_WARP, int MAX_ATOM>
//__launch_bounds__(128,NUM_CTA__BUILD_PAIRLIST)
__global__ void kern_build_pairlist(
    const REAL  *_coord_pbc,            // ( 1:MaxAtom, 1:3, 1:ncel_max ), AKA trans2
    const char  *_cell_move,            // ( 1:3, 1:ncel_max, 1:ncel_max )
    const int   *_natom,                // ( 1:ncel_max )
    const int   *_univ_cell_pairlist1,  // ( 1:2, 1:univ_maxcell )
    uchar       *_univ_ix_list,         // ( 1:MaxAtom, 1:univ_maxcell1? )
    uchar       *_univ_iy_list,         // ( 1:MaxAtom, 1:univ_maxcell1? )
    int         *_univ_ix_natom,        // ( 1:univ_maxcell1? )
    int         *_univ_iy_natom,        // ( 1:univ_maxcell1? )
    int         MaxAtom,
    int         ncel_local,
    int         ncel_bound,
    int         ncel_max,
    int         univ_maxcell,
    int         univ_maxcell1,
    REAL        pairdist2,
    REAL        cutoffdist2,
    REAL        system_x,
    REAL        system_y,
    REAL        system_z
    )
{
    const int  w_id = threadIdx.y;
    const int  univ_ij = 1 + w_id + (blockIdx.x * NUM_WARP);
    if ( univ_ij > univ_maxcell ) return;

    const int  i = univ_cell_pairlist1(1,univ_ij);
    const int  j = univ_cell_pairlist1(2,univ_ij);

    const int  t_id = threadIdx.x;
    const int  t_num = blockDim.x;

    const int  t_num_x = 8; /* do not change */
    const int  t_num_y = 4; /* do not change */
    const int  t_id_x = t_id / t_num_y;
    const int  t_id_y = t_id % t_num_y;

    __shared__ ushort  smem_ix_list[NUM_WARP][MAX_ATOM]; // hi  8-bit: counter
    __shared__ ushort  smem_iy_list[NUM_WARP][MAX_ATOM]; // low 8-bit: atom id
    __shared__ int  smem_max[NUM_WARP];

    int  ix_natom_p2 = ceil_pow2(natom(i));
    int  iy_natom_p2 = ceil_pow2(natom(j));

    for ( int ix = 1 + t_id; ix <= ix_natom_p2 ; ix += t_num ) {
	if (ix <= natom(i)) smem_ix_list[w_id][ix-1] = ix;
	else  	            smem_ix_list[w_id][ix-1] = 0;
    }
    for ( int iy = 1 + t_id; iy <= iy_natom_p2 ; iy += t_num ) {
	if (iy <= natom(j)) smem_iy_list[w_id][iy-1] = iy;
	else 	            smem_iy_list[w_id][iy-1] = 0;
    }

    for ( int ix = 1 + t_id_x; ix <= natom(i); ix += t_num_x ) {
	REAL rtmp1 = __ldg(&cell_move(1,j,i))*system_x + __ldg(&coord_pbc(ix,1,i));
	REAL rtmp2 = __ldg(&cell_move(2,j,i))*system_y + __ldg(&coord_pbc(ix,2,i));
	REAL rtmp3 = __ldg(&cell_move(3,j,i))*system_z + __ldg(&coord_pbc(ix,3,i));
	int  val_x = 0;
	for ( int iy = 1 + t_id_y; iy <= natom(j); iy += t_num_y ) {
	    REAL dij1 = rtmp1 - __ldg(&coord_pbc(iy,1,j));
	    REAL dij2 = rtmp2 - __ldg(&coord_pbc(iy,2,j));
	    REAL dij3 = rtmp3 - __ldg(&coord_pbc(iy,3,j));
	    REAL rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3;
	    int  val_y = 0;
	    if (rij2 < pairdist2) {
		val_x += (1 << 8);
		val_y += (1 << 8);
	    }
#if 0
	    if (rij2 < cutoffdist2) {
		val_x += (1 << 8);
		val_y += (1 << 8);
	    }
#endif
	    WARP_RSUM_345( val_y );
	    if ( t_id_x == 0 ) {
		smem_iy_list[w_id][iy-1] += val_y;
	    }
	}
	WARP_RSUM_12( val_x );
	if ( t_id_y == 0 ) {
	    smem_ix_list[w_id][ix-1] += val_x;
	}
    }

    /* sort list x (bitonic sort) */
    for ( int phase = 1; phase < ix_natom_p2; phase *= 2 ) {
	for ( int step = phase; step >= 1; step /= 2 ) {
	    int mask = step - 1;
	    int ofst = step;
	    if ( phase == step ) {
		ofst = step * 2 - 1;
	    }
	    for ( int i0 = t_id; i0 < ix_natom_p2 / 2; i0 += t_num ) {
		int ii = (i0 & mask) + ((i0 & ~mask) << 1);
		int jj = ii ^ ofst;
		ushort val_ii = smem_ix_list[w_id][ii];
		ushort val_jj = smem_ix_list[w_id][jj];
		if ( val_ii < val_jj ) {
		    smem_ix_list[w_id][ii] = val_jj;
		    smem_ix_list[w_id][jj] = val_ii;
		}
	    }
	}
    }
    int  ix_natom = 0;
    for ( int ix = 1 + t_id; ix <= natom(i); ix += t_num ) {
	if ( smem_ix_list[w_id][ix-1] >= (1<<8) ) {
	    univ_ix_list(ix, univ_ij) = smem_ix_list[w_id][ix-1] & 0xff;
	    ix_natom = ix;
	}
    }

    smem_max[w_id] = 0;
    if ( ix_natom > 0 ) {
	atomicMax( & smem_max[w_id], ix_natom );
    }
    if ( t_id == 0 ) {
	univ_ix_natom(univ_ij) = smem_max[w_id];
    }

    /* sort list y (bitonic sort) */
    for ( int phase = 1; phase < iy_natom_p2; phase *= 2 ) {
	for ( int step = phase; step > 0; step /= 2 ) {
	    int mask = step - 1;
	    int ofst = step;
	    if ( phase == step ) {
		ofst = step * 2 - 1;
	    }
	    for ( int i0 = t_id; i0 < iy_natom_p2 / 2; i0 += t_num ) {
		int ii = (i0 & mask) + ((i0 & ~mask) << 1);
		int jj = ii ^ ofst;
		ushort val_ii = smem_iy_list[w_id][ii];
		ushort val_jj = smem_iy_list[w_id][jj];
		if ( val_ii < val_jj ) {
		    smem_iy_list[w_id][ii] = val_jj;
		    smem_iy_list[w_id][jj] = val_ii;
		}
	    }
	}
    }
    int  iy_natom = 0;
    for ( int iy = 1 + t_id; iy <= natom(j); iy += t_num ) {
	if ( smem_iy_list[w_id][iy-1] >= (1<<8) ) {
	    univ_iy_list(iy, univ_ij) = smem_iy_list[w_id][iy-1] & 0xff;
	    iy_natom = iy;
	}
    }

    smem_max[w_id] = 0;
    if ( iy_natom > 0 ) {
	atomicMax( & smem_max[w_id], iy_natom );
    }
    if ( t_id == 0 ) {
	univ_iy_natom(univ_ij) = smem_max[w_id];
    }
}

/*
 *
 */
double cur_time()
{
    struct timeval  tv;
    gettimeofday( &tv, NULL );
    return( (double)tv.tv_sec + (double)tv.tv_usec * 1e-6 );
}

/*
 * JJ : energy calculation with linear lookup table
 */
void gpu_launch_compute_energy_nonbond_table_linear_univ(
    REAL    *_coord_pbc,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL    *_force,                   // ( 1:MaxAtom, 1:3, 1:ncell_all, 1:nthread )
    double  *_ene_virial,
    const REALI *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_ene,           // ( 1:6*cutoff_int )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:ncel_max, 1:ncel_max)
    int   MaxAtom,
    int   MaxAtomCls,
    int   num_atom_cls,
    int   ncel_local,
    int   ncel_bound,
    int   ncel_max,
    int   cutoff_int,
    int   univ_maxcell,
    int   univ_maxcell1,
    int   univ_ncell_nonzero,
    int   univ_ncell_near,
    int   univ_update,
    int   univ_mask2_size,
    int   univ_natom_max,
    int   check_virial,
    REAL  density,
    REAL  cutoff2,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{

#ifdef DEBUG
    printf( "[%s,%s,%d]\n", __FILE__, __func__, __LINE__ );
    printf( "MaxAtom = %d\n", MaxAtom );
    printf( "MaxAtomCls = %d\n", MaxAtomCls );
    printf( "num_atom_cls = %d\n", num_atom_cls );
    printf( "ncel_local = %d\n", ncel_local );
    printf( "ncel_bound = %d\n", ncel_bound );
    printf( "ncel_max   = %d\n", ncel_max );
    printf( "cutoff_int = %d\n", cutoff_int );
    printf( "univ_maxcell = %d\n", univ_maxcell );
    printf( "univ_maxcell1 = %d\n", univ_maxcell1 );
    printf( "univ_ncell_nonzero = %d\n", univ_ncell_nonzero );
    printf( "univ_update = %d\n", univ_update );
    printf( "check_virial = %d\n", check_virial );
    printf( "cutoff2 = %lf\n", cutoff2 );
#endif

    static int first = 1;

    // memory allocation
    if ( first ) {
        gpu_init_buffer(
            _coord_pbc,
            _force,
            _ene_virial,
            _coord,
            _trans1,
            _cell_move,
            _charge,
            _atmcls,
            _natom,
            _nonb_lj12,
            _nonb_lj6,
            _nonb_lj6_factor,
            _table_ene,
            _table_grad,
            _univ_cell_pairlist1,
            _univ_mask2,
            _univ_ix_natom,
            _univ_ix_list,
            _univ_iy_natom,
            _univ_iy_list,
            _virial_check,
            MaxAtom,
            MaxAtomCls,
            num_atom_cls,
            ncel_local,
            ncel_bound,
            ncel_max,
            cutoff_int,
            univ_maxcell,
            univ_maxcell1,
            univ_ncell_near,
            univ_mask2_size
            );
    }

    // send data to GPU
    gpu_memcpy_h2d_energy(
        _coord_pbc,
        _force,
        _ene_virial,
        _coord,
        _trans1,
        _cell_move,
        _charge,
        _atmcls,
        _natom,
        _nonb_lj12,
        _nonb_lj6,
        _nonb_lj6_factor,
        _table_ene,
        _table_grad,
        _univ_cell_pairlist1,
        _univ_mask2,
        _univ_ix_natom,
        _univ_ix_list,
        _univ_iy_natom,
        _univ_iy_list,
        _univ_ij_sort_list,
        _virial_check,
        MaxAtom,
        univ_maxcell,
        univ_ncell_near,
        univ_mask2_size,
        univ_update,
        check_virial,
        first
        );

    first = 0;

    // return; // debug

    /* */
    dim3  def_dim3(1,1,1);
    dim3  num_block;
    dim3  num_thread;
    int  index_start, index_end;

    /* */
    num_block = def_dim3;
    num_block.x = ncel_local + ncel_bound;
    num_thread = def_dim3;
    num_thread.x = 256;

    /* coordinates in periodic boundary condition is written */
    kern_compute_force_nonbond_table_linear_univ__coord_pbc<<< num_block, num_thread, 0, stream[0] >>>(
        dev_coord_pbc,
        dev_coord,
        dev_trans1,
        dev_natom,
        MaxAtom,
        ncel_local,
        ncel_bound,
        ncel_max );

    /* energy/force calculation within a cell */
    index_start = 1;
    index_end  = ncel_local;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_end - index_start + 1), num_thread.y);

    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[1], 0 ) );
    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[2], 0 ) ); /* univ_mask2 */

    // CUDA_CALL( cudaStreamSynchronize( stream[1] ) );

    /* energy calculation */
    kern_compute_energy_nonbond_table_linear_univ__energyforce_intra_cell<<< num_block, num_thread, 0, stream[0] >>>(
        dev_coord_pbc,
        dev_force,
        dev_ene_virial,
        dev_ene_viri_mid,
        dev_cell_move,
        dev_charge,
        dev_atmcls,
        dev_natom,
        dev_nonb_lj12,
        dev_nonb_lj6,
        dev_table_ene,
        dev_table_grad,
        dev_univ_cell_pairlist1,
        dev_univ_mask2,
        dev_univ_ix_natom,
        dev_univ_ix_list,
        dev_univ_iy_natom,
        dev_univ_iy_list,
        dev_univ_ij_sort_list,
        MaxAtom,
        MaxAtomCls,
        num_atom_cls,
        ncel_local,
        ncel_bound,
        ncel_max,
        cutoff_int,
        univ_maxcell,
        univ_maxcell1,
        univ_mask2_size,
        univ_natom_max,
        index_start, index_end,
        density,
        cutoff2 );

    /* energy/force calculation between different near cells */
    index_start = ncel_local + 1;
    index_end = univ_ncell_nonzero;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;  /* do not change */
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_end - index_start + 1), num_thread.y);

    int max_iy_natom = (SMEM_SIZE/sizeof(REAL)/3) / (num_thread.y*NUM_CTA__ENERGYFORCE_INTER_CELL);
    max_iy_natom -= max_iy_natom % 4;

    size_t  smem_size = num_thread.y * max_iy_natom * 3 * sizeof(REAL);

    kern_compute_energy_nonbond_table_linear_univ__energyforce_inter_cell<<< num_block, num_thread, smem_size, stream[0] >>>(
	dev_coord_pbc,
	dev_force,
	dev_ene_virial,
	dev_ene_viri_mid,
	dev_cell_move,
	dev_charge,
	dev_atmcls,
	dev_natom,
	dev_nonb_lj12,
	dev_nonb_lj6,
	dev_table_ene,
	dev_table_grad,
	dev_univ_cell_pairlist1,
        dev_univ_mask2,
	dev_univ_ix_natom,
	dev_univ_ix_list,
	dev_univ_iy_natom,
	dev_univ_iy_list,
	dev_univ_ij_sort_list,
	dev_virial_check,
	MaxAtom,
	MaxAtomCls,
	num_atom_cls,
	ncel_local,
	ncel_bound,
	ncel_max,
	cutoff_int,
	univ_maxcell,
	univ_maxcell1,
        univ_ncell_near,
        univ_mask2_size,
        univ_natom_max,
	index_start, index_end, max_iy_natom,
	density,
	cutoff2,
        system_x, system_y, system_z );
    /* */
    num_block = def_dim3;
    num_thread = def_dim3;
    num_thread.x = 32;
    num_thread.y = 32;
    kern_compute_energy_nonbond_table_linear_univ_energy_sum<<< num_block, num_thread, 0, stream[0] >>>(
	dev_ene_virial,
	dev_ene_viri_mid,
	ncel_local,
	ncel_max,
	univ_maxcell,
	univ_ncell_nonzero);

}


void gpu_launch_compute_energy_nonbond_table_ljpme_univ(
    REAL    *_coord_pbc,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL    *_force,                   // ( 1:MaxAtom, 1:3, 1:ncell_all, 1:nthread )
    double  *_ene_virial,
    const REALI *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_ene,           // ( 1:6*cutoff_int )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:ncel_max, 1:ncel_max)
    int   MaxAtom,
    int   MaxAtomCls,
    int   num_atom_cls,
    int   ncel_local,
    int   ncel_bound,
    int   ncel_max,
    int   cutoff_int,
    int   univ_maxcell,
    int   univ_maxcell1,
    int   univ_ncell_nonzero,
    int   univ_ncell_near,
    int   univ_update,
    int   univ_mask2_size,
    int   univ_natom_max,
    int   check_virial,
    REAL  density,
    REAL  cutoff2,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{

#ifdef DEBUG
    printf( "[%s,%s,%d]\n", __FILE__, __func__, __LINE__ );
    printf( "MaxAtom = %d\n", MaxAtom );
    printf( "MaxAtomCls = %d\n", MaxAtomCls );
    printf( "num_atom_cls = %d\n", num_atom_cls );
    printf( "ncel_local = %d\n", ncel_local );
    printf( "ncel_bound = %d\n", ncel_bound );
    printf( "ncel_max   = %d\n", ncel_max );
    printf( "cutoff_int = %d\n", cutoff_int );
    printf( "univ_maxcell = %d\n", univ_maxcell );
    printf( "univ_maxcell1 = %d\n", univ_maxcell1 );
    printf( "univ_ncell_nonzero = %d\n", univ_ncell_nonzero );
    printf( "univ_update = %d\n", univ_update );
    printf( "check_virial = %d\n", check_virial );
    printf( "cutoff2 = %lf\n", cutoff2 );
#endif

    static int first = 1;

    // memory allocation
    if ( first ) {
        gpu_init_buffer(
            _coord_pbc,
            _force,
            _ene_virial,
            _coord,
            _trans1,
            _cell_move,
            _charge,
            _atmcls,
            _natom,
            _nonb_lj12,
            _nonb_lj6,
            _nonb_lj6_factor,
            _table_ene,
            _table_grad,
            _univ_cell_pairlist1,
            _univ_mask2,
            _univ_ix_natom,
            _univ_ix_list,
            _univ_iy_natom,
            _univ_iy_list,
            _virial_check,
            MaxAtom,
            MaxAtomCls,
            num_atom_cls,
            ncel_local,
            ncel_bound,
            ncel_max,
            cutoff_int,
            univ_maxcell,
            univ_maxcell1,
            univ_ncell_near,
            univ_mask2_size
            );
    }

    // send data to GPU
    gpu_memcpy_h2d_energy(
        _coord_pbc,
        _force,
        _ene_virial,
        _coord,
        _trans1,
        _cell_move,
        _charge,
        _atmcls,
        _natom,
        _nonb_lj12,
        _nonb_lj6,
        _nonb_lj6_factor,
        _table_ene,
        _table_grad,
        _univ_cell_pairlist1,
        _univ_mask2,
        _univ_ix_natom,
        _univ_ix_list,
        _univ_iy_natom,
        _univ_iy_list,
        _univ_ij_sort_list,
        _virial_check,
        MaxAtom,
        univ_maxcell,
        univ_ncell_near,
        univ_mask2_size,
        univ_update,
        check_virial,
        first
        );

    first = 0;

    // return; // debug

    /* */
    dim3  def_dim3(1,1,1);
    dim3  num_block;
    dim3  num_thread;
    int  index_start, index_end;

    /* */
    num_block = def_dim3;
    num_block.x = ncel_local + ncel_bound;
    num_thread = def_dim3;
    num_thread.x = 256;

    /* coordinates in periodic boundary condition is written */
    kern_compute_force_nonbond_table_linear_univ__coord_pbc<<< num_block, num_thread, 0, stream[0] >>>(
        dev_coord_pbc,
        dev_coord,
        dev_trans1,
        dev_natom,
        MaxAtom,
        ncel_local,
        ncel_bound,
        ncel_max );

    /* energy/force calculation within a cell */
    index_start = 1;
    index_end  = ncel_local;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_end - index_start + 1), num_thread.y);

    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[1], 0 ) );
    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[2], 0 ) ); /* univ_mask2 */

    // CUDA_CALL( cudaStreamSynchronize( stream[1] ) );

    /* energy calculation */
    kern_compute_energy_nonbond_table_ljpme_univ__energyforce_intra_cell<<< num_block, num_thread, 0, stream[0] >>>(
        dev_coord_pbc,
        dev_force,
        dev_ene_virial,
        dev_ene_viri_mid,
        dev_cell_move,
        dev_charge,
        dev_atmcls,
        dev_natom,
        dev_nonb_lj12,
        dev_nonb_lj6,
        dev_nonb_lj6_factor,
        dev_table_ene,
        dev_table_grad,
        dev_univ_cell_pairlist1,
        dev_univ_mask2,
        dev_univ_ix_natom,
        dev_univ_ix_list,
        dev_univ_iy_natom,
        dev_univ_iy_list,
        dev_univ_ij_sort_list,
        MaxAtom,
        MaxAtomCls,
        num_atom_cls,
        ncel_local,
        ncel_bound,
        ncel_max,
        cutoff_int,
        univ_maxcell,
        univ_maxcell1,
        univ_mask2_size,
        univ_natom_max,
        index_start, index_end,
        density,
        cutoff2 );

    /* energy/force calculation between different near cells */
    index_start = ncel_local + 1;
    index_end = univ_ncell_nonzero;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;  /* do not change */
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_end - index_start + 1), num_thread.y);

    int max_iy_natom = (SMEM_SIZE/sizeof(REAL)/3) / (num_thread.y*NUM_CTA__ENERGYFORCE_INTER_CELL);
    max_iy_natom -= max_iy_natom % 4;

    size_t  smem_size = num_thread.y * max_iy_natom * 3 * sizeof(REAL);

    kern_compute_energy_nonbond_table_ljpme_univ__energyforce_inter_cell<<< num_block, num_thread, smem_size, stream[0] >>>(
	dev_coord_pbc,
	dev_force,
	dev_ene_virial,
	dev_ene_viri_mid,
	dev_cell_move,
	dev_charge,
	dev_atmcls,
	dev_natom,
	dev_nonb_lj12,
	dev_nonb_lj6,
	dev_nonb_lj6_factor,
	dev_table_ene,
	dev_table_grad,
	dev_univ_cell_pairlist1,
        dev_univ_mask2,
	dev_univ_ix_natom,
	dev_univ_ix_list,
	dev_univ_iy_natom,
	dev_univ_iy_list,
	dev_univ_ij_sort_list,
	dev_virial_check,
	MaxAtom,
	MaxAtomCls,
	num_atom_cls,
	ncel_local,
	ncel_bound,
	ncel_max,
	cutoff_int,
	univ_maxcell,
	univ_maxcell1,
        univ_ncell_near,
        univ_mask2_size,
        univ_natom_max,
	index_start, index_end, max_iy_natom,
	density,
	cutoff2,
        system_x, system_y, system_z );

    /* */
    num_block = def_dim3;
    num_thread = def_dim3;
    num_thread.x = 32;
    num_thread.y = 32;
    kern_compute_energy_nonbond_table_linear_univ_energy_sum<<< num_block, num_thread, 0, stream[0] >>>(
	dev_ene_virial,
	dev_ene_viri_mid,
	ncel_local,
	ncel_max,
	univ_maxcell,
	univ_ncell_nonzero);

}

/*
 * JJ : energy calculation without lookup table of vdw
 */
void gpu_launch_compute_energy_nonbond_notable_univ(
    REAL    *_coord_pbc,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL    *_force,                   // ( 1:MaxAtom, 1:3, 1:ncell_all, 1:nthread )
    double  *_ene_virial,
    const REALI *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_ene,           // ( 1:6*cutoff_int )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:ncel_max, 1:ncel_max)
    int   MaxAtom,
    int   MaxAtomCls,
    int   num_atom_cls,
    int   ncel_local,
    int   ncel_bound,
    int   ncel_max,
    int   cutoff_int,
    int   univ_maxcell,
    int   univ_maxcell1,
    int   univ_ncell_nonzero,
    int   univ_ncell_near,
    int   univ_update,
    int   univ_mask2_size,
    int   univ_natom_max,
    int   check_virial,
    REAL  density,
    REAL  cutoff2,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{

#ifdef DEBUG
    printf( "[%s,%s,%d]\n", __FILE__, __func__, __LINE__ );
    printf( "MaxAtom = %d\n", MaxAtom );
    printf( "MaxAtomCls = %d\n", MaxAtomCls );
    printf( "num_atom_cls = %d\n", num_atom_cls );
    printf( "ncel_local = %d\n", ncel_local );
    printf( "ncel_bound = %d\n", ncel_bound );
    printf( "ncel_max   = %d\n", ncel_max );
    printf( "cutoff_int = %d\n", cutoff_int );
    printf( "univ_maxcell = %d\n", univ_maxcell );
    printf( "univ_maxcell1 = %d\n", univ_maxcell1 );
    printf( "univ_ncell_nonzero = %d\n", univ_ncell_nonzero );
    printf( "univ_update = %d\n", univ_update );
    printf( "check_virial = %d\n", check_virial );
    printf( "cutoff2 = %lf\n", cutoff2 );
#endif

    static int first = 1;

    // memory allocation
    if ( first ) {
        gpu_init_buffer(
            _coord_pbc,
            _force,
            _ene_virial,
            _coord,
            _trans1,
            _cell_move,
            _charge,
            _atmcls,
            _natom,
            _nonb_lj12,
            _nonb_lj6,
            _nonb_lj6_factor,
            _table_ene,
            _table_grad,
            _univ_cell_pairlist1,
            _univ_mask2,
            _univ_ix_natom,
            _univ_ix_list,
            _univ_iy_natom,
            _univ_iy_list,
            _virial_check,
            MaxAtom,
            MaxAtomCls,
            num_atom_cls,
            ncel_local,
            ncel_bound,
            ncel_max,
            cutoff_int,
            univ_maxcell,
            univ_maxcell1,
            univ_ncell_near,
            univ_mask2_size
            );
    }

    // send data to GPU
    gpu_memcpy_h2d_energy(
        _coord_pbc,
        _force,
        _ene_virial,
        _coord,
        _trans1,
        _cell_move,
        _charge,
        _atmcls,
        _natom,
        _nonb_lj12,
        _nonb_lj6,
        _nonb_lj6_factor,
        _table_ene,
        _table_grad,
        _univ_cell_pairlist1,
        _univ_mask2,
        _univ_ix_natom,
        _univ_ix_list,
        _univ_iy_natom,
        _univ_iy_list,
        _univ_ij_sort_list,
        _virial_check,
        MaxAtom,
        univ_maxcell,
        univ_ncell_near,
        univ_mask2_size,
        univ_update,
        check_virial,
        first
        );

    first = 0;

    // return; // debug

    /* */
    dim3  def_dim3(1,1,1);
    dim3  num_block;
    dim3  num_thread;
    int  index_start, index_end;

    /* */
    num_block = def_dim3;
    num_block.x = ncel_local + ncel_bound;
    num_thread = def_dim3;
    num_thread.x = 256;

    /* coordinates in periodic boundary condition is written */
    kern_compute_force_nonbond_table_linear_univ__coord_pbc<<< num_block, num_thread, 0, stream[0] >>>(
        dev_coord_pbc,
        dev_coord,
        dev_trans1,
        dev_natom,
        MaxAtom,
        ncel_local,
        ncel_bound,
        ncel_max );

    /* energy/force calculation within a cell */
    index_start = 1;
    index_end  = ncel_local;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_end - index_start + 1), num_thread.y);

    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[1], 0 ) ); /* univ_mask2 */
    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[2], 0 ) ); /* univ_mask2 */

    // CUDA_CALL( cudaStreamSynchronize( stream[1] ) );

    /* energy calculation */
    kern_compute_energy_nonbond_notable_univ__energyforce_intra_cell<<< num_block, num_thread, 0, stream[0] >>>(
        dev_coord_pbc,
        dev_force,
        dev_ene_virial,
        dev_ene_viri_mid,
        dev_cell_move,
        dev_charge,
        dev_atmcls,
        dev_natom,
        dev_nonb_lj12,
        dev_nonb_lj6,
        dev_table_ene,
        dev_table_grad,
        dev_univ_cell_pairlist1,
        dev_univ_mask2,
        dev_univ_ix_natom,
        dev_univ_ix_list,
        dev_univ_iy_natom,
        dev_univ_iy_list,
        dev_univ_ij_sort_list,
        MaxAtom,
        MaxAtomCls,
        num_atom_cls,
        ncel_local,
        ncel_bound,
        ncel_max,
        cutoff_int,
        univ_maxcell,
        univ_maxcell1,
        univ_mask2_size,
        univ_natom_max,
        index_start, index_end,
        density,
        cutoff2 );

    /* energy/force calculation between different near cells */
    index_start = ncel_local + 1;
    index_end = univ_ncell_nonzero;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;  /* do not change */
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_end - index_start + 1), num_thread.y);

    int max_iy_natom = (SMEM_SIZE/sizeof(REAL)/3) / (num_thread.y*NUM_CTA__ENERGYFORCE_INTER_CELL);
    max_iy_natom -= max_iy_natom % 4;

    size_t  smem_size = num_thread.y * max_iy_natom * 3 * sizeof(REAL);

    kern_compute_energy_nonbond_notable_univ__energyforce_inter_cell<<< num_block, num_thread, smem_size, stream[0] >>>(
	dev_coord_pbc,
	dev_force,
	dev_ene_virial,
	dev_ene_viri_mid,
	dev_cell_move,
	dev_charge,
	dev_atmcls,
	dev_natom,
	dev_nonb_lj12,
	dev_nonb_lj6,
	dev_table_ene,
	dev_table_grad,
	dev_univ_cell_pairlist1,
        dev_univ_mask2,
	dev_univ_ix_natom,
	dev_univ_ix_list,
	dev_univ_iy_natom,
	dev_univ_iy_list,
	dev_univ_ij_sort_list,
	dev_virial_check,
	MaxAtom,
	MaxAtomCls,
	num_atom_cls,
	ncel_local,
	ncel_bound,
	ncel_max,
	cutoff_int,
	univ_maxcell,
	univ_maxcell1,
        univ_ncell_near,
        univ_mask2_size,
        univ_natom_max,
	index_start, index_end, max_iy_natom,
	density,
	cutoff2,
        system_x, system_y, system_z );
    /* */
    num_block = def_dim3;
    num_thread = def_dim3;
    num_thread.x = 32;
    num_thread.y = 32;
    kern_compute_energy_nonbond_table_linear_univ_energy_sum<<< num_block, num_thread, 0, stream[0] >>>(
	dev_ene_virial,
	dev_ene_viri_mid,
	ncel_local,
	ncel_max,
	univ_maxcell,
	univ_ncell_nonzero);

}


void gpu_launch_compute_force_nonbond_table_linear_univ(
    REAL    *_coord_pbc,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL    *_force,                   // ( 1:MaxAtom, 1:3, 1:ncell_all, 1:nthread )
    double  *_ene_virial,              // ( 5 )
    const REALI *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near )
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:univ_maxcell1 )
    int   MaxAtom,
    int   MaxAtomCls,
    int   num_atom_cls,
    int   ncel_local,
    int   ncel_bound,
    int   ncel_max,
    int   cutoff_int,
    int   univ_maxcell,
    int   univ_maxcell1,
    int	  univ_ncell_nonzero,
    int	  univ_ncell_near,
    int	  univ_update,
    int   univ_mask2_size,
    int   univ_natom_max,
    int   check_virial,
    int   cpu_calc,
    REAL  density,
    REAL  cutoff2,
    REAL  pairlistdist2,
    int   univ_gpu_start,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{
#ifdef DEBUG
    printf( "[%s,%s,%d]\n", __FILE__, __func__, __LINE__ );
    printf( "MaxAtom = %d\n", MaxAtom );
    printf( "MaxAtomCls = %d\n", MaxAtomCls );
    printf( "num_atom_cls = %d\n", num_atom_cls );
    printf( "ncel_local = %d\n", ncel_local );
    printf( "ncel_bound = %d\n", ncel_bound );
    printf( "ncel_max   = %d\n", ncel_max );
    printf( "cutoff_int = %d\n", cutoff_int );
    printf( "univ_maxcell = %d\n", univ_maxcell );
    printf( "univ_maxcell1 = %d\n", univ_maxcell1 );
    printf( "univ_ncell_nonzero = %d\n", univ_ncell_nonzero );
    printf( "univ_update = %d\n", univ_update );
    printf( "check_virial = %d\n", check_virial );
    printf( "cutoff2 = %lf\n", cutoff2 );
    printf( "pairlistdist2 = %lf\n", pairlistdist2 );
#endif

    gpu_memcpy_h2d_force(
	_coord_pbc,
	_force,
        _ene_virial,
	_coord,
	_trans1,
	_cell_move,
	_charge,
	_atmcls,
	_natom,
	_nonb_lj12,
	_nonb_lj6,
	_nonb_lj6_factor,
	_table_grad,
	_univ_cell_pairlist1,
        _univ_mask2,
	_univ_ix_natom,
	_univ_ix_list,
	_univ_iy_natom,
	_univ_iy_list,
	_univ_ij_sort_list,
        _virial_check,
	MaxAtom,
	univ_maxcell,
        univ_ncell_near,
        univ_mask2_size,
	univ_update,
	check_virial,
	stream[0]
	);

    // return; // debug

    /* */
    dim3  def_dim3(1,1,1);
    dim3  num_block;
    dim3  num_thread;
    int  index_s, index_e;

    /* */
    if ( flag_build_pairlist_on_GPU == 0 || univ_update == 0 ) {
	num_block = def_dim3;
	num_block.x = ncel_local + ncel_bound;
	num_thread = def_dim3;
	num_thread.x = 256;
	kern_compute_force_nonbond_table_linear_univ__coord_pbc<<< num_block, num_thread, 0, stream[0] >>>(
	    dev_coord_pbc,
	    dev_coord,
	    dev_trans1,
	    dev_natom,
	    MaxAtom,
	    ncel_local,
	    ncel_bound,
	    ncel_max );
    }

    /* force calculation within a cell */
//  if ( cpu_calc == 0) index_s = 1;
//  if ( cpu_calc != 0) index_s = ncel_local + 1;
    index_s = 1;
    index_e = ncel_local;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_e - index_s + 1), num_thread.y);

    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[1], 0 ) ); /* univ_mask2 */
    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[2], 0 ) ); /* univ_mask2 */
    // CUDA_CALL( cudaStreamSynchronize( stream[1] ) );


    kern_compute_force_nonbond_table_linear_univ__force_intra_cell<<< num_block, num_thread, 0, stream[0] >>>(
	dev_coord_pbc,
	dev_force,
        dev_ene_virial,
        dev_ene_viri_mid,
	dev_cell_move,
	dev_charge,
	dev_atmcls,
	dev_natom,
	dev_nonb_lj12,
	dev_nonb_lj6,
	dev_table_grad,
	dev_univ_cell_pairlist1,
        dev_univ_mask2,
	dev_univ_ix_natom,
	dev_univ_ix_list,
        dev_univ_iy_natom,
	dev_univ_iy_list,
	dev_univ_ij_sort_list,
	MaxAtom,
	MaxAtomCls,
	num_atom_cls,
	ncel_local,
	ncel_bound,
	ncel_max,
	cutoff_int,
	univ_maxcell,
	univ_maxcell1,
        univ_mask2_size,
        univ_natom_max,
	index_s, index_e,
        check_virial,
	density,
	cutoff2,
	pairlistdist2 );

    /* force calculation between different cells */
//  index_s = ncel_local + 1;
//  index_e = univ_ncell_nonzero;
    index_s = univ_gpu_start + 1;
    index_e = univ_ncell_nonzero;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;  /* do not change */
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_e - index_s + 1), num_thread.y);

    // int max_iy_natom = (SMEM_SIZE/sizeof(REAL)/3) / (num_thread.y*NUM_CTA__FORCE_INTER_CELL);
    int max_iy_natom = SMEM_SIZE / (num_thread.y*NUM_CTA__FORCE_INTER_CELL) / (sizeof(REAL)*6);
    max_iy_natom -= max_iy_natom % 4;
    // if ( max_iy_natom > 32 ) max_iy_natom = 32;

    // size_t  smem_size = num_thread.y * max_iy_natom * 3 * sizeof(REAL);
    size_t smem_size = num_thread.y * sizeof(REAL) * (max_iy_natom*6);
    // printf( "max_iy_natom: %d, smem_size: %d\n", max_iy_natom, smem_size ); /* debug */

    kern_compute_force_nonbond_table_linear_univ__force_inter_cell<<< num_block, num_thread, smem_size, stream[0] >>>(
	dev_coord_pbc,
	dev_force,
	dev_ene_virial,
	dev_ene_viri_mid,
	dev_cell_move,
	dev_charge,
	dev_atmcls,
	dev_natom,
	dev_nonb_lj12,
	dev_nonb_lj6,
	dev_table_grad,
	dev_univ_cell_pairlist1,
        dev_univ_mask2,
	dev_univ_ix_natom,
	dev_univ_ix_list,
	dev_univ_iy_natom,
	dev_univ_iy_list,
	dev_univ_ij_sort_list,
	dev_virial_check,
	MaxAtom,
	MaxAtomCls,
	num_atom_cls,
	ncel_local,
	ncel_bound,
	ncel_max,
	cutoff_int,
	univ_maxcell,
	univ_maxcell1,
        univ_ncell_near,
        univ_mask2_size,
        univ_natom_max,
	index_s, index_e, max_iy_natom,
	check_virial,
	density,
	cutoff2,
	pairlistdist2,
        system_x, system_y, system_z );

    if (check_virial != 0) {
       num_block = def_dim3;
       num_thread = def_dim3;
       num_thread.x = 32;
       num_thread.y = 32;

       kern_compute_force_nonbond_table_linear_univ_sum<<< num_block, num_thread, 0, stream[0] >>>(
               dev_ene_virial,
               dev_ene_viri_mid,
               ncel_local,
               ncel_max,
               univ_maxcell,
               univ_gpu_start,
               univ_ncell_nonzero);
    }

}

void gpu_launch_compute_force_nonbond_table_ljpme_univ(
    REAL    *_coord_pbc,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL    *_force,                   // ( 1:MaxAtom, 1:3, 1:ncell_all, 1:nthread )
    double  *_ene_virial,              // ( 5 )
    const REALI *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near )
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:univ_maxcell1 )
    int   MaxAtom,
    int   MaxAtomCls,
    int   num_atom_cls,
    int   ncel_local,
    int   ncel_bound,
    int   ncel_max,
    int   cutoff_int,
    int   univ_maxcell,
    int   univ_maxcell1,
    int	  univ_ncell_nonzero,
    int	  univ_ncell_near,
    int	  univ_update,
    int   univ_mask2_size,
    int   univ_natom_max,
    int   check_virial,
    int   cpu_calc,
    REAL  density,
    REAL  cutoff2,
    REAL  pairlistdist2,
    int   univ_gpu_start,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{
#ifdef DEBUG
    printf( "[%s,%s,%d]\n", __FILE__, __func__, __LINE__ );
    printf( "MaxAtom = %d\n", MaxAtom );
    printf( "MaxAtomCls = %d\n", MaxAtomCls );
    printf( "num_atom_cls = %d\n", num_atom_cls );
    printf( "ncel_local = %d\n", ncel_local );
    printf( "ncel_bound = %d\n", ncel_bound );
    printf( "ncel_max   = %d\n", ncel_max );
    printf( "cutoff_int = %d\n", cutoff_int );
    printf( "univ_maxcell = %d\n", univ_maxcell );
    printf( "univ_maxcell1 = %d\n", univ_maxcell1 );
    printf( "univ_ncell_nonzero = %d\n", univ_ncell_nonzero );
    printf( "univ_update = %d\n", univ_update );
    printf( "check_virial = %d\n", check_virial );
    printf( "cutoff2 = %lf\n", cutoff2 );
    printf( "pairlistdist2 = %lf\n", pairlistdist2 );
#endif

    gpu_memcpy_h2d_force(
	_coord_pbc,
	_force,
        _ene_virial,
	_coord,
	_trans1,
	_cell_move,
	_charge,
	_atmcls,
	_natom,
	_nonb_lj12,
	_nonb_lj6,
	_nonb_lj6_factor,
	_table_grad,
	_univ_cell_pairlist1,
        _univ_mask2,
	_univ_ix_natom,
	_univ_ix_list,
	_univ_iy_natom,
	_univ_iy_list,
	_univ_ij_sort_list,
        _virial_check,
	MaxAtom,
	univ_maxcell,
        univ_ncell_near,
        univ_mask2_size,
	univ_update,
	check_virial,
	stream[0]
	);

    // return; // debug

    /* */
    dim3  def_dim3(1,1,1);
    dim3  num_block;
    dim3  num_thread;
    int  index_s, index_e;

    /* */
    if ( flag_build_pairlist_on_GPU == 0 || univ_update == 0 ) {
	num_block = def_dim3;
	num_block.x = ncel_local + ncel_bound;
	num_thread = def_dim3;
	num_thread.x = 256;
	kern_compute_force_nonbond_table_linear_univ__coord_pbc<<< num_block, num_thread, 0, stream[0] >>>(
	    dev_coord_pbc,
	    dev_coord,
	    dev_trans1,
	    dev_natom,
	    MaxAtom,
	    ncel_local,
	    ncel_bound,
	    ncel_max );
    }

    /* force calculation within a cell */
//  if ( cpu_calc == 0) index_s = 1;
//  if ( cpu_calc != 0) index_s = ncel_local + 1;
    index_s = 1;
    index_e = ncel_local;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_e - index_s + 1), num_thread.y);

    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[1], 0 ) ); /* univ_mask2 */
    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[2], 0 ) ); /* univ_mask2 */
    // CUDA_CALL( cudaStreamSynchronize( stream[1] ) );

    kern_compute_force_nonbond_table_ljpme_univ__force_intra_cell<<< num_block, num_thread, 0, stream[0] >>>(
	dev_coord_pbc,
	dev_force,
        dev_ene_virial,
        dev_ene_viri_mid,
	dev_cell_move,
	dev_charge,
	dev_atmcls,
	dev_natom,
	dev_nonb_lj12,
	dev_nonb_lj6,
	dev_nonb_lj6_factor,
	dev_table_grad,
	dev_univ_cell_pairlist1,
        dev_univ_mask2,
	dev_univ_ix_natom,
	dev_univ_ix_list,
        dev_univ_iy_natom,
	dev_univ_iy_list,
	dev_univ_ij_sort_list,
	MaxAtom,
	MaxAtomCls,
	num_atom_cls,
	ncel_local,
	ncel_bound,
	ncel_max,
	cutoff_int,
	univ_maxcell,
	univ_maxcell1,
        univ_mask2_size,
        univ_natom_max,
	index_s, index_e,
        check_virial,
	density,
	cutoff2,
	pairlistdist2 );

    /* force calculation between different cells */
//  index_s = ncel_local + 1;
//  index_e = univ_ncell_nonzero;
    index_s = univ_gpu_start + 1;
    index_e = univ_ncell_nonzero;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;  /* do not change */
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_e - index_s + 1), num_thread.y);

    // int max_iy_natom = (SMEM_SIZE/sizeof(REAL)/3) / (num_thread.y*NUM_CTA__FORCE_INTER_CELL);
    int max_iy_natom = SMEM_SIZE / (num_thread.y*NUM_CTA__FORCE_INTER_CELL) / (sizeof(REAL)*6);
    max_iy_natom -= max_iy_natom % 4;
    // if ( max_iy_natom > 32 ) max_iy_natom = 32;

    // size_t  smem_size = num_thread.y * max_iy_natom * 3 * sizeof(REAL);
    size_t smem_size = num_thread.y * sizeof(REAL) * (max_iy_natom*6);
    // printf( "max_iy_natom: %d, smem_size: %d\n", max_iy_natom, smem_size ); /* debug */

    kern_compute_force_nonbond_table_ljpme_univ__force_inter_cell<<< num_block, num_thread, smem_size, stream[0] >>>(
	dev_coord_pbc,
	dev_force,
	dev_ene_virial,
	dev_ene_viri_mid,
	dev_cell_move,
	dev_charge,
	dev_atmcls,
	dev_natom,
	dev_nonb_lj12,
	dev_nonb_lj6,
	dev_nonb_lj6_factor,
	dev_table_grad,
	dev_univ_cell_pairlist1,
        dev_univ_mask2,
	dev_univ_ix_natom,
	dev_univ_ix_list,
	dev_univ_iy_natom,
	dev_univ_iy_list,
	dev_univ_ij_sort_list,
	dev_virial_check,
	MaxAtom,
	MaxAtomCls,
	num_atom_cls,
	ncel_local,
	ncel_bound,
	ncel_max,
	cutoff_int,
	univ_maxcell,
	univ_maxcell1,
        univ_ncell_near,
        univ_mask2_size,
        univ_natom_max,
	index_s, index_e, max_iy_natom,
	check_virial,
	density,
	cutoff2,
	pairlistdist2,
        system_x, system_y, system_z );

    if (check_virial != 0) {
       num_block = def_dim3;
       num_thread = def_dim3;
       num_thread.x = 32;
       num_thread.y = 32;

       kern_compute_force_nonbond_table_linear_univ_sum<<< num_block, num_thread, 0, stream[0] >>>(
               dev_ene_virial,
               dev_ene_viri_mid,
               ncel_local,
               ncel_max,
               univ_maxcell,
               univ_gpu_start,
               univ_ncell_nonzero);
    }

}

void gpu_launch_compute_force_nonbond_notable_univ(
    REAL    *_coord_pbc,               // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL    *_force,                   // ( 1:MaxAtom, 1:3, 1:ncell_all, 1:nthread )
    double  *_ene_virial,              // ( 5 )
    const REALI *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near )
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:univ_maxcell1 )
    int   MaxAtom,
    int   MaxAtomCls,
    int   num_atom_cls,
    int   ncel_local,
    int   ncel_bound,
    int   ncel_max,
    int   cutoff_int,
    int   univ_maxcell,
    int   univ_maxcell1,
    int	  univ_ncell_nonzero,
    int	  univ_ncell_near,
    int	  univ_update,
    int   univ_mask2_size,
    int   univ_natom_max,
    int   check_virial,
    int   cpu_calc,
    REAL  density,
    REAL  cutoff2,
    REAL  pairlistdist2,
    int   univ_gpu_start,
    REAL  system_x,
    REAL  system_y,
    REAL  system_z
    )
{
#ifdef DEBUG
    printf( "[%s,%s,%d]\n", __FILE__, __func__, __LINE__ );
    printf( "MaxAtom = %d\n", MaxAtom );
    printf( "MaxAtomCls = %d\n", MaxAtomCls );
    printf( "num_atom_cls = %d\n", num_atom_cls );
    printf( "ncel_local = %d\n", ncel_local );
    printf( "ncel_bound = %d\n", ncel_bound );
    printf( "ncel_max   = %d\n", ncel_max );
    printf( "cutoff_int = %d\n", cutoff_int );
    printf( "univ_maxcell = %d\n", univ_maxcell );
    printf( "univ_maxcell1 = %d\n", univ_maxcell1 );
    printf( "univ_ncell_nonzero = %d\n", univ_ncell_nonzero );
    printf( "univ_update = %d\n", univ_update );
    printf( "check_virial = %d\n", check_virial );
    printf( "cutoff2 = %lf\n", cutoff2 );
    printf( "pairlistdist2 = %lf\n", pairlistdist2 );
#endif

    gpu_memcpy_h2d_force(
	_coord_pbc,
	_force,
        _ene_virial,
	_coord,
	_trans1,
	_cell_move,
	_charge,
	_atmcls,
	_natom,
	_nonb_lj12,
	_nonb_lj6,
	_nonb_lj6_factor,
	_table_grad,
	_univ_cell_pairlist1,
        _univ_mask2,
	_univ_ix_natom,
	_univ_ix_list,
	_univ_iy_natom,
	_univ_iy_list,
	_univ_ij_sort_list,
        _virial_check,
	MaxAtom,
	univ_maxcell,
        univ_ncell_near,
        univ_mask2_size,
	univ_update,
	check_virial,
	stream[0]
	);

    // return; // debug

    /* */
    dim3  def_dim3(1,1,1);
    dim3  num_block;
    dim3  num_thread;
    int  index_s, index_e;

    /* */
    if ( flag_build_pairlist_on_GPU == 0 || univ_update == 0 ) {
	num_block = def_dim3;
	num_block.x = ncel_local + ncel_bound;
	num_thread = def_dim3;
	num_thread.x = 256;
	kern_compute_force_nonbond_table_linear_univ__coord_pbc<<< num_block, num_thread, 0, stream[0] >>>(
	    dev_coord_pbc,
	    dev_coord,
	    dev_trans1,
	    dev_natom,
	    MaxAtom,
	    ncel_local,
	    ncel_bound,
	    ncel_max );
    }

    /* force calculation within a cell */
//  if ( cpu_calc == 0) index_s = 1;
//  if ( cpu_calc != 0) index_s = ncel_local + 1;
    index_s = 1;
    index_e = ncel_local;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_e - index_s + 1), num_thread.y);

    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[1], 0 ) ); /* univ_mask2 */
    CUDA_CALL( cudaStreamWaitEvent( stream[0], event[2], 0 ) ); /* univ_mask2 */
    // CUDA_CALL( cudaStreamSynchronize( stream[1] ) );


    kern_compute_force_nonbond_notable_univ__force_intra_cell<<< num_block, num_thread, 0, stream[0] >>>(
	dev_coord_pbc,
	dev_force,
        dev_ene_virial,
        dev_ene_viri_mid,
	dev_cell_move,
	dev_charge,
	dev_atmcls,
	dev_natom,
	dev_nonb_lj12,
	dev_nonb_lj6,
	dev_table_grad,
	dev_univ_cell_pairlist1,
        dev_univ_mask2,
	dev_univ_ix_natom,
	dev_univ_ix_list,
        dev_univ_iy_natom,
	dev_univ_iy_list,
	dev_univ_ij_sort_list,
	MaxAtom,
	MaxAtomCls,
	num_atom_cls,
	ncel_local,
	ncel_bound,
	ncel_max,
	cutoff_int,
	univ_maxcell,
	univ_maxcell1,
        univ_mask2_size,
        univ_natom_max,
	index_s, index_e,
        check_virial,
	density,
	cutoff2,
	pairlistdist2 );

    /* force calculation between different cells */
//  index_s = ncel_local + 1;
//  index_e = univ_ncell_nonzero;
    index_s = univ_gpu_start + 1;
    index_e = univ_ncell_nonzero;
    num_thread = def_dim3;
    num_thread.x = 32; num_thread.y = 4;  /* do not change */
    assert( num_thread.x == 32 );
    assert( num_thread.y ==  4 );
    num_block = def_dim3;
    num_block.x = DIVCEIL((index_e - index_s + 1), num_thread.y);

    // int max_iy_natom = (SMEM_SIZE/sizeof(REAL)/3) / (num_thread.y*NUM_CTA__FORCE_INTER_CELL);
    int max_iy_natom = SMEM_SIZE / (num_thread.y*NUM_CTA__FORCE_INTER_CELL) / (sizeof(REAL)*6);
    max_iy_natom -= max_iy_natom % 4;
    // if ( max_iy_natom > 32 ) max_iy_natom = 32;

    // size_t  smem_size = num_thread.y * max_iy_natom * 3 * sizeof(REAL);
    size_t smem_size = num_thread.y * sizeof(REAL) * (max_iy_natom*6);
    // printf( "max_iy_natom: %d, smem_size: %d\n", max_iy_natom, smem_size ); /* debug */

    kern_compute_force_nonbond_notable_univ__force_inter_cell<<< num_block, num_thread, smem_size, stream[0] >>>(
	dev_coord_pbc,
	dev_force,
	dev_ene_virial,
	dev_ene_viri_mid,
	dev_cell_move,
	dev_charge,
	dev_atmcls,
	dev_natom,
	dev_nonb_lj12,
	dev_nonb_lj6,
	dev_table_grad,
	dev_univ_cell_pairlist1,
        dev_univ_mask2,
	dev_univ_ix_natom,
	dev_univ_ix_list,
	dev_univ_iy_natom,
	dev_univ_iy_list,
	dev_univ_ij_sort_list,
	dev_virial_check,
	MaxAtom,
	MaxAtomCls,
	num_atom_cls,
	ncel_local,
	ncel_bound,
	ncel_max,
	cutoff_int,
	univ_maxcell,
	univ_maxcell1,
        univ_ncell_near,
        univ_mask2_size,
        univ_natom_max,
	index_s, index_e, max_iy_natom,
	check_virial,
	density,
	cutoff2,
	pairlistdist2,
        system_x, system_y, system_z );

    if (check_virial != 0) {
       num_block = def_dim3;
       num_thread = def_dim3;
       num_thread.x = 32;
       num_thread.y = 32;

       kern_compute_force_nonbond_table_linear_univ_sum<<< num_block, num_thread, 0, stream[0] >>>(
               dev_ene_virial,
               dev_ene_viri_mid,
               ncel_local,
               ncel_max,
               univ_maxcell,
               univ_gpu_start,
               univ_ncell_nonzero);
    }

}

/*
 * wapper function called from fortran subroutine
 */
extern "C"
void gpu_launch_compute_energy_nonbond_table_linear_univ_(
    REAL        *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        *_force,               // ( 1:MaxAtom, 1:3, 1:ncell_all )
    double      *_ene_virial,
    const REALI *_coord,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_ene,           // ( 1:6*cutoff_int )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:ncel_max, 1:ncel_max)
    int  *_MaxAtom,
    int  *_MaxAtomCls,
    int  *_num_atom_cls,
    int  *_ncel_local,
    int  *_ncel_bound,
    int  *_ncel_max,
    int  *_cutoff_int,
    int  *_univ_maxcell,
    int  *_univ_maxcell1,
    int  *_univ_ncell_nonzero,
    int  *_univ_ncell_near,
    int  *_univ_update,
    int  *_univ_mask2_size,
    int  *_univ_natom_max,
    int  *_maxcell,
    REAL *_density,
    REAL *_cutoff2,
    REAL *_system_x,
    REAL *_system_y,
    REAL *_system_z
    )
{
    gpu_init();

    gpu_launch_compute_energy_nonbond_table_linear_univ(
        _coord_pbc,
        _force,
        _ene_virial,
        _coord,
        _trans1,
        _cell_move,
        _charge,
        _atmcls,
        _natom,
        _nonb_lj12,
        _nonb_lj6,
        _nonb_lj6_factor,
        _table_ene,
        _table_grad,
        _univ_cell_pairlist1,
        _univ_mask2,
        _univ_ix_natom,
        _univ_ix_list,
        _univ_iy_natom,
        _univ_iy_list,
        _univ_ij_sort_list,
        _virial_check,
        *_MaxAtom,
        *_MaxAtomCls,
        *_num_atom_cls,
        *_ncel_local,
        *_ncel_bound,
        *_ncel_max,
        *_cutoff_int,
        *_univ_maxcell,
        *_univ_maxcell1,
        *_univ_ncell_nonzero,
        *_univ_ncell_near,
        *_univ_update,
        *_univ_mask2_size,
        *_univ_natom_max,
        *_maxcell,
        *_density,
        *_cutoff2,
        *_system_x, *_system_y, *_system_z
        );
}

extern "C"
void gpu_launch_compute_energy_nonbond_table_ljpme_univ_(
    REAL        *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        *_force,               // ( 1:MaxAtom, 1:3, 1:ncell_all )
    double      *_ene_virial,
    const REALI *_coord,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_ene,           // ( 1:6*cutoff_int )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:ncel_max, 1:ncel_max)
    int  *_MaxAtom,
    int  *_MaxAtomCls,
    int  *_num_atom_cls,
    int  *_ncel_local,
    int  *_ncel_bound,
    int  *_ncel_max,
    int  *_cutoff_int,
    int  *_univ_maxcell,
    int  *_univ_maxcell1,
    int  *_univ_ncell_nonzero,
    int  *_univ_ncell_near,
    int  *_univ_update,
    int  *_univ_mask2_size,
    int  *_univ_natom_max,
    int  *_maxcell,
    REAL *_density,
    REAL *_cutoff2,
    REAL *_system_x,
    REAL *_system_y,
    REAL *_system_z
    )
{
    gpu_init();

    gpu_launch_compute_energy_nonbond_table_ljpme_univ(
        _coord_pbc,
        _force,
        _ene_virial,
        _coord,
        _trans1,
        _cell_move,
        _charge,
        _atmcls,
        _natom,
        _nonb_lj12,
        _nonb_lj6,
        _nonb_lj6_factor,
        _table_ene,
        _table_grad,
        _univ_cell_pairlist1,
        _univ_mask2,
        _univ_ix_natom,
        _univ_ix_list,
        _univ_iy_natom,
        _univ_iy_list,
        _univ_ij_sort_list,
        _virial_check,
        *_MaxAtom,
        *_MaxAtomCls,
        *_num_atom_cls,
        *_ncel_local,
        *_ncel_bound,
        *_ncel_max,
        *_cutoff_int,
        *_univ_maxcell,
        *_univ_maxcell1,
        *_univ_ncell_nonzero,
        *_univ_ncell_near,
        *_univ_update,
        *_univ_mask2_size,
        *_univ_natom_max,
        *_maxcell,
        *_density,
        *_cutoff2,
        *_system_x, *_system_y, *_system_z
        );
}


extern "C"
void gpu_launch_compute_energy_nonbond_notable_univ_(
    REAL        *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        *_force,               // ( 1:MaxAtom, 1:3, 1:ncell_all )
    double      *_ene_virial,
    const REALI *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_ene,           // ( 1:6*cutoff_int )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:ncel_max, 1:ncel_max)
    int  *_MaxAtom,
    int  *_MaxAtomCls,
    int  *_num_atom_cls,
    int  *_ncel_local,
    int  *_ncel_bound,
    int  *_ncel_max,
    int  *_cutoff_int,
    int  *_univ_maxcell,
    int  *_univ_maxcell1,
    int  *_univ_ncell_nonzero,
    int  *_univ_ncell_near,
    int  *_univ_update,
    int  *_univ_mask2_size,
    int  *_univ_natom_max,
    int  *_maxcell,
    REAL *_density,
    REAL *_cutoff2,
    REAL *_system_x,
    REAL *_system_y,
    REAL *_system_z
    )
{
    gpu_init();

    gpu_launch_compute_energy_nonbond_notable_univ(
        _coord_pbc,
        _force,
        _ene_virial,
        _coord,
        _trans1,
        _cell_move,
        _charge,
        _atmcls,
        _natom,
        _nonb_lj12,
        _nonb_lj6,
        _nonb_lj6_factor,
        _table_ene,
        _table_grad,
        _univ_cell_pairlist1,
        _univ_mask2,
        _univ_ix_natom,
        _univ_ix_list,
        _univ_iy_natom,
        _univ_iy_list,
        _univ_ij_sort_list,
        _virial_check,
        *_MaxAtom,
        *_MaxAtomCls,
        *_num_atom_cls,
        *_ncel_local,
        *_ncel_bound,
        *_ncel_max,
        *_cutoff_int,
        *_univ_maxcell,
        *_univ_maxcell1,
        *_univ_ncell_nonzero,
        *_univ_ncell_near,
        *_univ_update,
        *_univ_mask2_size,
        *_univ_natom_max,
        *_maxcell,
        *_density,
        *_cutoff2,
        *_system_x, *_system_y, *_system_z
        );
}


/*
 * wapper function called from fortran subroutine
 */
extern "C"
void gpu_launch_compute_force_nonbond_table_linear_univ_(
    REAL        *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        *_force,               // ( 1:MaxAtom, 1:3, 1:ncell_all )
    double      *_ene_virial,          // ( 5 )
    const REALI *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:univ_maxcell1 )
    int  *_MaxAtom,
    int  *_MaxAtomCls,
    int  *_num_atom_cls,
    int  *_ncel_local,
    int  *_ncel_bound,
    int  *_ncel_max,
    int  *_cutoff_int,
    int  *_univ_maxcell,
    int  *_univ_maxcell1,
    int	 *_univ_ncell_nonzero,
    int	 *_univ_ncell_near,
    int	 *_univ_update,
    int  *_univ_mask2_size,
    int  *_univ_natom_max,
    int  *_check_virial,
    int  *_cpu_calc,
    REAL *_density,
    REAL *_cutoff2,
    REAL *_pairlistdist2,
    int  *_univ_gpu_start,
    REAL *_system_x,
    REAL *_system_y,
    REAL *_system_z
    )
{
    gpu_init();

    gpu_launch_compute_force_nonbond_table_linear_univ(
	_coord_pbc,
	_force,
        _ene_virial,
	_coord,
	_trans1,
	_cell_move,
	_charge,
	_atmcls,
	_natom,
	_nonb_lj12,
	_nonb_lj6,
	_nonb_lj6_factor,
	_table_grad,
	_univ_cell_pairlist1,
        _univ_mask2,
	_univ_ix_natom,
	_univ_ix_list,
        _univ_iy_natom,
	_univ_iy_list,
	_univ_ij_sort_list,
	_virial_check,
	*_MaxAtom,
	*_MaxAtomCls,
	*_num_atom_cls,
	*_ncel_local,
	*_ncel_bound,
	*_ncel_max,
	*_cutoff_int,
	*_univ_maxcell,
	*_univ_maxcell1,
	*_univ_ncell_nonzero,
	*_univ_ncell_near,
	*_univ_update,
        *_univ_mask2_size,
        *_univ_natom_max,
	*_check_virial,
        *_cpu_calc,
	*_density,
	*_cutoff2,
	*_pairlistdist2,
        *_univ_gpu_start,
        *_system_x, *_system_y, *_system_z
	);
}

extern "C"
void gpu_launch_compute_force_nonbond_table_ljpme_univ_(
    REAL        *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        *_force,               // ( 1:MaxAtom, 1:3, 1:ncell_all )
    double      *_ene_virial,          // ( 5 )
    const REALI *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:univ_maxcell1 )
    int  *_MaxAtom,
    int  *_MaxAtomCls,
    int  *_num_atom_cls,
    int  *_ncel_local,
    int  *_ncel_bound,
    int  *_ncel_max,
    int  *_cutoff_int,
    int  *_univ_maxcell,
    int  *_univ_maxcell1,
    int	 *_univ_ncell_nonzero,
    int	 *_univ_ncell_near,
    int	 *_univ_update,
    int  *_univ_mask2_size,
    int  *_univ_natom_max,
    int  *_check_virial,
    int  *_cpu_calc,
    REAL *_density,
    REAL *_cutoff2,
    REAL *_pairlistdist2,
    int  *_univ_gpu_start,
    REAL *_system_x,
    REAL *_system_y,
    REAL *_system_z
    )
{
    gpu_init();

    gpu_launch_compute_force_nonbond_table_ljpme_univ(
	_coord_pbc,
	_force,
        _ene_virial,
	_coord,
	_trans1,
	_cell_move,
	_charge,
	_atmcls,
	_natom,
	_nonb_lj12,
	_nonb_lj6,
	_nonb_lj6_factor,
	_table_grad,
	_univ_cell_pairlist1,
        _univ_mask2,
	_univ_ix_natom,
	_univ_ix_list,
        _univ_iy_natom,
	_univ_iy_list,
	_univ_ij_sort_list,
	_virial_check,
	*_MaxAtom,
	*_MaxAtomCls,
	*_num_atom_cls,
	*_ncel_local,
	*_ncel_bound,
	*_ncel_max,
	*_cutoff_int,
	*_univ_maxcell,
	*_univ_maxcell1,
	*_univ_ncell_nonzero,
	*_univ_ncell_near,
	*_univ_update,
        *_univ_mask2_size,
        *_univ_natom_max,
	*_check_virial,
        *_cpu_calc,
	*_density,
	*_cutoff2,
	*_pairlistdist2,
        *_univ_gpu_start,
        *_system_x, *_system_y, *_system_z
	);
}
/*
 *
 */

extern "C"
void gpu_launch_compute_force_nonbond_notable_univ_(
    REAL        *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        *_force,               // ( 1:MaxAtom, 1:3, 1:ncell_all )
    double      *_ene_virial,          // ( 5 )
    const REALI *_coord,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,              // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,           // ( 1:3, 1:ncel_max, 1:ncel_max )
    const REAL  *_charge,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_atmcls,              // ( 1:MaxAtom, 1:ncel_max )
    const int   *_natom,               // ( 1:ncel_max )
    const REAL  *_nonb_lj12,           // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6,            // ( 1:num_atom_cls, 1:num_atom_cls )
    const REAL  *_nonb_lj6_factor,     // ( 1:num_atom_cls )
    const REAL  *_table_grad,          // ( 1:6*cutoff_int )
    const int   *_univ_cell_pairlist1, // ( 1:2, 1:univ_maxcell )
    const char  *_univ_mask2,          // ( 1:univ_mask2_size, 1:univ_ncell_near)
    const int   *_univ_ix_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_ix_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_iy_natom,       // ( 1:univ_maxcell1 )
    const uchar *_univ_iy_list,        // ( 1:MaxAtom, 1:univ_maxcell1 )
    const int   *_univ_ij_sort_list,   // ( 1:univ_maxcell1 )
    const char  *_virial_check,        // ( 1:univ_maxcell1 )
    int  *_MaxAtom,
    int  *_MaxAtomCls,
    int  *_num_atom_cls,
    int  *_ncel_local,
    int  *_ncel_bound,
    int  *_ncel_max,
    int  *_cutoff_int,
    int  *_univ_maxcell,
    int  *_univ_maxcell1,
    int	 *_univ_ncell_nonzero,
    int	 *_univ_ncell_near,
    int	 *_univ_update,
    int  *_univ_mask2_size,
    int  *_univ_natom_max,
    int  *_check_virial,
    int  *_cpu_calc,
    REAL *_density,
    REAL *_cutoff2,
    REAL *_pairlistdist2,
    int  *_univ_gpu_start,
    REAL *_system_x,
    REAL *_system_y,
    REAL *_system_z
    )
{
    gpu_init();

    gpu_launch_compute_force_nonbond_notable_univ(
	_coord_pbc,
	_force,
        _ene_virial,
	_coord,
	_trans1,
	_cell_move,
	_charge,
	_atmcls,
	_natom,
	_nonb_lj12,
	_nonb_lj6,
	_nonb_lj6_factor,
	_table_grad,
	_univ_cell_pairlist1,
        _univ_mask2,
	_univ_ix_natom,
	_univ_ix_list,
        _univ_iy_natom,
	_univ_iy_list,
	_univ_ij_sort_list,
	_virial_check,
	*_MaxAtom,
	*_MaxAtomCls,
	*_num_atom_cls,
	*_ncel_local,
	*_ncel_bound,
	*_ncel_max,
	*_cutoff_int,
	*_univ_maxcell,
	*_univ_maxcell1,
	*_univ_ncell_nonzero,
	*_univ_ncell_near,
	*_univ_update,
        *_univ_mask2_size,
        *_univ_natom_max,
	*_check_virial,
        *_cpu_calc,
	*_density,
	*_cutoff2,
	*_pairlistdist2,
        *_univ_gpu_start,
        *_system_x, *_system_y, *_system_z
	);
}


void gpu_wait_compute_energy_nonbond_table_linear_univ(
    REAL        *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        *_force,               // ( 1:MaxAtom, 1:3, 1:ncell_all )
    double      *_ene_virial
    )
{
#if 1 /* for debug */
    cudaMemsetAsync( dummy_buf, 0, 1, p_stream[0] );
    cudaStreamSynchronize( p_stream[0] );
#endif

    for( int n = 0; n < NUM_STREAM; n++ ) {
        cudaStreamSynchronize( stream[n] );
    }
    gpu_memcpy_d2h_energy(
        _coord_pbc,
        _force,
        _ene_virial
        );
}



void gpu_wait_compute_force_nonbond_table_linear_univ(
    REAL        *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        *_force,               // ( 1:MaxAtom, 1:3, 1:ncell_all )
    double      *_ene_virial,          // ( 1:3, 1:MaxAtom, 1:ncell_all )
    int         check_virial
    )
{
#if 1 /* for debug */
    cudaMemsetAsync( dummy_buf, 0, 1, p_stream[0] );
    cudaStreamSynchronize( p_stream[0] );
#endif

    // cudaStreamSynchronize( stream[0] );

    gpu_memcpy_d2h_force(
	_coord_pbc,
	_force,
	_ene_virial,
        check_virial ,
	stream[0]
	);

    cudaStreamSynchronize( stream[0] );
}

/*
 * wapper function called from fortran subroutine
 */
extern "C"
void gpu_wait_compute_energy_nonbond_table_linear_univ_(
    REAL        *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        *_force,               // ( 1:MaxAtom, 1:3, 1:ncell_all )
    double      *_ene_virial
    )
{
    gpu_wait_compute_energy_nonbond_table_linear_univ(
        _coord_pbc,
        _force,
        _ene_virial );
}


/*
 * wapper function called from fortran subroutine
 */
extern "C"
void gpu_wait_compute_force_nonbond_table_linear_univ_(
    REAL        *_coord_pbc,           // ( 1:MaxAtom, 1:3, 1:ncel_max )
    REAL        *_force,               // ( 1:MaxAtom, 1:3, 1:ncell_all )
    double      *_ene_virial,
    int         *_check_virial
    )
{
    gpu_wait_compute_force_nonbond_table_linear_univ(
	_coord_pbc,
	_force,
        _ene_virial,
       *_check_virial );
}

/*
 *
 */
void gpu_launch_build_pairlist(
    REAL        *_coord_pbc,            // ( 1:MaxAtom, 1:3, 1:ncel_max ), AKA trans2
    const REALI *_coord,                // ( 1:MaxAtom, 1:3, 1:ncel_max )
    const REAL  *_trans1,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,            // ( 1:3, 1:ncel_max, 1:ncel_max )
    const int   *_natom,                // ( 1:ncel_max )
    const int   *_univ_cell_pairlist1,  // ( 1:2, 1:univ_maxcell )
    uchar       *_univ_ix_list,         // ( 1:MaxAtom, 1:univ_maxcell1? )
    uchar       *_univ_iy_list,         // ( 1:MaxAtom, 1:univ_maxcell1? )
    int         *_univ_ix_natom,        // ( 1:univ_maxcell1? )
    int         *_univ_iy_natom,        // ( 1:univ_maxcell1? )
    int         MaxAtom,
    int         ncel_local,
    int         ncel_bound,
    int         ncel_max,
    int         univ_maxcell,
    int         univ_maxcell1,
    REAL        pairdist2,
    REAL        cutoffdist2,
    REAL        system_x,
    REAL        system_y,
    REAL        system_z
)
{
#ifdef DEBUG
    printf( "[%s,%s,%d]\n", __FILE__, __func__, __LINE__ );
    printf( "_coord_pbc = %p\n", _coord_pbc );
    printf( "_coord = %p\n", _coord );
    printf( "_trans1 = %p\n", _trans1 );
    printf( "_cell_move = %p\n", _cell_move );
    printf( "MaxAtom = %d\n", MaxAtom );
    printf( "ncel_local = %d\n", ncel_local );
    printf( "ncel_bound = %d\n", ncel_bound );
    printf( "ncel_max   = %d\n", ncel_max );
    printf( "univ_maxcell = %d\n", univ_maxcell );
    printf( "univ_maxcell1 = %d\n", univ_maxcell1 );
    printf( "pairdist2 = %f\n", pairdist2 );
    printf( "cutoffdist2 = %f\n", cutoffdist2 );
#endif

    static int first = 1;

    if ( first ) {
	gpu_init_buffer_pairlist(
	    _coord_pbc,
	    _coord,
	    _trans1,
	    _cell_move,
	    _natom,
	    _univ_cell_pairlist1,
	    _univ_ix_list,
	    _univ_iy_list,
	    _univ_ix_natom,
	    _univ_iy_natom,
	    MaxAtom,
	    ncel_local,
	    ncel_bound,
	    ncel_max,
	    univ_maxcell,
	    univ_maxcell1 );
    }

    gpu_memcpy_h2d_pairlist(
	_coord_pbc,
	_coord,
	_trans1,
	_cell_move,
	_natom,
	_univ_cell_pairlist1,
	_univ_ix_list,
	_univ_iy_list,
	_univ_ix_natom,
	_univ_iy_natom,
	first,
	stream[0] );

    first = 0;

    dim3  def_dim3(1,1,1);
    dim3  num_block;
    dim3  num_thread;

    /* */
    num_block = def_dim3;
    num_block.x = ncel_local + ncel_bound;
    num_thread = def_dim3;
    num_thread.x = 256;

    kern_compute_force_nonbond_table_linear_univ__coord_pbc<<< num_block, num_thread, 0, stream[0] >>>(
        dev_coord_pbc,
        dev_coord,
        dev_trans1,
        dev_natom,
        MaxAtom,
        ncel_local,
        ncel_bound,
        ncel_max );
    /* */
#define MAX_ATOM 256
//#define MAX_ATOM 128
#define NUM_WARP_BUILD_PAIRLIST 4
    num_block = def_dim3;
    num_block.x = (univ_maxcell + NUM_WARP_BUILD_PAIRLIST - 1) / NUM_WARP_BUILD_PAIRLIST;
    num_thread = def_dim3;
    num_thread.x = 32; /* do not change */
    num_thread.y = NUM_WARP_BUILD_PAIRLIST;
    assert( MaxAtom < MAX_ATOM );

    kern_build_pairlist<NUM_WARP_BUILD_PAIRLIST, MAX_ATOM><<< num_block, num_thread, 0, stream[0] >>>(
	dev_coord_pbc,
	dev_cell_move,
	dev_natom,
	dev_univ_cell_pairlist1,
	dev_univ_ix_list,
	dev_univ_iy_list,
	dev_univ_ix_natom,
	dev_univ_iy_natom,
	MaxAtom,
	ncel_local,
	ncel_bound,
	ncel_max,
	univ_maxcell,
	univ_maxcell1,
	pairdist2,
	cutoffdist2,
        system_x, system_y, system_z );

    flag_build_pairlist_on_GPU = 1;
}

/* called from frotran */
extern "C"
void gpu_launch_build_pairlist_(
    REAL        *_coord_pbc,            // ( 1:MaxAtom, 1:3, 1:ncel_max ), AKA trans2
    const REALI *_coord,                // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const REAL  *_trans1,               // ( 1:3, 1:MaxAtom, 1:ncel_max )
    const char  *_cell_move,            // ( 1:3, 1:ncel_max, 1:ncel_max )
    const int   *_natom,                // ( 1:ncel_max )
    const int   *_univ_cell_pairlist1,  // ( 1:2, 1:univ_maxcell )
    uchar       *_univ_ix_list,         // ( 1:MaxAtom, 1:univ_maxcell1? )
    uchar       *_univ_iy_list,         // ( 1:MaxAtom, 1:univ_maxcell1? )
    int         *_univ_ix_natom,        // ( 1:univ_maxcell1? )
    int         *_univ_iy_natom,        // ( 1:univ_maxcell1? )
    int         *_MaxAtom,
    int         *_ncel_local,
    int         *_ncel_bound,
    int         *_ncel_max,
    int         *_univ_maxcell,
    int         *_univ_maxcell1,
    REAL        *_pairdist2,
    REAL        *_cutoffdist2,
    REAL        *_system_x,
    REAL        *_system_y,
    REAL        *_system_z,
    void        *_dummy
    )
{
    gpu_launch_build_pairlist(
	_coord_pbc,
	_coord,
	_trans1,
	_cell_move,
	_natom,
	_univ_cell_pairlist1,
	_univ_ix_list,
	_univ_iy_list,
	_univ_ix_natom,
	_univ_iy_natom,
	*_MaxAtom,
	*_ncel_local,
	*_ncel_bound,
	*_ncel_max,
	*_univ_maxcell,
	*_univ_maxcell1,
	*_pairdist2,
	*_cutoffdist2,
        *_system_x, *_system_y, *_system_z );
}

void gpu_wait_build_pairlist(
    uchar       *_univ_ix_list,
    uchar       *_univ_iy_list,
    int         *_univ_ix_natom,
    int         *_univ_iy_natom,
    uchar       *_univ_mask2,
    int         univ_mask2_size,
    int         univ_ncell_near
    )
{
#if 1 /* for debug */
    cudaMemsetAsync( dummy_buf, 0, 1, p_stream[0] );
    cudaStreamSynchronize( p_stream[0] );
#endif

    /* transfer univ_mask2 to GPU */
    /*
    size_univ_mask2 = sizeof(char) * univ_mask2_size * univ_ncell_near;
    if (max_size_univ_mask2 < size_univ_mask2){
	max_size_univ_mask2 = size_univ_mask2;
	if ( dev_univ_mask2 ) {
	    CUDA_CALL( cudaFree( dev_univ_mask2 ) );
	    dev_univ_mask2 = NULL;
	}
	CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_mask2, max_size_univ_mask2 ) );
    }
    */
//  CUDA_CALL( cudaMemcpyAsync( dev_univ_mask2, _univ_mask2, size_univ_mask2, cudaMemcpyHostToDevice, stream[1] ) );

//  cudaEventRecord( event[1], stream[1] ); /* dev_univ_mask2 */

    /* */
    // CUDA_CALL( cudaStreamSynchronize( stream[0] ) );

    gpu_memcpy_d2h_pairlist(
	_univ_ix_list,
	_univ_iy_list,
	_univ_ix_natom,
	_univ_iy_natom,
	stream[0],
	0 );

    CUDA_CALL( cudaStreamSynchronize( stream[0] ) );
}

/* called from frotran */
extern "C"
void gpu_wait_build_pairlist_(
    uchar       *_univ_ix_list,
    uchar       *_univ_iy_list,
    int         *_univ_ix_natom,
    int         *_univ_iy_natom,
    uchar       *_univ_mask2,
    int         *_univ_mask2_size,
    int         *_univ_ncell_near
    )
{
    gpu_init();

    gpu_wait_build_pairlist(
	_univ_ix_list,
	_univ_iy_list,
	_univ_ix_natom,
	_univ_iy_natom,
	_univ_mask2,
	*_univ_mask2_size,
	*_univ_ncell_near );
}

/* */
#define blockidx(X)   _blockidx[(X)-1]
#define blockidy(X)   _blockidy[(X)-1]
#define blockidz(X)   _blockidz[(X)-1]
#define blockid(X,Y)  _blockid[CALIDX2((X)-1,MaxAtom, (Y)-1,ncell)]

#define tasktable(X,Y)  _tasktable[CALIDX2((X),bnum, (Y)-1,ncell)]
#define tasktable_compact(X,Y)  _tasktable_compact[CALIDX2((X),bnum, (Y)-1,ncell)]
#define cellist(X,Y)    _cellist[CALIDX2((X),bnum, (Y),27)]

#define  COMPACT_TABLE_CNUM  16
// #define  COMPACT_TABLE_CNUM  32

/*
 *
 */
template <int BLENX, int BLENY, int BLENZ>
__global__ void kern_pme_recip_build_qdf_step1(
    int        *_vi,
    REAL       *_bsc,
    REAL       *_bscd,
    const REALI*_coord,
    const int  *_natom,
    const REAL *_r_scale,
    const int  *_ngrid,
    const int  MaxAtom,
    const int  ncell,
    const int  nlocalx,
    const int  nlocaly,
    const int  nlocalz,
    const int  *_blockidx,
    const int  *_blockidy,
    const int  *_blockidz,
    int        *_blockid,
    char       *_tasktable,
    char       *_tasktable_compact )
{
    int i = blockIdx.x + 1;
    if ( i > ncell ) return;

    for ( int ix_k = threadIdx.x; ix_k < (natom(i) * 3); ix_k += blockDim.x ) {
	int  ix = (ix_k / 3) + 1;
	int  k  = (ix_k % 3) + 1;

	REAL  vr = coord(k,ix,i) * r_scale(k);
	vr = vr + (REAL)(ngrid(k)/2) - (REAL)ngrid(k) * (REAL)NINT(vr/(REAL)ngrid(k));

	int  vi_tmp = (int)vr;
	if ( vi_tmp >= ngrid(k) ) {
	    vr = (REAL)ngrid(k) - 1e-5;
	    vi_tmp = ngrid(k)-1;
	}

	REAL  u = vr - (REAL)vi_tmp;
	REAL  u_2 = u * u;
	REAL  u_3 = u_2 * u;
	bsc(1,k,ix,i) =      u_3;
	bsc(2,k,ix,i) = -3.0*u_3 + 3.0*u_2 + 3.0*u + 1.0;
	bsc(3,k,ix,i) =  3.0*u_3 - 6.0*u_2         + 4.0;
	bsc(4,k,ix,i) = (1.0-u)*(1.0-u)*(1.0-u);
	bscd(1,k,ix,i) =      u_2;
	bscd(2,k,ix,i) = -3.0*u_2 + 2.0*u + 1.0;
	bscd(3,k,ix,i) =  3.0*u_2 - 4.0*u;
	bscd(4,k,ix,i) = -(1.0-u)*(1.0-u);

	vi(k,ix,i) = vi_tmp;
    }

    __syncthreads();

    int  bnumx = DIVCEIL(nlocalx, BLENX) + 1;
    int  bnumy = DIVCEIL(nlocaly, BLENY) + 1;
    int  bnumz = DIVCEIL(nlocalz, BLENZ) + 1;
    int  bnum = bnumx * bnumy * bnumz;

    int  i_compact = ((i-1)/COMPACT_TABLE_CNUM)+1;

    for ( int ix = 1 + threadIdx.x; ix <= natom(i); ix += blockDim.x ) {

	int gx = vi(1,ix,i) + 1;
	int gy = vi(2,ix,i) + 1;
	int gz = vi(3,ix,i) + 1;

	int bx = blockidx(gx);  // start from 0
	int by = blockidy(gy);  // start from 0
	int bz = blockidz(gz);  // start from 0

	int bid = -1;
	if ( bx >= 0 && by >= 0 && bz >= 0 ) {
	    bid = bx + bnumx * (by + bnumy * (bz));  // start from 0
	}

	if ( bid >= 0 ) {
	    blockid(ix,i) = bid;  // start from 0
	    tasktable(bid,i) = 1;
	    tasktable_compact(bid,i_compact) = 1;
	}
    }

}

/* */
template <int BLOCKDIMX, int SNUM>
__global__ void kern_pme_recip_build_qdf_step2(
    const int   bnum,
    const int   ncell,
    const char  *_tasktable,
    const char  *_tasktable_compact,
    short       *_cellist )
{
    int bid = threadIdx.x + blockDim.x * blockIdx.x;
    if ( bid >= bnum ) return;

    __shared__ short smem_cellist[SNUM][BLOCKDIMX];

    int  n = 0;
    for ( int i0 = 0; i0 < ncell; i0 += COMPACT_TABLE_CNUM ) {
	if ( tasktable_compact(bid, (i0/COMPACT_TABLE_CNUM)+1) == 0 ) continue;

	int i1_max = MIN(COMPACT_TABLE_CNUM, ncell-i0);
	for ( int i1 = 0; i1 < i1_max; i1++ ) {
	    int icel = i0+i1+1;
	    if ( tasktable(bid, icel) == 0 ) continue;

	    if ( n < SNUM ) {
	     	smem_cellist[n][threadIdx.x] = (short)icel;
	    }
	    else {
		cellist(bid, n) = (short)icel;
	    }
	    n++;
	}
    }

    int max_n = MIN(n, SNUM);
    for ( n = 0; n < max_n; n++ ) {
     	cellist(bid, n) = smem_cellist[n][threadIdx.x];
    }
}

/*
 *
 */
template <int BLENX, int BLENY, int BLENZ, int BLOCKDIMY>
__global__ void kern_pme_recip_build_qdf_step3(
    REAL       *_qdf_real,
    const int  *_vi,
    const REAL *_bsc,
    const REAL *_charge,
    const int  *_natom,
    const int  *_ngrid,
    int        MaxAtom,
    int        ncell,
    int        nlocalx,
    int        nlocaly,
    int        nlocalz,
    int        x_start,
    int        x_end,
    int        y_start,
    int        y_end,
    int        z_start,
    int        z_end,
    REAL       bs_fact3,
    const int  *_blockid,
    const char *_tasktable,
    const short *_cellist )
{
    const int  bnumx = DIVCEIL(nlocalx, BLENX) + 1;
    const int  bnumy = DIVCEIL(nlocaly, BLENY) + 1;
    const int  bnumz = DIVCEIL(nlocalz, BLENZ) + 1;
    const int  bnum = bnumx * bnumy * bnumz;

    int  bid = threadIdx.y + BLOCKDIMY * blockIdx.x;
    if ( bid >= bnum ) return;

    int  bidx =  bid % bnumx;
    int  bidy = (bid / bnumx) % bnumy;
    int  bidz = (bid / bnumx) / bnumy;

    if ( (bidz%2) == 1 ) bidy = bnumy-1 - bidy;
    if ( (bidy%2) == 1 ) bidx = bnumx-1 - bidx;
    bid = bidx + bnumx * (bidy + bnumy * bidz);

    if ( cellist(bid, 0) <= 0 ) return; /* */

    int  our_x_start = x_start + bidx * BLENX;
    int  our_y_start = y_start + bidy * BLENY;
    int  our_z_start = z_start + bidz * BLENZ;
    int  our_x_end = MIN(our_x_start+BLENX-1, x_end);
    int  our_y_end = MIN(our_y_start+BLENY-1, y_end);
    int  our_z_end = MIN(our_z_start+BLENZ-1, z_end);

    if ( bidx == bnumx-1 && nlocalx != ngrid(1) ) {
	our_x_start = ((x_end + 1) - 1) % ngrid(1) + 1;
	our_x_end   = ((x_end + 3) - 1) % ngrid(1) + 1;
    }
    if ( bidy == bnumy-1 && nlocaly != ngrid(2) ) {
	our_y_start = ((y_end + 1) - 1) % ngrid(2) + 1;
	our_y_end   = ((y_end + 3) - 1) % ngrid(2) + 1;
    }
    if ( bidz == bnumz-1 && nlocalz != ngrid(3) ) {
	our_z_start = ((z_end + 1) - 1) % ngrid(3) + 1;
	our_z_end   = ((z_end + 3) - 1) % ngrid(3) + 1;
    }

    if ( our_x_start > our_x_end ) return;
    if ( our_y_start > our_y_end ) return;
    if ( our_z_start > our_z_end ) return;

    __shared__ REAL smem_qdf[BLOCKDIMY][3+BLENZ][3+BLENY][3+BLENX];
    for ( int i = threadIdx.x; i < (3+BLENZ)*(3+BLENY)*(3+BLENX); i+= blockDim.x ) {
	int  sx =  i % (3+BLENX);
	int  sy = (i / (3+BLENX)) % (3+BLENY);
	int  sz = (i / (3+BLENX)) / (3+BLENY);
	smem_qdf[threadIdx.y][sz][sy][sx] = 0.0;
    }

    const int  my_ix = 1 +  threadIdx.x % 4;       // 1..4
    const int  my_iy = 1 + (threadIdx.x / 4) % 4;  // 1..4
    const int  my_iz = 1 + (threadIdx.x / 4) / 4;  // 1..2

    for ( int ic = 0; ; ic++ ) {

	int icel = (int) cellist(bid, ic);
	if ( icel <= 0 ) break;

	for ( int i0 = 1; i0 <= natom(icel); i0 += 32 ) {

	    int  ratom = natom(icel) - i0 + 1;
	    int  tmp_bid = -1;
	    if ( threadIdx.x < ratom ) {
		tmp_bid = blockid(i0+threadIdx.x, icel);
	    }
	    unsigned int vote = __ballot( (int)(tmp_bid == bid) );

	    int  i = i0;
	    while ( vote ) {
		int pos = __ffs(vote);
		i += (pos-1);
		vote = vote >> pos;

		int gx = vi(1,i,icel) + 1;  // 1..ngrid(1)
		int gy = vi(2,i,icel) + 1;  // 1..ngrid(2)
		int gz = vi(3,i,icel) + 1;  // 1..ngrid(3)
		gx = gx - (my_ix-1);
		gy = gy - (my_iy-1);
		gz = gz - (my_iz-1);

		int sx = gx - our_x_start + 3;
		int sy = gy - our_y_start + 3;
		int sz = gz - our_z_start + 3;

		REAL val = bs_fact3 * charge(i,icel) * bsc(my_ix,1,i,icel) * bsc(my_iy,2,i,icel);
		smem_qdf[threadIdx.y][sz  ][sy][sx] += val * bsc(my_iz  ,3,i,icel);
		smem_qdf[threadIdx.y][sz-2][sy][sx] += val * bsc(my_iz+2,3,i,icel);

		i++;
	    }
	}
    }

    for ( int i = threadIdx.x; i < (3+BLENZ)*(3+BLENY)*(3+BLENX); i+= blockDim.x ) {
	int sx =  i % (3+BLENX);
	int sy = (i / (3+BLENX)) % (3+BLENY);
	int sz = (i / (3+BLENX)) / (3+BLENY);

	int gx = sx + our_x_start - 3;
	int gy = sy + our_y_start - 3;
	int gz = sz + our_z_start - 3;
	if ( gx <= 0 ) gx += ngrid(1);
	if ( gy <= 0 ) gy += ngrid(2);
	if ( gz <= 0 ) gz += ngrid(3);

	if ( gx < x_start || gx > x_end ) continue;
	if ( gy < y_start || gy > y_end ) continue;
	if ( gz < z_start || gz > z_end ) continue;

	int ixs = gx - x_start + 1;
	int iys = gy - y_start + 1;
	int izs = gz - z_start + 1;
	// if ( ixs < 1 || ixs > nlocalx ) continue;
	// if ( iys < 1 || iys > nlocaly ) continue;
	// if ( izs < 1 || izs > nlocalz ) continue;

	REAL val = smem_qdf[threadIdx.y][sz][sy][sx];
	atomicAdd( &(qdf(ixs,iys,izs)), val );
    }
}

/*
 *
 */
void gpu_pme_recip_build_qdf(
    REAL        *_qdf_real, // ( 1:nlocalx*nlocaly*nlocalz )
    const REALI *_coord,    // ( 1:3, 1:MaxAtom, 1:ncell )
    const REAL  *_charge,   // ( 1:MaxAtom, 1:ncell )
    const int   *_natom,    // ( 1:ncell )
    int         MaxAtom,
    int         ncell,
    int         nlocalx,
    int         nlocaly,
    int         nlocalz,
    const int   *_ngrid,
    int         x_start,
    int         x_end,
    int         y_start,
    int         y_end,
    int         z_start,
    int         z_end,
    const REAL  *_r_scale,
    REAL        bs_fact3,
    const int   *_vi,
    const REAL  *_bsc,
    const REAL  *_bscd
    )
{
#ifdef DEBUG
// #if 1
    printf( "[%s,%s,%d]\n", __FILE__, __func__, __LINE__ );
    printf( "_qdf_real = %p\n", _qdf_real );
    printf( "_coord = %p\n", _coord );
    printf( "_charge = %p\n", _charge );
    printf( "_natom = %p\n", _natom );
    printf( "MaxAtom = %d\n", MaxAtom );
    printf( "ncell = %d\n", ncell );
    printf( "nlocalx/y/z = %d/%d/%d\n", nlocalx, nlocaly, nlocalz );
    printf( "_ngrid = %p\n", _ngrid );
    printf( "x_stat/end = %d/%d\n", x_start, x_end );
    printf( "y_stat/end = %d/%d\n", y_start, y_end );
    printf( "z_stat/end = %d/%d\n", z_start, z_end );
    printf( "_scale = %p\n", _r_scale );
    printf( "bs_fact3 = %f\n", bs_fact3 );
    printf( "_vi = %p\n", _vi );
    printf( "_bsc = %p\n", _bsc );
    printf( "_bscd = %p\n", _bscd );
#endif

    static int first = 1;

    static int  *dev_blockidx;
    static int  *dev_blockidy;
    static int  *dev_blockidz;
    static int  *dev_blockid;
    static char *dev_tasktable;
    static char *dev_tasktable_compact;
    static short *dev_cellist;

    static int  size_blockid;
    static int  size_tasktable;
    static int  size_tasktable_compact;
    static int  size_cellist;

#define BLENX 5
#define BLENY 5
#define BLENZ 3

    int  bnumx = DIVCEIL(nlocalx, BLENX) + 1;
    int  bnumy = DIVCEIL(nlocaly, BLENY) + 1;
    int  bnumz = DIVCEIL(nlocalz, BLENZ) + 1;
    int  bnum = bnumx * bnumy * bnumz;

    if ( first ) {
	size_qdf_real = sizeof(REAL) * nlocalx * nlocaly * nlocalz;
	size_vi = sizeof(int) * 3 * MaxAtom * ncell;
	size_bsc = sizeof(REAL) * 4 * 3 * MaxAtom * ncell;
	size_bscd = sizeof(REAL) * 4 * 3 * MaxAtom * ncell;
	size_r_scale = sizeof(REAL) * 3;
	size_ngrid = sizeof(int) * 3;

	tmp_force = (REAL*)malloc( size_force );
	tmp_qdf_real = (REAL*) malloc( size_qdf_real );
	tmp_vi = (int*) malloc( size_vi );
	tmp_bsc = (REAL*) malloc( size_bsc );
	tmp_bscd = (REAL*) malloc( size_bscd );

	CUDA_CALL( cudaMalloc_WN( (void**) &dev_qdf_real, size_qdf_real ) );
	CUDA_CALL( cudaMalloc_WN( (void**) &dev_vi, size_vi ) );
	CUDA_CALL( cudaMalloc_WN( (void**) &dev_bsc, size_bsc ) );
	CUDA_CALL( cudaMalloc_WN( (void**) &dev_bscd, size_bscd ) );
	CUDA_CALL( cudaMalloc_WN( (void**) &dev_r_scale, size_r_scale ) );
	CUDA_CALL( cudaMalloc_WN( (void**) &dev_ngrid, size_ngrid ) );

	CUDA_CALL( cudaMemcpy( dev_r_scale, _r_scale, size_r_scale, cudaMemcpyHostToDevice ) );
	CUDA_CALL( cudaMemcpy( dev_ngrid, _ngrid, size_ngrid, cudaMemcpyHostToDevice ) );

	size_t  size_blockidx = sizeof(int) * ngrid(1);
	size_t  size_blockidy = sizeof(int) * ngrid(2);
	size_t  size_blockidz = sizeof(int) * ngrid(3);
	int  *_blockidx = (int*) malloc( size_blockidx );
	int  *_blockidy = (int*) malloc( size_blockidy );
	int  *_blockidz = (int*) malloc( size_blockidz );
	CUDA_CALL( cudaMalloc( (void**) &dev_blockidx, size_blockidx ) );
	CUDA_CALL( cudaMalloc( (void**) &dev_blockidy, size_blockidy ) );
	CUDA_CALL( cudaMalloc( (void**) &dev_blockidz, size_blockidz ) );

	memset( _blockidx, -1, sizeof(int) * ngrid(1) );
	for ( int bidx = 0; bidx < bnumx; bidx++ ) {
	    int  our_x_start = x_start + bidx * BLENX;
	    int  our_x_end = MIN(our_x_start+BLENX-1, x_end);
	    if ( bidx == bnumx-1 && nlocalx != ngrid(1) ) {
		our_x_start = ((x_end + 1) - 1) % ngrid(1) + 1;
		our_x_end   = ((x_end + 3) - 1) % ngrid(1) + 1;
	    }
	    for ( int gx = our_x_start; gx <= our_x_end; gx++ ) {
		blockidx(gx) = bidx;
	    }
	}
	CUDA_CALL( cudaMemcpy( dev_blockidx, _blockidx, size_blockidx, cudaMemcpyHostToDevice ) );
#if 0
	printf( "blockidx:" );
	for ( int gx = 1; gx <= ngrid(1); gx++ ) {
	    printf( " %d", blockidx(gx) );
	}
	printf( "\n" );
	exit( 0 );
#endif

	memset( _blockidy, -1, sizeof(int) * ngrid(2) );
	for ( int bidy = 0; bidy < bnumy; bidy++ ) {
	    int  our_y_start = y_start + bidy * BLENY;
	    int  our_y_end = MIN(our_y_start+BLENY-1, y_end);
	    if ( bidy == bnumy-1 && nlocaly != ngrid(2) ) {
		our_y_start = ((y_end + 1) - 1) % ngrid(2) + 1;
		our_y_end   = ((y_end + 3) - 1) % ngrid(2) + 1;
	    }
	    for ( int gy = our_y_start; gy <= our_y_end; gy++ ) {
		blockidy(gy) = bidy;
	    }
	}
	CUDA_CALL( cudaMemcpy( dev_blockidy, _blockidy, size_blockidy, cudaMemcpyHostToDevice ) );
#if 0
	printf( "blockidy:" );
	for ( int gy = 1; gy <= ngrid(2); gy++ ) {
	    printf( " %d", blockidy(gy) );
	}
	printf( "\n" );
	exit( 0 );
#endif

	memset( _blockidz, -1, sizeof(int) * ngrid(3) );
	for ( int bidz = 0; bidz < bnumz; bidz++ ) {
	    int  our_z_start = z_start + bidz * BLENZ;
	    int  our_z_end = MIN(our_z_start+BLENZ-1, z_end);
	    if ( bidz == bnumz-1 && nlocalz != ngrid(3) ) {
		our_z_start = ((z_end + 1) - 1) % ngrid(3) + 1;
		our_z_end   = ((z_end + 3) - 1) % ngrid(3) + 1;
	    }
	    for ( int gz = our_z_start; gz <= our_z_end; gz++ ) {
		blockidz(gz) = bidz;
	    }
	}
	CUDA_CALL( cudaMemcpy( dev_blockidz, _blockidz, size_blockidz, cudaMemcpyHostToDevice ) );
#if 0
	printf( "blockidz:" );
	for ( int gz = 1; gz <= ngrid(3); gz++ ) {
	    printf( " %d", blockidz(gz) );
	}
	printf( "\n" );
	exit( 0 );
#endif

	free( _blockidx );
	free( _blockidy );
	free( _blockidz );

	size_blockid = sizeof(int) * MaxAtom * ncell;
	CUDA_CALL( cudaMalloc( (void**) &dev_blockid, size_blockid ) );

	size_tasktable = sizeof(char) * bnum * ncell;
	CUDA_CALL( cudaMalloc( (void**) &dev_tasktable, size_tasktable ) );

	size_tasktable_compact = sizeof(char) * bnum * DIVCEIL(ncell,COMPACT_TABLE_CNUM);
	CUDA_CALL( cudaMalloc( (void**) &dev_tasktable_compact, size_tasktable_compact ) );

	size_cellist = sizeof(short) * 27 * bnum;
	CUDA_CALL( cudaMalloc( (void**) &dev_cellist, size_cellist ) );
    }

    dim3  def_dim3(1,1,1);
    dim3  num_block;
    dim3  num_thread;

    /* */
    CUDA_CALL( cudaMemsetAsync( dev_blockid, -1, size_blockid, p_stream[0] ) );
    CUDA_CALL( cudaMemsetAsync( dev_tasktable, 0, size_tasktable, p_stream[0] ) );
    CUDA_CALL( cudaMemsetAsync( dev_tasktable_compact, 0, size_tasktable_compact, p_stream[0] ) );
    CUDA_CALL( cudaMemsetAsync( dev_cellist, 0, size_cellist, p_stream[0] ) );

    /* */
    num_block = def_dim3;
    num_block.x = ncell;
    num_thread = def_dim3;
    num_thread.x = 128;

    CUDA_CALL( cudaStreamWaitEvent( p_stream[0], event[0], 0 ) ); /* dev_coord, dev_natom */

    kern_pme_recip_build_qdf_step1<BLENX,BLENY,BLENZ><<< num_block, num_thread, 0, p_stream[0] >>>(
	dev_vi,
	dev_bsc,
	dev_bscd,
	dev_coord,
	dev_natom,
	dev_r_scale,
	dev_ngrid,
	MaxAtom,
	ncell,
	nlocalx, nlocaly, nlocalz,
	dev_blockidx, dev_blockidy, dev_blockidz, dev_blockid,
	dev_tasktable, dev_tasktable_compact );
    // CUDA_CALL( cudaStreamSynchronize( p_stream[0] ) );
    // CUDA_CALL( cudaMemcpy( tmp_vi, dev_vi, size_vi, cudaMemcpyDeviceToHost ) );
    // CUDA_CALL( cudaMemcpy( tmp_bsc, dev_bsc, size_bsc, cudaMemcpyDeviceToHost ) );
    // CUDA_CALL( cudaMemcpy( tmp_bscd, dev_bscd, size_bscd, cudaMemcpyDeviceToHost ) );

    /* */
    num_thread = def_dim3;
    num_thread.x = 128; /* do not change */
    num_block = def_dim3;
    num_block.x = DIVCEIL(bnum, num_thread.x);
    kern_pme_recip_build_qdf_step2<128, 8><<< num_block, num_thread, 0, p_stream[0] >>>(
	bnum, ncell, dev_tasktable, dev_tasktable_compact, dev_cellist );
    // CUDA_CALL( cudaStreamSynchronize( p_stream[0] ) );

#if 0
    int  *_tasktable = (int*) malloc( size_tasktable );
    CUDA_CALL( cudaMemcpy( _tasktable, dev_tasktable, size_tasktable, cudaMemcpyDeviceToHost ) );
    for ( int bid = 1; bid <= bnum; bid++ ) {
	printf( "bid:%d, icel:", bid );
	for ( int icel = 1; icel <= ncell; icel++ ) {
	    if ( tasktable(bid, icel) ) {
		printf( " %d", icel );
	    }
	}
	printf("\n");
    }
    free( _tasktable );
    exit( 0 );
#endif

    CUDA_CALL( cudaMemsetAsync( dev_qdf_real, 0, size_qdf_real, p_stream[0] ) );
    // CUDA_CALL( cudaStreamSynchronize( p_stream[0] ) );

    /* */
    num_thread = def_dim3;
    num_thread.x = 32; /* do not change */
    num_thread.y =  2;
    num_block = def_dim3;
    num_block.x = DIVCEIL(bnumx*bnumy*bnumz, num_thread.y);
    kern_pme_recip_build_qdf_step3<BLENX,BLENY,BLENZ, 2><<< num_block, num_thread, 0, p_stream[0] >>>(
	dev_qdf_real,
	dev_vi,
	dev_bsc,
	dev_charge,
	dev_natom,
	dev_ngrid,
	MaxAtom,
	ncell,
	nlocalx, nlocaly, nlocalz,
	x_start, x_end,
	y_start, y_end,
	z_start, z_end,
	bs_fact3,
	dev_blockid,
	dev_tasktable,
	dev_cellist );
    // CUDA_CALL( cudaStreamSynchronize( p_stream[0] ) );
    // CUDA_CALL( cudaMemcpy( tmp_qdf_real, dev_qdf_real, size_qdf_real, cudaMemcpyDeviceToHost ) );
    CUDA_CALL( cudaMemcpyAsync( _qdf_real, dev_qdf_real, size_qdf_real, cudaMemcpyDeviceToHost, p_stream[0] ) );

    CUDA_CALL( cudaStreamSynchronize( p_stream[0] ) );

    first = 0;

    /* */

    /*
     * CPU code for checking and reference
     */
    int MAX_NG = 200;
    int num_ng = 0;

#if 0
    for ( int i = 1; i <= ncell; i++ ) {
	for ( int ix = 1; ix <= natom(i); ix++ ) {

	    for ( int k = 1; k <= 3; k++ ) {
		REAL  vr = coord(k,ix,i) * r_scale(k);
		vr = vr + (REAL)(ngrid(k)/2) - (REAL)ngrid(k) * (REAL)NINT(vr/(REAL)ngrid(k));

		TMP_vi(k,ix,i) = (int)vr;
		if ( TMP_vi(k,ix,i) >= ngrid(k) ) {
		    vr = (REAL)ngrid(k) - 1e-5;
		    TMP_vi(k,ix,i) = ngrid(k)-1;
		}

		REAL  u = vr - (REAL)TMP_vi(k,ix,i);
		REAL  u_2 = u * u;
		REAL  u_3 = u_2 * u;
		TMP_bsc(1,k,ix,i) =      u_3;
		TMP_bsc(2,k,ix,i) = -3.0*u_3 + 3.0*u_2 + 3.0*u + 1.0;
		TMP_bsc(3,k,ix,i) =  3.0*u_3 - 6.0*u_2         + 4.0;
		TMP_bsc(4,k,ix,i) = (1.0-u)*(1.0-u)*(1.0-u);
		TMP_bscd(1,k,ix,i) =      u_2;
		TMP_bscd(2,k,ix,i) = -3.0*u_2 + 2.0*u + 1.0;
		TMP_bscd(3,k,ix,i) =  3.0*u_2 - 4.0*u;
		TMP_bscd(4,k,ix,i) = -(1.0-u)*(1.0-u);

	    }
	}
    }
#endif

#if 0
    /* check */
    for ( int i = 1; i <= ncell; i++ ) {
	for ( int ix = 1; ix <= natom(i); ix++ ) {
	    for ( int k = 1; k <= 3; k++ ) {

		if ( TMP_vi(k,ix,i) != vi(k,ix,i) ) {
		    printf( "mismatch: vi(%d,%d,%d) = %d, %d\n",
			     k, ix, i, vi(k,ix,i), TMP_vi(k,ix,i) );
		    ++num_ng;
		    if ( num_ng > MAX_NG ) exit( -1 );
		}

		for ( int n = 1; n <= 4; n++ ) {
		    if ( fabs(TMP_bsc(n,k,ix,i) - bsc(n,k,ix,i)) > 1e-5 ) {
			printf( "mismatch: bsc(%d,%d,%d,%d) = %f, %f\n",
				 n, k, ix, i, bsc(n,k,ix,i), TMP_bsc(n,k,ix,i) );
			++num_ng;
			if ( num_ng > MAX_NG ) exit( -1 );
		    }
		    if ( fabs(TMP_bscd(n,k,ix,i) - bscd(n,k,ix,i)) > 1e-5 ) {
			printf( "mismatch: bscd(%d,%d,%d,%d) = %f, %f\n",
				 n, k, ix, i, bscd(n,k,ix,i), TMP_bscd(n,k,ix,i) );
			++num_ng;
			if ( num_ng > MAX_NG ) exit( -1 );
		    }
		}
	    }
	}
    }
#endif

#if 0
    for ( int izs = 1; izs <= nlocalz; izs++ ) {
    for ( int iys = 1; iys <= nlocaly; iys++ ) {
    for ( int ixs = 1; ixs <= nlocalz; ixs++ ) {
	TMP_qdf(ixs,iys,izs) = 0.0;
    }}}

    for ( int icel = 1; icel <= ncell; icel++ ) {
	for ( int i = 1; i <= natom(icel); i++ ) {

	    for ( int iz = 1; iz <= 4; iz++ ) {
		int izs = vi(3,i,icel) - (iz-1) + 1;
		if (izs <= 0) izs = izs + ngrid(3);
		if (izs < z_start || izs > z_end) continue;
		izs = izs - z_start + 1;

	    for ( int iy = 1; iy <= 4; iy++ ) {
		int iys = vi(2,i,icel) - (iy-1) + 1;
		if (iys <= 0) iys = iys + ngrid(2);
		if (iys < y_start || iys > y_end) continue;
		iys = iys - y_start + 1;

	    for ( int ix = 1; ix <= 4; ix++ ) {
		int ixs = vi(1,i,icel) - (ix-1) + 1;
		if (ixs <= 0) ixs = ixs + ngrid(1);
		if (ixs < x_start || ixs > x_end) continue;
		ixs = ixs - x_start + 1;

		TMP_qdf(ixs,iys,izs) = TMP_qdf(ixs,iys,izs) + bs_fact3
		    * charge(i,icel) * bsc(ix,1,i,icel) * bsc(iy,2,i,icel) * bsc(iz,3,i,icel);

	    }}}
	}
    }
#endif

#if 0
    /* check */
    for ( int izs = 1; izs <= nlocalz; izs++ ) {
    for ( int iys = 1; iys <= nlocaly; iys++ ) {
    for ( int ixs = 1; ixs <= nlocalz; ixs++ ) {
	if ( fabs(TMP_qdf(ixs,iys,izs) - qdf(ixs,iys,izs)) > 1e-5 ) {
	    printf( "mismatch: qdf(%d,%d,%d) = %f, %f\n",
		    x_start+ixs-1, y_start+iys-1, z_start+izs-1,
		    TMP_qdf(ixs,iys,izs), qdf(ixs,iys,izs) );
	    ++num_ng;
	    if ( num_ng > MAX_NG ) exit( -1 );
	}
	qdf(ixs,iys,izs) = TMP_qdf(ixs,iys,izs);
    }}}
#endif

    if ( num_ng ) exit( -1 );
}

/* called from fortran */
extern "C"
void gpu_pme_recip_build_qdf_(
    REAL        *_qdf_real, // ( 1:nlocalx*nlocaly*nlocalz )
    const REALI *_coord,    // ( 1:3, 1:MaxAtom, 1:ncell )
    const REAL  *_charge,   // ( 1:MaxAtom, 1:ncell )
    const int   *_natom,    // ( 1:ncell )
    int         *_MaxAtom,
    int         *_ncell,
    int         *_nlocalx,
    int         *_nlocaly,
    int         *_nlocalz,
    const int   *_ngrid,    // ( 1:3 )
    int         *_x_start,
    int         *_x_end,
    int         *_y_start,
    int         *_y_end,
    int         *_z_start,
    int         *_z_end,
    const REAL  *_r_scale,  // ( 1:3 )
    REAL        *_bs_fact3,
    const int   *_vi,
    const REAL  *_bsc,
    const REAL  *_bscd,
    void        *_dummy
    )
{
    gpu_pme_recip_build_qdf(
	_qdf_real,
	_coord,
	_charge,
	_natom,
	*_MaxAtom,
	*_ncell,
	*_nlocalx, *_nlocaly, *_nlocalz,
	_ngrid,
	*_x_start, *_x_end,
	*_y_start, *_y_end,
	*_z_start, *_z_end,
	_r_scale,
	*_bs_fact3,
	_vi, _bsc, _bscd );
}

/*
 *
 */
__global__ void kern_pme_recip_interpolate_force_naive(
    REAL        *_force,    // ( 1:3, 1:MaxAtom, 1:ncell )
    const REAL  *_qdf_real, // ( 1:nlocalx*nlocaly*nlocalz )
    const REAL  *_charge,   // ( 1:MaxAtom, 1:ncell )
    const int   *_natom,    // ( 1:ncell )
    const REAL  *_r_scale,  // ( 1:3 )
    const int   *_ngrid,    // ( 1:3 )
    const int   *_vi,
    const REAL  *_bsc,
    const REAL  *_bscd,
    int         MaxAtom,
    int         ncell,
    int         nlocalx,
    int         nlocaly,
    int         nlocalz,
    int         x_start,
    int         x_end,
    int         y_start,
    int         y_end,
    int         z_start,
    int         z_end,
    REAL        bs_fact3d,
    REAL        vol_fact4 )
{
    const int icel = 1 + blockIdx.x;
    if ( icel > ncell ) return;

    for ( int i = 1 + threadIdx.x; i <= natom(icel); i += blockDim.x  ) {
	REAL f_1 = 0.0;
	REAL f_2 = 0.0;
	REAL f_3 = 0.0;

	for ( int iz = 1; iz <= 4; iz++ ) {
	    int izs = vi(3,i,icel) - (iz-1) + 1;
	    if (izs <= 0) izs = izs + ngrid(3);
	    if (izs < z_start || izs > z_end) continue;
	    izs = izs - z_start + 1;

	for ( int iy = 1; iy <= 4; iy++ ) {
	    int iys = vi(2,i,icel) - (iy-1) + 1;
	    if (iys <= 0) iys = iys + ngrid(2);
	    if (iys < y_start || iys > y_end) continue;
	    iys = iys - y_start + 1;

	for ( int ix = 1; ix <= 4; ix++ ) {
	    int ixs = vi(1,i,icel) - (ix-1) + 1;
	    if (ixs <= 0) ixs = ixs + ngrid(1);
	    if (ixs < x_start || ixs > x_end) continue;
	    ixs = ixs - x_start + 1;

	    f_1 = f_1 + qdf(ixs,iys,izs)
		* bscd(ix,1,i,icel) * bsc(iy,2,i,icel) * bsc(iz,3,i,icel);
	    f_2 = f_2 + qdf(ixs,iys,izs)
		* bsc(ix,1,i,icel) * bscd(iy,2,i,icel) * bsc(iz,3,i,icel);
	    f_3 = f_3 + qdf(ixs,iys,izs)
		* bsc(ix,1,i,icel) * bsc(iy,2,i,icel) * bscd(iz,3,i,icel);
	}}}

	force(1,i,icel) = force(1,i,icel)
	    - f_1 * charge(i,icel) * r_scale(1) * vol_fact4 * bs_fact3d;
	force(2,i,icel) = force(2,i,icel)
	    - f_2 * charge(i,icel) * r_scale(2) * vol_fact4 * bs_fact3d;
	force(3,i,icel) = force(3,i,icel)
	    - f_3 * charge(i,icel) * r_scale(3) * vol_fact4 * bs_fact3d;
    }
}

/*
 *
 */
template <int NATOM>
__global__ void kern_pme_recip_interpolate_force(
    REAL        *_force,    // ( 1:3, 1:MaxAtom, 1:ncell )
    const REAL  *_qdf_real, // ( 1:nlocalx*nlocaly*nlocalz )
    const REAL  *_charge,   // ( 1:MaxAtom, 1:ncell )
    const int   *_natom,    // ( 1:ncell )
    const REAL  *_r_scale,  // ( 1:3 )
    const int   *_ngrid,    // ( 1:3 )
    const int   *_vi,
    const REAL  *_bsc,
    const REAL  *_bscd,
    const int   MaxAtom,
    const int   ncell,
    const int   nlocalx,
    const int   nlocaly,
    const int   nlocalz,
    const int   x_start,
    const int   x_end,
    const int   y_start,
    const int   y_end,
    const int   z_start,
    const int   z_end,
    const REAL  bs_fact3d,
    const REAL  vol_fact4 )
{
    const int icel = 1 + blockIdx.x;
    if ( icel > ncell ) return;

    const int  ix = 1 +  threadIdx.x % 4;       // 1..4
    const int  iy = 1 + (threadIdx.x / 4) % 4;  // 1..4
    const int  iz = 1 + (threadIdx.x / 4) / 4;  // 1..2

    __shared__ REAL  smem_f[NATOM][3];

    for ( int i = 1 + threadIdx.y; i <= natom(icel); i += blockDim.y ) {

	REAL f_1 = 0.0;
	REAL f_2 = 0.0;
	REAL f_3 = 0.0;

	int gx = vi(1,i,icel) - (ix-1) + 1;  // 1..ngrid(1)
	int gy = vi(2,i,icel) - (iy-1) + 1;  // 1..ngrid(2)
	int gz = vi(3,i,icel) - (iz-1) + 1;  // 1..ngrid(3)

	if (gx <= 0) gx += ngrid(1);
	if (gx >= x_start && gx <= x_end) {
	    int ixs = gx - x_start + 1;

	    if (gy <= 0) gy += ngrid(2);
	    if (gy >= y_start && gy <= y_end) {
		int iys = gy - y_start + 1;

		REAL bsc_x_y  = bsc(ix,1,i,icel) * bsc(iy,2,i,icel);
		REAL bsc_dx_y = bscd(ix,1,i,icel) * bsc(iy,2,i,icel);
		REAL bsc_x_dy = bsc(ix,1,i,icel) * bscd(iy,2,i,icel);

		if (gz <= 0) gz += ngrid(3);
		if (gz >= z_start && gz <= z_end) {
		    int izs = gz - z_start + 1;

		    REAL qdf_tmp = __ldg(&qdf(ixs,iys,izs));
		    f_1 = f_1 + qdf_tmp	* bsc_dx_y * bsc(iz,3,i,icel);
		    f_2 = f_2 + qdf_tmp	* bsc_x_dy * bsc(iz,3,i,icel);
		    f_3 = f_3 + qdf_tmp	* bsc_x_y * bscd(iz,3,i,icel);
		}

		gz = gz - 2;
		if (gz <= 0) gz += ngrid(3);
		if (gz >= z_start && gz <= z_end) {
		    int izs = gz - z_start + 1;

		    REAL qdf_tmp = __ldg(&qdf(ixs,iys,izs));
		    f_1 = f_1 + qdf_tmp	* bsc_dx_y * bsc(iz+2,3,i,icel);
		    f_2 = f_2 + qdf_tmp	* bsc_x_dy * bsc(iz+2,3,i,icel);
		    f_3 = f_3 + qdf_tmp	* bsc_x_y * bscd(iz+2,3,i,icel);
		}
	    }
	}

	WARP_RSUM_12345( f_1 );
	if ( threadIdx.x == 0 ) smem_f[i-1][1-1] = f_1;
	WARP_RSUM_12345( f_2 );
	if ( threadIdx.x == 0 ) smem_f[i-1][2-1] = f_2;
	WARP_RSUM_12345( f_3 );
	if ( threadIdx.x == 0 ) smem_f[i-1][3-1] = f_3;
    }

    __syncthreads();

    {
	const int  tid = threadIdx.x + threadIdx.y * blockDim.x;
	const int  tnum = blockDim.x * blockDim.y;

	for ( int i_k = tid; i_k < natom(icel) * 3; i_k += tnum ) {
	    int i = (i_k / 3) + 1;
	    int k = (i_k % 3) + 1;
	    force(k,i,icel) = force(k,i,icel)
		- smem_f[i-1][k-1] * charge(i,icel) * r_scale(k) * vol_fact4 * bs_fact3d;
	}
    }
}

/*
 *
 */
void gpu_pme_recip_interpolate_force(
    REAL        *_force,    // ( 1:3, 1:MaxAtom, 1:ncell )
    const REAL  *_qdf_real, // ( 1:nlocalx*nlocaly*nlocalz )
    const REAL  *_charge,   // ( 1:MaxAtom, 1:ncell )
    const int   *_natom,    // ( 1:ncell )
    int         MaxAtom,
    int         ncell,
    int         nlocalx,
    int         nlocaly,
    int         nlocalz,
    const int   *_ngrid,    // ( 1:3 )
    int         x_start,
    int         x_end,
    int         y_start,
    int         y_end,
    int         z_start,
    int         z_end,
    const REAL  *_r_scale,  // ( 1:3 )
    REAL        bs_fact3d,
    REAL        vol_fact4,
    const int   *_vi,
    const REAL  *_bsc,
    const REAL  *_bscd
    )
{
#ifdef DEBUG
    printf( "[%s,%s,%d]\n", __FILE__, __func__, __LINE__ );
    printf( "_force = %p\n", _force );
    printf( "_qdf_real = %p\n", _qdf_real );
    printf( "_charge = %p\n", _charge );
    printf( "_natom = %p\n", _natom );
    printf( "MaxAtom = %d\n", MaxAtom );
    printf( "ncell = %d\n", ncell );
    printf( "nlocalx/y/z = %d/%d/%d\n", nlocalx, nlocaly, nlocalz );
    printf( "_ngrid = %p\n", _ngrid );
    printf( "ngridx/y/z = %d/%d/%d\n", _ngrid[0], _ngrid[1], _ngrid[2] );
    printf( "x_stat/end = %d/%d\n", x_start, x_end );
    printf( "y_stat/end = %d/%d\n", y_start, y_end );
    printf( "z_stat/end = %d/%d\n", z_start, z_end );
    printf( "_r_scale = %p\n", _r_scale );
    printf( "r_scalex/y/z = %f/%f/%f\n", _r_scale[0], _r_scale[1], _r_scale[2] );
    printf( "bs_fact3d = %f\n", bs_fact3d );
    printf( "vol_fact4 = %f\n", vol_fact4 );
#endif

#undef  DO_VERIFY
// #define DO_VERIFY

#ifdef DO_VERIFY
    /* for debug */
    for ( int icel = 1; icel <= ncell; icel++ ) {
	for ( int i = 1; i <= natom(icel); i++ ) {
	    for ( int k = 1; k <= 3; k++ ) {
		TMP_force(k,i,icel) = force(k,i,icel);
	    }
	}
    }
#endif

#undef  MAX_ATOM
#define MAX_ATOM (128+64)

    assert( MAX_ATOM > MaxAtom );

    dim3  def_dim3(1,1,1);
    dim3  num_block;
    dim3  num_thread;

    /* */
    CUDA_CALL( cudaMemcpyAsync( dev_qdf_real, _qdf_real, size_qdf_real, cudaMemcpyHostToDevice, p_stream[0] ) );
    CUDA_CALL( cudaMemsetAsync( dev_force_recip, 0, size_force, p_stream[0] ) );

    /* */
    num_thread = def_dim3;
    num_thread.x = 32; /* do not change */
    num_thread.y = 4;
    num_block = def_dim3;
    num_block.x = ncell;
    kern_pme_recip_interpolate_force<MAX_ATOM><<< num_block, num_thread, 0, p_stream[0] >>>(
	dev_force_recip,
	dev_qdf_real, dev_charge, dev_natom, dev_r_scale, dev_ngrid, dev_vi, dev_bsc, dev_bscd,
	MaxAtom, ncell, nlocalx, nlocaly, nlocalz,
	x_start, x_end, y_start, y_end, z_start, z_end,
	bs_fact3d, vol_fact4 );

    /* */
    CUDA_CALL( cudaMemcpyAsync( _force, dev_force_recip, size_force, cudaMemcpyDeviceToHost, p_stream[0] ) );

    CUDA_CALL( cudaStreamSynchronize( p_stream[0] ) );

#ifdef DO_VERIFY
    /*
     * reference on CPU
     */
    int max_NG = 200;
    int num_NG = 0;
    for ( int icel = 1; icel <= ncell; icel++ ) {
	for ( int i = 1; i <= natom(icel); i++ ) {
	    REAL f_1 = 0.0;
	    REAL f_2 = 0.0;
	    REAL f_3 = 0.0;

	    for ( int iz = 1; iz <= 4; iz++ ) {
		int izs = vi(3,i,icel) - (iz-1) + 1;
		if (izs <= 0) izs = izs + ngrid(3);
		if (izs < z_start || izs > z_end) continue;
		izs = izs - z_start + 1;

	    for ( int iy = 1; iy <= 4; iy++ ) {
		int iys = vi(2,i,icel) - (iy-1) + 1;
		if (iys <= 0) iys = iys + ngrid(2);
		if (iys < y_start || iys > y_end) continue;
		iys = iys - y_start + 1;

	    for ( int ix = 1; ix <= 4; ix++ ) {
		int ixs = vi(1,i,icel) - (ix-1) + 1;
		if (ixs <= 0) ixs = ixs + ngrid(1);
		if (ixs < x_start || ixs > x_end) continue;
		ixs = ixs - x_start + 1;

		f_1 = f_1 + qdf(ixs,iys,izs)
		    * bscd(ix,1,i,icel) * bsc(iy,2,i,icel) * bsc(iz,3,i,icel);
		f_2 = f_2 + qdf(ixs,iys,izs)
		    * bsc(ix,1,i,icel) * bscd(iy,2,i,icel) * bsc(iz,3,i,icel);
		f_3 = f_3 + qdf(ixs,iys,izs)
		    * bsc(ix,1,i,icel) * bsc(iy,2,i,icel) * bscd(iz,3,i,icel);
	    }}}

	    TMP_force(1,i,icel) = TMP_force(1,i,icel)
		- f_1 * charge(i,icel) * r_scale(1) * vol_fact4 * bs_fact3d;
	    TMP_force(2,i,icel) = TMP_force(2,i,icel)
		- f_2 * charge(i,icel) * r_scale(2) * vol_fact4 * bs_fact3d;
	    TMP_force(3,i,icel) = TMP_force(3,i,icel)
		- f_3 * charge(i,icel) * r_scale(3) * vol_fact4 * bs_fact3d;
	}
    }

    /* check force() */
    for ( int icel = 1; icel <= ncell; icel++ ) {
	for ( int i = 1; i <= natom(icel); i++ ) {
	    for ( int k = 1; k <= 3; k++ ) {
		if ( fabs(force(k,i,icel) - TMP_force(k,i,icel)) > 1e-4 ) {
		    if ( fabs((force(k,i,icel) - TMP_force(k,i,icel)) / force(k,i,icel)) > 1e-4 ) {
			printf( "mismatch: force(%d,%d,%d) = %f, %f\n",
				k, i, icel,
				force(k,i,icel), TMP_force(k,i,icel) );
			++num_NG;
			if ( num_NG > max_NG ) exit( -1 );
		    }
		}
	    }
	}
    }
    if ( num_NG > 0 ) exit( -1 );
#endif
}

/* */
extern "C"
void gpu_pme_recip_interpolate_force_(
    REAL        *_force,    // ( 1:3, 1:MaxAtom, 1:ncell )
    const REAL  *_qdf_real, // ( 1:nlocalx*nlocaly*nlocalz )
    const REAL  *_charge,   // ( 1:MaxAtom, 1:ncell )
    const int   *_natom,    // ( 1:ncell )
    int         *_MaxAtom,
    int         *_ncell,
    int         *_nlocalx,
    int         *_nlocaly,
    int         *_nlocalz,
    const int   *_ngrid,    // ( 1:3 )
    int         *_x_start,
    int         *_x_end,
    int         *_y_start,
    int         *_y_end,
    int         *_z_start,
    int         *_z_end,
    const REAL  *_r_scale,  // ( 1:3 )
    REAL        *_bs_fact3d,
    REAL        *_vol_fact4,
    const int   *_vi,
    const REAL  *_bsc,
    const REAL  *_bscd,
    void        *_dummy
    )
{
    gpu_pme_recip_interpolate_force(
	_force,    // ( 1:3, 1:MaxAtom, 1:ncell )
	_qdf_real, // ( 1:nlocalx*nlocaly*nlocalz )
	_charge,   // ( 1:MaxAtom, 1:ncell )
	_natom,    // ( 1:ncell )
	*_MaxAtom,
	*_ncell,
	*_nlocalx,
	*_nlocaly,
	*_nlocalz,
	_ngrid,    // ( 1:3 )
	*_x_start,
	*_x_end,
	*_y_start,
	*_y_end,
	*_z_start,
	*_z_end,
	_r_scale,  // ( 1:3 )
	*_bs_fact3d,
	*_vol_fact4,
	_vi,
	_bsc,
	_bscd
	);
}

extern "C"
void gpu_pme_recip_keep_force_(
    REAL        *_force,    // ( 1:3, 1:MaxAtom, 1:ncell )
    const int   *_natom,    // ( 1:ncell )
    int         *_MaxAtom,
    int         *_ncell )
{
    const int  MaxAtom = *_MaxAtom;
    const int  ncell = *_ncell;

    for ( int icel = 1; icel <= ncell; icel++ ) {
	for ( int i = 1; i <= natom(icel); i++ ) {
	    for ( int k = 1; k <= 3; k++ ) {
		TMP_force(k,i,icel) = force(k,i,icel);
	    }
	}
    }
}

/*
 * JJ : cudaMalloc and cudaGet..., and cudaStream.. are performed
        at the first execution of gpu_init
 */
void gpu_init()
{
    static int first = 1;
    if ( first ) {
	first = 0;
	cudaMalloc( &dummy_buf, 1 );

	cudaGetDeviceProperties( &prop, 0 );

	int n;
	for ( n = 0; n < NUM_STREAM; n++ ) {
	    cudaStreamCreate( & (stream[n]) );
	}
	for ( n = 0; n < NUM_EVENT; n++ ) {
	    cudaEventCreate( & (event[n]) );
	}

	int  p_low, p_high;
	cudaDeviceGetStreamPriorityRange ( &p_low, &p_high );
	for ( n = 0; n < NUM_P_STREAM; n++ ) {
	    cudaStreamCreateWithPriority( &(p_stream[n]), cudaStreamDefault, p_high );
	}
    }
}

/*
 * JJ : gpu_init_ is called from frotran
 */
extern "C"
void gpu_init_( )
{
    gpu_init();
}

void gpu_upload_lj_coeffs( const int  &num_atom_cls,
                           const REAL *nonb_lj12,
                           const REAL *nonb_lj6 ) {
  if ( dev_nonb_lj12 != NULL ) {
    // upload only if the array is allocated
    CUDA_CALL( cudaMemcpy( dev_nonb_lj12, nonb_lj12, size_nonb_lj12,
                           cudaMemcpyHostToDevice ) );
  }

  if ( dev_nonb_lj6 != NULL ) {
    // upload only if the array is allocated
    CUDA_CALL( cudaMemcpy( dev_nonb_lj6, nonb_lj6, size_nonb_lj6,
                           cudaMemcpyHostToDevice ) );
  }
}

// MK : upload LJ coeffs if necessary
extern "C"
void gpu_upload_lj_coeffs_( const int  &num_atom_cls,
                            const REAL *nonb_lj12,
                            const REAL *nonb_lj6 ) {
  gpu_init();

  gpu_upload_lj_coeffs( num_atom_cls, nonb_lj12, nonb_lj6 );
}

void gpu_upload_charge( const REAL *charge ) {
  if ( dev_charge != NULL ) {
    // upload only if the array is allocated
    CUDA_CALL( cudaMemcpy( dev_charge, charge, size_charge,
                           cudaMemcpyHostToDevice ) );
  }

}

// JJ : upload charge if necessary
extern "C"
void gpu_upload_charge_( const REAL *charge) {
  gpu_init();

  gpu_upload_charge( charge );
}

extern "C"
int set_pinned_memory_(
        void *ptr,
        int *insize
        )
{
    size_t size = *insize;
    cudaHostRegister( ptr, size, cudaHostRegisterPortable );
    return 0;
}

extern "C"
int unset_pinned_memory_(
        void *ptr
        )
{
    cudaHostUnregister( ptr );
    return 0;
}

#define BLOCK_SIZE 256
__global__ void gpu_unpack_data(unsigned int *do_data, char *di_data, size_t size) {
    static const unsigned int table[16] =
          {0x00000000, 0x00000001, 0x00000100, 0x00000101, 0x00010000, 0x00010001, 0x00010100, 0x00010101,
           0x01000000, 0x01000001, 0x01000100, 0x01000101, 0x01010000, 0x01010001, 0x01010100, 0x01010101};

    size_t idx = blockIdx.x*BLOCK_SIZE + threadIdx.x;

    unsigned char bufchar, loadchar;
    unsigned int outint=0;

    __shared__ unsigned int outbuff[BLOCK_SIZE*2];

    if( idx < size ){
        loadchar = di_data[idx];

        bufchar = loadchar&0x0f;
        outint = table[bufchar];
        outbuff[threadIdx.x*2] = outint;

        bufchar = loadchar>>4;
        outint = table[bufchar];
        outbuff[threadIdx.x*2+1] = outint;
    }
    __syncthreads();

    if( (blockIdx.x+1)*BLOCK_SIZE <= size ){
        do_data[blockIdx.x*BLOCK_SIZE*2 + threadIdx.x] = outbuff[threadIdx.x];
        do_data[blockIdx.x*BLOCK_SIZE*2 + BLOCK_SIZE + threadIdx.x] = outbuff[BLOCK_SIZE + threadIdx.x];
    }else{
        int runthread = (size - blockIdx.x*BLOCK_SIZE)*2;
        if( threadIdx.x < runthread ){
            do_data[blockIdx.x*BLOCK_SIZE*2 + threadIdx.x] = outbuff[threadIdx.x];
            runthread -= BLOCK_SIZE;
        }
        if( runthread > 0 ){
            if ( threadIdx.x < runthread ){
                do_data[blockIdx.x*BLOCK_SIZE*2 + BLOCK_SIZE + threadIdx.x] = outbuff[BLOCK_SIZE + threadIdx.x];
            }
        }
    }
}

void gpu_copy_mask2(
    uchar       *_pack_univ_mask2,
    int         univ_mask2_size,
    int         univ_ncell_near,
    int         start_pack,
    int         end_pack,
    int         in_stream_no
    )
{
    int stream_no=1;

    if( in_stream_no < 3 ){
        stream_no = in_stream_no;
    }else{
        stream_no = in_stream_no%3;
    }

    size_t copysize = (end_pack - start_pack + 1)*sizeof(char);
    CUDA_CALL( cudaMemcpyAsync( dev_pack_univ_mask2+start_pack,
                                _pack_univ_mask2+start_pack,
                                copysize, cudaMemcpyHostToDevice, stream[in_stream_no] ) );

    const int blockCount = (copysize+BLOCK_SIZE-1)/BLOCK_SIZE;
    gpu_unpack_data<<<blockCount, BLOCK_SIZE, 0, stream[in_stream_no]>>>
                                 ((unsigned int *)(dev_univ_mask2+start_pack*8),
                                  dev_pack_univ_mask2+start_pack,
                                  (size_t)(end_pack - start_pack + 1));

    cudaEventRecord( event[in_stream_no], stream[in_stream_no] ); /* dev_univ_mask2 */
}

extern "C"
void gpu_copy_mask2_(
    uchar       *_pack_univ_mask2,
    int         *univ_mask2_size,
    int         *univ_ncell_near,
    int         *start_pack,
    int         *end_pack,
    int         *in_stream_no
    )
{
    gpu_copy_mask2(
        _pack_univ_mask2,
        *univ_mask2_size,
        *univ_ncell_near,
        *start_pack,
        *end_pack,
        *in_stream_no );
}


void gpu_allocate_packdata(
    uchar       *_pack_univ_mask2,
    int         univ_mask2_size,
    int         univ_ncell_near
    )
{
    size_t size_pack_univ_mask2;
    size_pack_univ_mask2 = sizeof(char)*(univ_mask2_size * univ_ncell_near+7)/8;
    size_univ_mask2      = sizeof(char)* univ_mask2_size * univ_ncell_near;

    if (max_size_univ_mask2 < size_univ_mask2){
        max_size_univ_mask2 = size_univ_mask2;
        if ( dev_univ_mask2 ) {
            CUDA_CALL( cudaFree( dev_univ_mask2 ) );
            dev_univ_mask2 = NULL;
        }
        size_t malloc_size = size_pack_univ_mask2 * 8;
        CUDA_CALL( cudaMalloc_WN( (void**) &dev_univ_mask2, malloc_size ) );

        if ( dev_pack_univ_mask2 ) {
            CUDA_CALL( cudaFree( dev_pack_univ_mask2 ) );
            dev_pack_univ_mask2 = NULL;
        }
        CUDA_CALL( cudaMalloc_WN( (void**) &dev_pack_univ_mask2, size_pack_univ_mask2 ) );
    }
}

extern "C"
void gpu_allocate_packdata_(
    uchar       *_pack_univ_mask2,
    int         *univ_mask2_size,
    int         *univ_ncell_near
    )
{
    gpu_allocate_packdata(
        _pack_univ_mask2,
        *univ_mask2_size,
        *univ_ncell_near);
}

/*
 * JJ : NVTX Range Push
 */
//extern "C"
//void nvtx_range_push_id_(
//    int         *_id
//    )
//{
//    char rangid[128];
//    sprintf(rangid,"id=%d",*_id);
//    nvtxRangePushA(rangid);
//}

/*
 * JJ : NVTX Range Pop
 */
//extern "C"
//void nvtx_range_pop_()
//{
//    nvtxRangePop();
//}



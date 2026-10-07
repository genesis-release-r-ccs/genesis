#!/usr/bin/env python3
"""
Calculate unnormalized Tiwary/WT reweighting factors for ABMD.

Designed for:
  - Normal ABMD: bias-history .bp chunks + GENESIS standard-output log chunks
  - MW-ABMD: shared bias-history .bp chunks + one walker's remlog chunks
  - BE-ABMD: one bias-state .bp history + the matching/sorted remlog chunks

Inputs are supplied through two meta files:

  --bias-meta-file
      Text file listing time-contiguous .bp files, one path per line.

  --log-meta-file
      Text file listing time-contiguous GENESIS log or remlog files,
      one path per line. Each listed file is inspected by its contents.
      Both normal GENESIS output logs and REMD/remlog-style files are
      accepted as long as they contain INFO: records with STEP, TIME,
      and ABMD_BIAS.

Example meta files:

  meta_bias.out
      run_abmd1.bp
      run_abmd2.bp
      run_abmd3.bp

  meta_log.out
      run_abmd1.out
      run_abmd2.out
      run_abmd3.out

Typical usage:
  python calc_reweight.py \
      --bias-meta-file meta_bias.out \
      --log-meta-file meta_log.out \
      -t 300 -wt 3000 \
      -o weight.dat

For MW/BE-ABMD, meta_log.out may instead list remlog files. The same parser
is used for both normal logs and remlogs.

IMPORTANT:
  This script DOES NOT normalize the weights to sum to 1.
  It writes log_weight and relative_weight for this walker/replica.
  Global normalization across multiple walkers/replicas should be done later,
  after all unnormalized weights have been collected.

Because exp(log_weight) may overflow/underflow, relative_weight is shifted by
the maximum log_weight within THIS output file:
    relative_weight = exp(log_weight - max(log_weight))
This preserves all relative weights but is NOT globally normalized.
The output also records log_shift=max(log_weight), so different files can be
combined later in log space without loss of scale.
"""

import argparse
from pathlib import Path
import numpy as np


def logsumexp(x):
    x = np.asarray(x, dtype=float)
    m = np.max(x)
    if not np.isfinite(m):
        return m
    return m + np.log(np.sum(np.exp(x - m)))


def read_path_list(path):
    p = Path(path)
    base = p.parent
    paths = []

    with p.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            q = Path(s)
            if not q.is_absolute():
                q = base / q
            paths.append(q)

    if not paths:
        raise ValueError(f"No paths found in meta file: {path}")

    missing = [str(q) for q in paths if not q.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing files:\n  " + "\n  ".join(missing)
        )

    return paths


def validate_meta_paths(path, expected_kind):
    """Read a meta file and validate the listed input files."""
    files = read_path_list(path)

    if expected_kind == "bp":
        bad = [str(fp) for fp in files if fp.suffix.lower() != ".bp"]
        if bad:
            print(
                "WARNING: the following files listed in --bias-meta-file do not "
                "have a .bp suffix:\n  " + "\n  ".join(bad)
            )

    return files


def parse_bp_file(path, fix=0.0):
    records = []
    cur_step = None
    cur_total = None
    cur_grid = []

    nrc = None
    nbins = []

    def flush():
        nonlocal cur_step, cur_total, cur_grid
        if cur_step is not None:
            records.append(
                {
                    "step": cur_step,
                    "total_bias": cur_total,
                    "grid": np.asarray(cur_grid, dtype=float),
                }
            )
        cur_step = None
        cur_total = None
        cur_grid = []

    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if s.startswith("#"):
                parts = s.replace(",", " ").split()

                if len(parts) >= 3 and parts[1] == "nreac_coord":
                    nrc = int(parts[2])

                if len(parts) >= 10 and parts[1] == "i" and parts[2] == "kind":
                    nbins.append(int(parts[9]))

                if len(parts) >= 5 and parts[1] == "STEP" and parts[2] == "U":
                    if cur_step is not None:
                        flush()
                    cur_step = int(float(parts[-2]))
                    cur_total = float(parts[-1]) - fix

                continue

            if cur_step is not None:
                cur_grid.append(float(s.split()[-1]) - fix)

    if cur_step is not None:
        flush()

    if not records:
        raise ValueError(f"No bias snapshots found in {path}")

    if nbins:
        expected = int(np.prod(nbins))
        bad = [i + 1 for i, r in enumerate(records) if len(r["grid"]) != expected]
        if bad:
            raise ValueError(
                f"{path}: grid-size mismatch. Header says {nbins} -> {expected} bins; "
                f"mismatch in blocks {bad[:5]}"
            )

    return records, nrc, nbins


def load_bp_series(path, fix=0.0):
    files = validate_meta_paths(path, "bp")

    all_records = []
    ref_nrc = None
    ref_bins = None

    for fp in files:
        records, nrc, bins = parse_bp_file(fp, fix=fix)

        if ref_nrc is None:
            ref_nrc = nrc
            ref_bins = bins
        else:
            if nrc is not None and ref_nrc is not None and nrc != ref_nrc:
                raise ValueError(f"nreac_coord changed between BP chunks: {fp}")
            if bins and ref_bins and bins != ref_bins:
                raise ValueError(f"Grid dimensions changed between BP chunks: {fp}")

        all_records.extend(records)

    return all_records, files, ref_nrc, ref_bins


def parse_genesis_log_file(path):
    records = []
    header_names = None

    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if s.startswith("INFO:") and "ABMD_BIAS" in s:
                header_names = s.split()[1:]
                continue

            if not s.startswith("INFO:") or header_names is None:
                continue

            vals = s.split()[1:]
            if len(vals) < len(header_names):
                continue

            try:
                idx_step = header_names.index("STEP")
                idx_time = header_names.index("TIME")
                idx_bias = header_names.index("ABMD_BIAS")

                step = int(float(vals[idx_step]))

                # Skip the initial INFO record at STEP 0.
                # Some remlog chunks contain this extra record while the BP
                # history starts from the first nonzero output step.
                if step == 0:
                    continue

                records.append(
                    {
                        "step": step,
                        "time": float(vals[idx_time]),
                        "abmd_bias": float(vals[idx_bias]),
                    }
                )
            except (ValueError, IndexError):
                continue

    if not records:
        raise ValueError(f"No INFO records containing STEP/TIME/ABMD_BIAS found in {path}")

    return records


def load_log_series(path):
    # --log-meta-file is always a meta file. Each listed file may be either
    # a normal GENESIS output log or a remlog; parse by contents, not filename.
    files = validate_meta_paths(path, "log")
    records = []
    for fp in files:
        records.extend(parse_genesis_log_file(fp))
    return records, files


def compute_log_weight(grid, applied_bias, beta, alpha):
    """
    Tiwary/WT weight used by the previous scripts:

      w(t) =
        sum_j exp[ beta * alpha       * (V_j(t)-V(s_t,t)) ]
        ---------------------------------------------------
        sum_j exp[ beta * (1+alpha)   * (V_j(t)-V(s_t,t)) ]

    with alpha = T / DeltaT.

    Returned in log form.
    """
    dv = grid - applied_bias
    log_num = logsumexp(beta * alpha * dv)
    log_den = logsumexp(beta * (1.0 + alpha) * dv)
    return log_num - log_den


def main():
    ap = argparse.ArgumentParser(
        description="Calculate unnormalized WT/Tiwary reweighting factors for one walker/replica."
    )
    ap.add_argument(
        "--bias-meta-file",
        required=True,
        help="meta file listing time-contiguous .bp files, one path per line",
    )
    ap.add_argument(
        "--log-meta-file",
        required=True,
        help=(
            "meta file listing time-contiguous GENESIS log or remlog files; "
            "each listed file is parsed by contents for INFO: STEP/TIME/ABMD_BIAS"
        ),
    )
    ap.add_argument("-t", "--temperature", type=float, default=300.0)
    ap.add_argument(
        "-wt",
        "--welltemp",
        type=float,
        default=10000000.0,
        help="GENESIS well-tempered temperature DeltaT [K]",
    )
    ap.add_argument("--fix", "-f", type=float, default=0.0)
    ap.add_argument(
        "--skip",
        type=int,
        default=0,
        help="skip this many aligned snapshots from the beginning",
    )
    ap.add_argument(
        "--ndata",
        type=int,
        default=None,
        help="use at most this many aligned snapshots after --skip",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help=(
            "write detailed columns (snapshot, BP/log STEP, time, ABMD_BIAS, "
            "log_weight, weight). Without --verbose, write only snapshot and log_weight."
        ),
    )
    ap.add_argument(
        "-o",
        "--output",
        default="weight_unnorm.dat",
        help="output weight file",
    )
    args = ap.parse_args()

    beta = (1000.0 * 4.184) / (8.314 * args.temperature)
    alpha = args.temperature / args.welltemp

    bp_records, bp_files, nrc, bins = load_bp_series(
        args.bias_meta_file, fix=args.fix
    )
    log_records, log_files = load_log_series(args.log_meta_file)

    start = args.skip
    stop = None if args.ndata is None else start + args.ndata

    bp_records = bp_records[start:stop]
    log_records = log_records[start:stop]

    if len(bp_records) != len(log_records):
        raise ValueError(
            f"Snapshot-count mismatch after slicing: "
            f"{len(bp_records)} BP snapshots vs {len(log_records)} log snapshots.\n"
            "Make sure the BP/log output frequencies and the ordering in both meta files match."
        )

    if not bp_records:
        raise ValueError("No snapshots remain after slicing.")

    step_mismatch = sum(
        1 for b, r in zip(bp_records, log_records) if b["step"] != r["step"]
    )
    if step_mismatch:
        print(
            f"WARNING: STEP differs in {step_mismatch}/{len(bp_records)} aligned snapshots. "
            "Alignment is by concatenated snapshot order."
        )

    log_weights = np.asarray(
        [
            compute_log_weight(
                b["grid"],
                r["abmd_bias"] - args.fix,
                beta,
                alpha,
            )
            for b, r in zip(bp_records, log_records)
        ],
        dtype=float,
    )

    # Do NOT normalize to sum=1 here.
    # Use a local shift only to make an exportable positive relative weight.
    log_shift = float(np.max(log_weights))
    rel_weights = np.exp(log_weights - log_shift)

    out = Path(args.output)
    if out.parent != Path("."):
        out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w") as f:
        f.write(f"# temperature_K {args.temperature:.12g}\n")
        f.write(f"# welltemp_DeltaT_K {args.welltemp:.12g}\n")
        f.write(f"# alpha_T_over_DeltaT {alpha:.12g}\n")
        f.write(f"# log_shift {log_shift:.17g}\n")
        f.write("# NOTE: weights are NOT normalized to sum=1.\n")
        f.write("# weight = exp(log_weight - log_shift)\n")

        if args.verbose:
            f.write(
                "# snapshot step_bp step_log time_ps ABMD_BIAS "
                "log_weight weight\n"
            )
            for i, (b, r, lw, rw) in enumerate(
                zip(bp_records, log_records, log_weights, rel_weights), 1
            ):
                f.write(
                    f"{i:d} {b['step']:d} {r['step']:d} {r['time']:.8f} "
                    f"{(r['abmd_bias'] - args.fix):.12g} "
                    f"{lw:.17g} {rw:.17g}\n"
                )
        else:
            f.write("# snapshot log_weight\n")
            for i, lw in enumerate(log_weights, 1):
                f.write(f"{i:d} {lw:.17g}\n")

    print(f"Wrote {len(log_weights)} unnormalized weights to {out}")
    print(f"log_shift = {log_shift:.12g}")
    print(f"sum(relative_weight) = {np.sum(rel_weights):.12g}  (not normalized)")
    print(f"BP files: {len(bp_files)}")
    print(f"log/remlog files: {len(log_files)}")
    print(f"nreac_coord: {nrc}")
    print(f"grid bins: {bins}")
    print(f"alpha=T/DeltaT: {alpha:.12g}")


if __name__ == "__main__":
    main()

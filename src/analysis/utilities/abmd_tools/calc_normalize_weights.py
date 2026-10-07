#!/usr/bin/env python3
"""
Globally normalize reweighting factors from multiple walkers/replicas.

Input:
  meta.txt
    weight_w1.dat
    weight_w2.dat
    ...
    weight_w10.dat

Each weight file is assumed to be produced by calc_reweight.py and
must contain a column named "log_weight" in its header.

The script:
  1) reads all walkers listed in meta.txt,
  2) collects every log_weight,
  3) computes ONE global normalization constant over all walkers/snapshots,
  4) writes one globally normalized output file per walker.

For every output walker k:
    global_weight(k,t) = exp(log_weight(k,t) - log_Z_global)

Therefore:
    sum_over_all_walkers_and_snapshots(global_weight) = 1

The sum within one walker is generally NOT 1; it is that walker's total
contribution to the globally normalized ensemble.



Usage examples
--------------
Use all records from every walker (default):
    python calc_normalize_weights.py meta_global.dat -o global_weight

Use an explicit record range:
    python calc_normalize_weights.py meta_global.dat \
        --start 1 --end 100000 -o global_weight

Use cumulative blocks with the default block size (5000 records/block):
    python calc_normalize_weights.py meta_global.dat \
        --block 20 -o global_weight
    -> uses records 1-100000
    -> outputs weight_1-20_w1.dat, weight_1-20_w2.dat, ...

Use 2000 records as one block:
    python calc_normalize_weights.py meta_global.dat \
        --block 20 --block-size 2000 -o global_weight
    -> uses records 1-40000
    -> outputs weight_1-20_w1.dat, weight_1-20_w2.dat, ...

If --block is omitted, all records are used unless --start/--end are specified.
--block cannot be combined with --start/--end.
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


def read_meta(path):
    p = Path(path)
    base = p.parent
    files = []

    with p.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            q = Path(s)
            if not q.is_absolute():
                q = base / q
            files.append(q)

    if not files:
        raise ValueError(f"No weight files found in meta file: {path}")

    missing = [str(q) for q in files if not q.exists()]
    if missing:
        raise FileNotFoundError(
            "The following weight files do not exist:\n  "
            + "\n  ".join(missing)
        )

    return files


def read_weight_file(path):
    """
    Read a calc_reweight_onewalker.py output file.

    Returns
    -------
    comments : list[str]
        Original comment lines.
    columns : list[str]
        Column names parsed from the last '#' header line that contains
        'log_weight'.
    rows : list[list[str]]
        Data rows as strings.
    logw : np.ndarray
        log_weight values.
    """
    comments = []
    columns = None
    rows = []

    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if s.startswith("#"):
                comments.append(s)
                names = s[1:].strip().split()
                if "log_weight" in names:
                    columns = names
                continue

            rows.append(s.split())

    if columns is None:
        raise ValueError(
            f"{path}: could not find a header containing a 'log_weight' column."
        )

    idx = columns.index("log_weight")

    bad = [i + 1 for i, r in enumerate(rows) if len(r) <= idx]
    if bad:
        raise ValueError(
            f"{path}: malformed data rows at {bad[:5]}"
        )

    logw = np.asarray([float(r[idx]) for r in rows], dtype=float)

    if not np.all(np.isfinite(logw)):
        raise ValueError(f"{path}: non-finite log_weight values found.")

    return comments, columns, rows, logw


def output_path_for(input_path, outdir, suffix):
    p = Path(input_path)
    stem = p.stem
    ext = p.suffix if p.suffix else ".dat"
    return outdir / f"{stem}{suffix}{ext}"


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Globally normalize log weights from multiple walker weight files "
            "listed in a meta file, and write one output file per walker."
        )
    )
    ap.add_argument(
        "meta",
        help="meta file listing one weight file per walker/replica",
    )
    ap.add_argument(
        "-o",
        "--outdir",
        default="global_weights",
        help="output directory [default: global_weights]",
    )
    ap.add_argument(
        "--suffix",
        default="_global",
        help="suffix added to each output filename [default: _global]",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help=(
            "write all original columns plus global_weight. "
            "Without --verbose, write only snapshot and global_weight."
        ),
    )
    ap.add_argument(
        "--start",
        type=int,
        default=1,
        help="first data record to use from each weight file (1-based, inclusive) [default: 1]",
    )
    ap.add_argument(
        "--end",
        type=int,
        default=None,
        help="last data record to use from each weight file (1-based, inclusive) [default: last record]",
    )
    ap.add_argument(
        "--block",
        type=int,
        default=None,
        help=(
            "number of cumulative blocks to use from the beginning. "
            "The number of records per block is set by --block-size. "
            "Cannot be combined with --start or --end."
        ),
    )
    ap.add_argument(
        "--block-size",
        type=int,
        default=5000,
        help="number of data records per block [default: 5000]",
    )
    args = ap.parse_args()

    if args.block_size < 1:
        raise ValueError("--block-size must be >= 1")

    if args.block is not None:
        if args.block < 1:
            raise ValueError("--block must be >= 1")
        if args.start != 1 or args.end is not None:
            raise ValueError("--block cannot be combined with --start or --end")
        args.start = 1
        args.end = args.block * args.block_size

    if args.start < 1:
        raise ValueError("--start must be >= 1")
    if args.end is not None and args.end < args.start:
        raise ValueError("--end must be >= --start")

    if args.block is not None:
        print(
            f"Block mode: 1-{args.block}, "
            f"{args.block_size} records/block, "
            f"{args.end} records/walker"
        )

    files = read_meta(args.meta)

    datasets = []
    all_logw = []

    for fp in files:
        comments, columns, rows, logw = read_weight_file(fp)

        # Apply the same 1-based inclusive record range to every walker file.
        i0 = args.start - 1
        i1 = None if args.end is None else args.end

        rows = rows[i0:i1]
        logw = logw[i0:i1]

        if len(logw) == 0:
            raise ValueError(
                f"{fp}: no data records remain after applying "
                f"--start {args.start}"
                + ("" if args.end is None else f" --end {args.end}")
            )

        datasets.append(
            {
                "path": fp,
                "comments": comments,
                "columns": columns,
                "rows": rows,
                "logw": logw,
            }
        )
        all_logw.append(logw)

    all_logw = np.concatenate(all_logw)
    log_Z_global = logsumexp(all_logw)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    total = 0.0
    contributions = []

    for iw, ds in enumerate(datasets, 1):
        gw = np.exp(ds["logw"] - log_Z_global)
        contrib = float(np.sum(gw))
        total += contrib
        contributions.append((iw, ds["path"], contrib))

        if args.block is not None:
            # For block analysis, use a simple walker-based filename independent
            # of the local-weight input filename.
            outpath = outdir / f"weight_1-{args.block}_w{iw}_gobal.dat"
        else:
            outpath = output_path_for(ds["path"], outdir, args.suffix)

        with outpath.open("w") as f:
            # Pure numerical output only: no header/comment lines.
            if args.verbose:
                for row, w in zip(ds["rows"], gw):
                    f.write(" ".join(row) + f" {w:.17g}\n")
            else:
                # 1st column: snapshot index
                # 2nd column: globally normalized weight
                for row, w in zip(ds["rows"], gw):
                    f.write(f"{row[0]} {w:.17g}\n")

        print(
            f"[{iw:02d}] {ds['path']} -> {outpath} "
            f"(sum global_weight = {contrib:.12g})"
        )

    summary = outdir / "global_weight_summary.dat"
    with summary.open("w") as f:
        for iw, fp, contrib in contributions:
            f.write(f"{iw:d} {fp} {contrib:.17g}\n")

    print(f"\nTotal global weight = {total:.17g}")
    print(f"Summary written to {summary}")

    # Numerical sanity check.
    if not np.isclose(total, 1.0, rtol=1e-12, atol=1e-12):
        raise RuntimeError(
            f"Global normalization check failed: total weight = {total}"
        )


if __name__ == "__main__":
    main()

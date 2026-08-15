#!/usr/bin/env python3
"""
Collect per-target predictor runtime and memory into one committable CSV.

Each benchmark run writes predictor_runtime_stats.csv from Nextflow's trace,
holding wall-clock and memory for every predictor task.  Those trees are
gitignored and HPC-side, so this pulls the numbers into a single small file the
analysis can read without HPC access.

Rows are mapped to the same (model, msa) convention as combine_metrics.py so
they join against combined_metrics.csv.  Only COMPLETED tasks are kept, since a
failed task's elapsed time measures how long it took to die.

Usage:
    build_runtime_csv.py --tier 1 --output runtime_stats.csv
"""

import argparse
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from combine_metrics import map_model, to_float  # noqa: E402

KEEP = ["status", "queue", "elapsed_s", "standalone_elapsed_s", "pct_cpu",
        "rss_gb", "peak_rss_gb", "peak_vmem_gb"]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmarks-dir", default=os.path.normpath(
        os.path.join(here, "..", "experiments", "benchmarks")))
    ap.add_argument("--manifest", default=os.path.normpath(
        os.path.join(here, "..", "data", "benchmark_complexes.tsv")))
    ap.add_argument("--tier", type=int, default=None)
    ap.add_argument("--output", default="runtime_stats.csv")
    args = ap.parse_args()

    wanted = None
    if args.tier is not None:
        with open(args.manifest, newline="") as fh:
            wanted = {r["pdb"] for r in csv.DictReader(fh, delimiter="\t")
                      if int(r["tier"]) == args.tier}

    # Two sources.  The constrained variants were re-run separately after their
    # restraints were fixed, so their timings live under rerun_constrained/ and
    # must supersede the originals; every other model comes from the main tree.
    rows, n_files, n_skipped = [], 0, 0
    patterns = [
        os.path.join(args.benchmarks_dir, "*", "*_benchmark_results",
                     "predictor_runtime_stats.csv"),
        os.path.join(args.benchmarks_dir, "*", "rerun_constrained",
                     "*_constrained_rerun_results", "predictor_runtime_stats.csv"),
    ]
    seen = set()          # (pdb, predictor) already taken from the re-run
    paths = sorted(glob.glob(patterns[1])) + sorted(glob.glob(patterns[0]))
    for path in paths:
        is_rerun = "rerun_constrained" in path
        pdb = path.split(os.sep)[-4] if is_rerun else path.split(os.sep)[-3]
        if wanted is not None and pdb not in wanted:
            continue
        n_files += 1
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh):
                predictor = (r.get("model") or "").strip()
                if not predictor:
                    continue
                if (r.get("status") or "").strip() != "COMPLETED":
                    n_skipped += 1
                    continue
                if to_float(r.get("elapsed_s")) is None:
                    n_skipped += 1
                    continue
                key = (pdb, predictor)
                if is_rerun:
                    seen.add(key)
                elif key in seen:
                    continue          # superseded by the re-run
                model, msa = map_model(predictor)
                out = {"model": model, "msa": msa, "pdb": pdb,
                       "predictor": predictor}
                out.update({c: r.get(c, "") for c in KEEP})
                rows.append(out)

    if not rows:
        sys.exit("ERROR: no runtime rows found.")

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["model", "msa", "pdb", "predictor"] + KEEP)
        w.writeheader()
        w.writerows(rows)

    n_pdb = len({r["pdb"] for r in rows})
    print(f"Read {n_files} runtime files, skipped {n_skipped} non-COMPLETED rows")
    print(f"Wrote {len(rows)} rows across {n_pdb} targets -> {args.output}")


if __name__ == "__main__":
    main()

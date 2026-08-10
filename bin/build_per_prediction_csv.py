#!/usr/bin/env python3
"""
Emit EVERY prediction from the benchmark runs, not just the best one per combo.

`combine_metrics.py` keeps the single lowest-ra_eff prediction for each
(model, msa, pdb).  That is the right unit for "how good is this method", but it
hides the spread across the 5 seeds x 5 diffusion samples the pipeline runs, and
it selects for good poses.  Analyses that ask "how much does running one
prediction cost me versus running 25" or "how well does confidence separate good
from bad poses" need the unselected distribution.

This script walks the same benchmark result trees and writes one row per
prediction, with the same derived (model, msa, pdb) columns so the output joins
cleanly against combined_metrics.csv.

Only a small set of columns is carried through, to keep the output committable
(the full 30-column dump is several MB).  AF3's aggregate rows are dropped: they
are a summary over samples, not an independent prediction, so including them
would bias the per-prediction spread.

Usage:
    build_per_prediction_csv.py [--benchmarks-dir DIR] [--output FILE]
                                [--tier N --manifest FILE]
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from combine_metrics import MODEL_MAP, map_model, to_float  # noqa: E402

RA_EFF_COL = "rmsd_effector_receptor_aligned"

# Columns kept in the output, beyond the derived model/msa/pdb/predictor.
KEEP = [
    "model_name",
    "avg_plddt", "ptm", "iptm", "pae_mean", "ipae",
    "ipsae_ab", "ipsae_ba", "ipsae_min", "actifptm",
    "rmsd_receptor", "rmsd_effector_independent",
    RA_EFF_COL,
]
DERIVED = ["model", "msa", "pdb", "predictor"]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    default_bench = os.path.normpath(
        os.path.join(here, "..", "experiments", "benchmarks"))
    default_manifest = os.path.normpath(
        os.path.join(here, "..", "data", "benchmark_complexes.tsv"))

    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmarks-dir", default=default_bench)
    ap.add_argument("--output", default="per_prediction_metrics.csv")
    ap.add_argument("--manifest", default=default_manifest,
                    help="benchmark_complexes.tsv, used only with --tier")
    ap.add_argument("--tier", type=int, default=None,
                    help="restrict to one tier (default: all)")
    args = ap.parse_args()

    if not os.path.isdir(args.benchmarks_dir):
        sys.exit(f"ERROR: benchmarks dir not found: {args.benchmarks_dir}")

    wanted = None
    if args.tier is not None:
        if not os.path.isfile(args.manifest):
            sys.exit(f"ERROR: manifest not found: {args.manifest}")
        wanted = set()
        with open(args.manifest, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if int(row["tier"]) == args.tier:
                    wanted.add(row["pdb"])
        print(f"Tier {args.tier}: {len(wanted)} complexes in manifest")

    out_rows = []
    n_files = n_seen = n_agg = n_noscore = 0
    unknown = set()
    missing = []

    for pdb in sorted(os.listdir(args.benchmarks_dir)):
        if wanted is not None and pdb not in wanted:
            continue
        sub = os.path.join(args.benchmarks_dir, pdb)
        if not os.path.isdir(sub):
            continue
        csv_path = os.path.join(sub, f"{pdb}_benchmark_results", "all_metrics.csv")
        if not os.path.isfile(csv_path):
            missing.append(pdb)
            continue
        n_files += 1

        with open(csv_path, newline="") as fh:
            for row in csv.DictReader(fh):
                n_seen += 1
                raw_model = (row.get("model") or "").strip()
                if not raw_model:
                    continue
                if raw_model not in MODEL_MAP:
                    unknown.add(raw_model)

                # AF3 emits an aggregate-over-samples row; not an independent
                # prediction, so it must not enter a per-prediction spread.
                if "aggregate" in str(row.get("model_name", "")).lower():
                    n_agg += 1
                    continue

                if to_float(row.get(RA_EFF_COL)) is None:
                    n_noscore += 1
                    continue

                model, msa = map_model(raw_model)
                out = {"model": model, "msa": msa, "pdb": pdb,
                       "predictor": raw_model}
                for c in KEEP:
                    out[c] = row.get(c, "")
                out_rows.append(out)

    if not out_rows:
        sys.exit("ERROR: no predictions found.")

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=DERIVED + KEEP)
        w.writeheader()
        w.writerows(out_rows)

    n_pdb = len({r["pdb"] for r in out_rows})
    print(f"Read {n_files} result trees, {n_seen} rows")
    print(f"  dropped {n_agg} aggregate rows, {n_noscore} rows with no ra_eff")
    if unknown:
        print(f"  WARNING: unmapped predictor names: {sorted(unknown)}")
    if missing:
        print(f"  no all_metrics.csv for: {sorted(missing)}")
    print(f"Wrote {len(out_rows)} predictions across {n_pdb} targets "
          f"-> {args.output}")


if __name__ == "__main__":
    main()

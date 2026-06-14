#!/usr/bin/env python3
"""
Combine per-benchmark all_metrics.csv files into one summary CSV.

For every (model, msa, pdb) combination it keeps the single prediction with the
LOWEST ra_eff (receptor-aligned effector RMSD) — i.e. the best effector pose.

Output columns:
    model, msa, pdb, plddt, ra_eff, rec_rmsd, eff_rmsd

Source columns (from <PDB>_benchmark_results/all_metrics.csv):
    plddt    <- avg_plddt
    ra_eff   <- rmsd_effector_receptor_aligned
    rec_rmsd <- rmsd_receptor
    eff_rmsd <- rmsd_effector_independent

The "model" column in all_metrics.csv is the predictor directory name; we split
it into a base model + an msa flag ("msa" / "no_msa"):

    af2m -> af2m/msa            af3 -> af3/msa          af3_nomsa -> af3/no_msa
    boltz1 -> boltz1/no_msa     boltz1_msa -> boltz1/msa
    boltz2 -> boltz2/no_msa     boltz2_msa -> boltz2/msa
    boltz{1,2}_constrained -> kept as their own model, no_msa
    chai1 -> chai1/no_msa       esmfold2 -> esmfold2/no_msa
    colabfold -> colabfold/msa  colabfold_nomsa -> colabfold/no_msa

Usage:
    combine_metrics.py [--benchmarks-dir DIR] [--output FILE]
"""

import argparse
import csv
import os
import sys

# predictor dir name -> (model, msa_flag)
MODEL_MAP = {
    "af2m":               ("af2m",               "msa"),
    "af3":                ("af3",                "msa"),
    "af3_nomsa":          ("af3",                "no_msa"),
    "boltz1":             ("boltz1",             "no_msa"),
    "boltz1_msa":         ("boltz1",             "msa"),
    "boltz1_constrained": ("boltz1_constrained", "no_msa"),
    "boltz2":             ("boltz2",             "no_msa"),
    "boltz2_msa":         ("boltz2",             "msa"),
    "boltz2_constrained": ("boltz2_constrained", "no_msa"),
    "chai1":              ("chai1",              "no_msa"),
    "colabfold":          ("colabfold",          "msa"),
    "colabfold_nomsa":    ("colabfold",          "no_msa"),
    "esmfold2":           ("esmfold2",           "no_msa"),
}

SRC = {  # output field -> source column
    "plddt":    "avg_plddt",
    "ra_eff":   "rmsd_effector_receptor_aligned",
    "rec_rmsd": "rmsd_receptor",
    "eff_rmsd": "rmsd_effector_independent",
}
OUT_COLS = ["model", "msa", "pdb", "plddt", "ra_eff", "rec_rmsd", "eff_rmsd"]


def to_float(v):
    """Parse a float, or None if empty / non-numeric (NA, nan, '')."""
    if v is None:
        return None
    s = str(v).strip()
    if s == "" or s.lower() in ("na", "nan", "none", "null"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def map_model(name):
    """Return (model, msa). Fall back to a suffix heuristic for unknown names."""
    if name in MODEL_MAP:
        return MODEL_MAP[name]
    if name.endswith("_nomsa"):
        return (name[: -len("_nomsa")], "no_msa")
    if name.endswith("_msa"):
        return (name[: -len("_msa")], "msa")
    return (name, "no_msa")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    default_bench = os.path.normpath(os.path.join(here, "..", "experiments", "benchmarks"))

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmarks-dir", default=default_bench,
                    help="dir containing <PDB>/<PDB>_benchmark_results/ (default: %(default)s)")
    ap.add_argument("--output", default="combined_metrics.csv",
                    help="output CSV path (default: %(default)s)")
    args = ap.parse_args()

    bench_dir = args.benchmarks_dir
    if not os.path.isdir(bench_dir):
        sys.exit(f"ERROR: benchmarks dir not found: {bench_dir}")

    # best[(model, msa, pdb)] = (ra_eff_value, row_dict)
    best = {}
    n_files = 0
    n_rows = 0
    unknown_models = set()
    no_score = set()           # combos seen but with no valid ra_eff anywhere
    missing = []               # benchmark dirs without an all_metrics.csv

    for pdb in sorted(os.listdir(bench_dir)):
        sub = os.path.join(bench_dir, pdb)
        if not os.path.isdir(sub):
            continue
        csv_path = os.path.join(sub, f"{pdb}_benchmark_results", "all_metrics.csv")
        if not os.path.isfile(csv_path):
            missing.append(pdb)
            continue
        n_files += 1

        with open(csv_path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                n_rows += 1
                raw_model = (row.get("model") or "").strip()
                if not raw_model:
                    continue
                if raw_model not in MODEL_MAP:
                    unknown_models.add(raw_model)
                model, msa = map_model(raw_model)
                key = (model, msa, pdb)

                ra = to_float(row.get(SRC["ra_eff"]))
                if ra is None:
                    no_score.add(key)
                    continue
                if key not in best or ra < best[key][0]:
                    best[key] = (ra, {
                        "model":    model,
                        "msa":      msa,
                        "pdb":      pdb,
                        "plddt":    (row.get(SRC["plddt"])    or "").strip(),
                        "ra_eff":   (row.get(SRC["ra_eff"])   or "").strip(),
                        "rec_rmsd": (row.get(SRC["rec_rmsd"]) or "").strip(),
                        "eff_rmsd": (row.get(SRC["eff_rmsd"]) or "").strip(),
                    })

    rows = [best[k][1] for k in sorted(best)]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUT_COLS)
        writer.writeheader()
        writer.writerows(rows)

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"Read {n_rows} predictions from {n_files} benchmark(s).")
    print(f"Wrote {len(rows)} (model, msa, pdb) best-by-ra_eff rows -> {args.output}")
    if missing:
        print(f"NOTE: {len(missing)} benchmark dir(s) had no all_metrics.csv: "
              f"{', '.join(missing)}")
    if unknown_models:
        print(f"NOTE: unrecognised model names (mapped by suffix heuristic): "
              f"{', '.join(sorted(unknown_models))}")
    only_unscored = sorted(k for k in no_score if k not in best)
    if only_unscored:
        print(f"NOTE: {len(only_unscored)} (model, msa, pdb) combo(s) had no valid "
              f"ra_eff and were dropped:")
        for model, msa, pdb in only_unscored:
            print(f"        {pdb}  {model}  {msa}")


if __name__ == "__main__":
    main()

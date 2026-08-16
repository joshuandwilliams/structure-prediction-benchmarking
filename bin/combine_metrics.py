#!/usr/bin/env python3
"""Combine per-target all_metrics.csv files into one table.

For each (model, msa, pdb) it keeps the prediction with the highest avg_plddt,
the model's own most confident structure, and carries every metric column
through. Selection must be reference-free: picking the lowest RMSD needs the
answer and reports an oracle no pipeline can reach on a novel target.

The model column is split into a base model and an msa flag via MODEL_MAP.
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
    # Run 2: effector-template variants. These must be listed explicitly — the
    # suffix heuristic below keys off a trailing _msa/_nomsa, so
    # "boltz2_msa_template" would otherwise be mislabelled no_msa.
    "boltz2_template":             ("boltz2_template",             "no_msa"),
    "boltz2_msa_template":         ("boltz2_template",             "msa"),
    "boltz2_constrained_template": ("boltz2_constrained_template", "no_msa"),
}

RA_EFF_COL = "rmsd_effector_receptor_aligned"   # carried through, not the selection key
SELECT_COL = "avg_plddt"                        # selection key: HIGHEST wins
DERIVED = ["model", "msa", "pdb", "predictor"]  # prepended; "predictor" = orig "model"


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

    # best[(model, msa, pdb)] = (avg_plddt_value, full_source_row_dict)
    best = {}
    n_files = 0
    n_rows = 0
    src_fieldnames = None       # column order from the first all_metrics.csv read
    unknown_models = set()
    no_score = set()           # combos with no usable confidence or ra_eff
    missing = []               # benchmark dirs without an all_metrics.csv

    for pdb in sorted(os.listdir(bench_dir)):
        sub = os.path.join(bench_dir, pdb)
        if not os.path.isdir(sub):   # pragma: no cover
            continue
        csv_path = os.path.join(sub, f"{pdb}_benchmark_results", "all_metrics.csv")
        if not os.path.isfile(csv_path):
            missing.append(pdb)
            continue
        n_files += 1

        with open(csv_path, newline="") as fh:
            reader = csv.DictReader(fh)
            if src_fieldnames is None and reader.fieldnames:
                src_fieldnames = list(reader.fieldnames)
            for row in reader:
                n_rows += 1
                raw_model = (row.get("model") or "").strip()
                if not raw_model:   # pragma: no cover
                    continue
                if raw_model not in MODEL_MAP:   # pragma: no cover
                    unknown_models.add(raw_model)
                model, msa = map_model(raw_model)
                key = (model, msa, pdb)

                score = to_float(row.get(SELECT_COL))
                if score is None or to_float(row.get(RA_EFF_COL)) is None:
                    # No confidence to select on, or no ra_eff to score with.
                    no_score.add(key)
                    continue
                if key not in best or score > best[key][0]:
                    best[key] = (score, dict(row))   # keep the whole prediction

    if src_fieldnames is None:   # pragma: no cover
        src_fieldnames = []

    # Output header: derived cols + every source col, with the original "model"
    # column renamed to "predictor" (our derived "model" is the base name).
    out_fieldnames = list(DERIVED) + ["predictor" if c == "model" else c
                                      for c in src_fieldnames if c != "model"]

    rows = []
    for (model, msa, pdb) in sorted(best):
        src_row = best[(model, msa, pdb)][1]
        out_row = {"model": model, "msa": msa, "pdb": pdb,
                   "predictor": (src_row.get("model") or "").strip()}
        for c in src_fieldnames:
            if c == "model":
                continue
            out_row[c] = src_row.get(c, "")
        rows.append(out_row)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"Read {n_rows} predictions from {n_files} benchmark(s).")
    print(f"Wrote {len(rows)} (model, msa, pdb) best-by-confidence rows -> {args.output}")
    if missing:
        print(f"NOTE: {len(missing)} benchmark dir(s) had no all_metrics.csv: "
              f"{', '.join(missing)}")
    if unknown_models:   # pragma: no cover
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

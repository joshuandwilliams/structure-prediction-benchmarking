#!/usr/bin/env python3
"""
Fold the corrected constrained re-run into the committed analysis CSVs.

The published boltz1_constrained / boltz2_constrained rows were produced with
constraints written in the reference PDB's author numbering, which Boltz cannot
resolve (see bin/_constraint_geometry.read_ca_indexed).  Boltz-2 discarded them
outright; Boltz-1's were equally mis-targeted.  Those rows measure nothing and
must be replaced, not merged alongside.

This reads the re-run trees under
``experiments/benchmarks/<PDB>/rerun_constrained/<PDB>_constrained_rerun_results/``
and rewrites the constrained rows in both

    analysis/structure_prediction_benchmark/combined_metrics.csv
    analysis/structure_prediction_benchmark/per_prediction_metrics.csv

leaving every other model untouched.  Idempotent: re-running replaces the
constrained rows again rather than duplicating them.

Usage:
    merge_constrained_rerun.py [--dry-run]
"""

import argparse
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from combine_metrics import MODEL_MAP  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSIS = os.path.join(REPO, "analysis", "structure_prediction_benchmark")
RA = "rmsd_effector_receptor_aligned"
CONSTRAINED = ["boltz1_constrained", "boltz2_constrained"]
KEEP = ["model_name", "avg_plddt", "ptm", "iptm", "pae_mean", "ipae",
        "ipsae_ab", "ipsae_ba", "ipsae_min", "actifptm",
        "rmsd_receptor", "rmsd_effector_independent", RA]


def load_rerun():
    """Every prediction from the re-run trees, tagged (model, msa, pdb)."""
    rows = []
    pattern = os.path.join(REPO, "experiments", "benchmarks", "*",
                           "rerun_constrained", "*_constrained_rerun_results",
                           "all_metrics.csv")
    for path in sorted(glob.glob(pattern)):
        pdb = path.split(os.sep)[-4]
        d = pd.read_csv(path)
        d = d[d["model"].isin(CONSTRAINED)].copy()
        if d.empty:
            continue
        d[RA] = pd.to_numeric(d[RA], errors="coerce")
        d = d.dropna(subset=[RA])
        # AF3 emits an aggregate row; Boltz does not, but guard anyway.
        d = d[~d["model_name"].astype(str).str.contains("aggregate", case=False)]
        d["pdb"] = pdb
        d["predictor"] = d["model"]
        d["msa"] = d["model"].map(lambda m: MODEL_MAP[m][1])
        rows.append(d)
    if not rows:
        sys.exit("ERROR: no re-run results found. Pull them back first.")
    return pd.concat(rows, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    new = load_rerun()
    print(f"Re-run: {len(new)} predictions, {new['pdb'].nunique()} targets, "
          f"models {sorted(new['model'].unique())}")

    # ---- per_prediction_metrics.csv ----------------------------------------
    pp_path = os.path.join(ANALYSIS, "per_prediction_metrics.csv")
    pp = pd.read_csv(pp_path)
    before = len(pp)
    pp = pp[~pp["model"].isin(CONSTRAINED)]
    add = new[["model", "msa", "pdb", "predictor"] + KEEP]
    pp_out = pd.concat([pp, add], ignore_index=True)
    print(f"per_prediction_metrics.csv: {before} -> {len(pp_out)} "
          f"(dropped {before - len(pp)} stale constrained rows, added {len(add)})")

    # ---- combined_metrics.csv (highest avg_plddt per model/msa/pdb) --------
    # Must match the selection key in combine_metrics.py. Selecting these rows
    # on ra_eff while every other model is selected on confidence would leave
    # the constrained arms as an oracle inside an otherwise top-1 benchmark.
    cm_path = os.path.join(ANALYSIS, "combined_metrics.csv")
    cm = pd.read_csv(cm_path)
    cols = list(cm.columns)
    before = len(cm)
    cm = cm[~cm["model"].isin(CONSTRAINED)]
    best = (new.sort_values("avg_plddt", ascending=False)
               .groupby(["model", "msa", "pdb"], as_index=False)
               .first())
    best = best.reindex(columns=cols)
    cm_out = pd.concat([cm, best], ignore_index=True)
    cm_out = cm_out.sort_values(["model", "msa", "pdb"]).reset_index(drop=True)
    print(f"combined_metrics.csv:       {before} -> {len(cm_out)} "
          f"(dropped {before - len(cm)} stale, added {len(best)})")

    if args.dry_run:
        print("\nDRY RUN — nothing written.")
        return
    pp_out.to_csv(pp_path, index=False)
    cm_out.to_csv(cm_path, index=False)
    print("\nWritten. Re-render the analyses to pick the change up.")


if __name__ == "__main__":
    main()

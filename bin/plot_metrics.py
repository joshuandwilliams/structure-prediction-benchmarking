#!/usr/bin/env python3
"""
Plots from combined_metrics.csv (output of combine_metrics.py).

Generates:
  a) ra_eff_boxplot.png      — ra_eff distribution per model/msa/constraint combo,
                               No-MSA combos grouped left, MSA combos grouped right.
  b) success_rate_barplot.png — number of PDBs per combo with ra_eff < threshold
                               (default 5 Å) i.e. the "posed correctly" rate.

ra_eff = rmsd_effector_receptor_aligned (receptor-aligned effector RMSD; lower
is better). Each combo's distribution is across the PDBs (one best-by-ra_eff
prediction per PDB, as produced by combine_metrics.py).

Usage:
    plot_metrics.py [--csv FILE] [--outdir DIR] [--threshold 5.0]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

RA_CANDIDATES = ["ra_eff", "rmsd_effector_receptor_aligned"]
NO_MSA_COLOR = "#4C72B0"   # blue
MSA_COLOR    = "#C44E52"   # red


def pretty(model):
    return model.replace("_constrained", "\n(constr)")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    default_csv = os.path.normpath(
        os.path.join(here, "..", "experiments", "combined_metrics.csv"))

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=default_csv, help="combined_metrics.csv (default: %(default)s)")
    ap.add_argument("--outdir", default=None,
                    help="output dir for PNGs (default: alongside the csv, in plots/)")
    ap.add_argument("--threshold", type=float, default=5.0,
                    help="ra_eff success threshold in Angstrom (default: %(default)s)")
    args = ap.parse_args()

    if not os.path.isfile(args.csv):
        sys.exit(f"ERROR: csv not found: {args.csv}")
    outdir = args.outdir or os.path.join(os.path.dirname(os.path.abspath(args.csv)), "plots")
    os.makedirs(outdir, exist_ok=True)

    df = pd.read_csv(args.csv)
    ra_col = next((c for c in RA_CANDIDATES if c in df.columns), None)
    if ra_col is None:
        sys.exit(f"ERROR: no ra_eff column found (looked for {RA_CANDIDATES})")
    df = df[df["msa"].isin(["msa", "no_msa"])].copy()
    df[ra_col] = pd.to_numeric(df[ra_col], errors="coerce")
    df = df.dropna(subset=[ra_col])

    # Order combos: No-MSA block then MSA block, each sorted by median ra_eff.
    def ordered(msa_flag):
        sub = df[df["msa"] == msa_flag]
        meds = sub.groupby("model")[ra_col].median().sort_values()
        return list(meds.index)

    no_msa = [("no_msa", m) for m in ordered("no_msa")]
    msa    = [("msa",    m) for m in ordered("msa")]
    combos = no_msa + msa
    n_no = len(no_msa)

    series = {(f, m): df[(df["msa"] == f) & (df["model"] == m)][ra_col].values
              for (f, m) in combos}
    labels = [pretty(m) for (_, m) in combos]
    colors = [NO_MSA_COLOR] * n_no + [MSA_COLOR] * len(msa)

    # x positions with a gap between the two groups.
    xs = list(range(n_no)) + [i + 1 for i in range(n_no, len(combos))]
    gap_x = n_no + 0.0  # centre of the gap

    # ── (a) boxplot ───────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(max(10, 0.7 * len(combos) + 3), 6))
    data = [series[c] for c in combos]
    bp = ax.boxplot(data, positions=xs, widths=0.6, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black"))
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col); patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for x, c in zip(xs, combos):
        y = series[c]
        ax.scatter(x + rng.uniform(-0.15, 0.15, len(y)), y, s=14,
                   color="black", alpha=0.45, zorder=3)
    ax.axhline(args.threshold, ls="--", color="grey", lw=1)
    ax.text(xs[-1], args.threshold, f"  ra_eff = {args.threshold:g} Å",
            va="bottom", ha="right", color="grey", fontsize=9)
    ax.axvline(gap_x + 0.0, color="lightgrey", lw=1)
    ax.set_xticks(xs); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("ra_eff  (receptor-aligned effector RMSD, Å)")
    ax.set_title("Effector pose accuracy (ra_eff) per model / MSA / constraint")
    # group headers
    ax.text(np.mean(xs[:n_no]), ax.get_ylim()[1], "No MSA", ha="center",
            va="top", fontsize=11, fontweight="bold", color=NO_MSA_COLOR)
    ax.text(np.mean(xs[n_no:]), ax.get_ylim()[1], "MSA", ha="center",
            va="top", fontsize=11, fontweight="bold", color=MSA_COLOR)
    fig.tight_layout()
    box_path = os.path.join(outdir, "ra_eff_boxplot.png")
    fig.savefig(box_path, dpi=150); plt.close(fig)

    # ── (b) success-rate bar ──────────────────────────────────────────────
    counts = [int((series[c] < args.threshold).sum()) for c in combos]
    ns     = [len(series[c]) for c in combos]
    fig, ax = plt.subplots(figsize=(max(10, 0.7 * len(combos) + 3), 6))
    ax.bar(xs, counts, width=0.6, color=colors, alpha=0.85)
    for x, k, n in zip(xs, counts, ns):
        ax.text(x, k + 0.15, f"{k}/{n}", ha="center", va="bottom", fontsize=8)
    ax.axvline(gap_x, color="lightgrey", lw=1)
    ax.set_xticks(xs); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel(f"# PDBs with ra_eff < {args.threshold:g} Å")
    ax.set_ylim(0, max(ns) + 1)
    ax.set_title(f"Successful-pose rate (ra_eff < {args.threshold:g} Å) per model / MSA / constraint")
    ax.text(np.mean(xs[:n_no]), max(ns) + 0.6, "No MSA", ha="center",
            fontsize=11, fontweight="bold", color=NO_MSA_COLOR)
    ax.text(np.mean(xs[n_no:]), max(ns) + 0.6, "MSA", ha="center",
            fontsize=11, fontweight="bold", color=MSA_COLOR)
    fig.tight_layout()
    bar_path = os.path.join(outdir, "success_rate_barplot.png")
    fig.savefig(bar_path, dpi=150); plt.close(fig)

    # ── text summary ──────────────────────────────────────────────────────
    print(f"ra_eff column: {ra_col} | threshold: {args.threshold} Å")
    print(f"{'combo':<28} {'n':>3} {'median':>7} {'success':>8}")
    for c, k, n in zip(combos, counts, ns):
        f, m = c
        print(f"{m+' ['+f+']':<28} {n:>3} {np.median(series[c]):>7.2f} {k:>3}/{n}")
    print(f"\nWrote:\n  {box_path}\n  {bar_path}")


if __name__ == "__main__":
    main()

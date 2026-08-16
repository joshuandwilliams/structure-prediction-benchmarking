#!/usr/bin/env python3
"""Plot how many benchmark structures fall inside each model's training set.

A structure counts as exposed when its PDB initial release date is on or before
the model's published training cutoff. Release date rather than deposit date,
because the AlphaFold, Boltz and Chai pipelines filter the PDB by release.

Release dates come from a cached CSV, with any missing PDB fetched from the RCSB
API and appended.
"""

import argparse
import csv
import json
import os
import sys
import urllib.request

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# model -> (cutoff ISO date, one-line source). Cutoff = PDB *release* date filter.
# Keys are the full display names (no abbreviations).
MODEL_CUTOFFS = {
    "AlphaFold2-Multimer": ("2018-04-30", "AlphaFold2/-Multimer (Jumper 2021; Evans 2021)"),
    "ColabFold":           ("2018-04-30", "ColabFold = AF2/AF2-Multimer weights"),
    "Chai-1":              ("2021-01-12", "Chai-1 technical report (Chai Discovery 2024)"),
    "AlphaFold3":          ("2021-09-30", "Abramson et al. 2024 Nature (AlphaFold3)"),
    "Boltz-1":             ("2021-09-30", "Wohlwend et al. 2024 (same cutoff as AF3)"),
    "Boltz-2":             ("2023-06-01", "Boltz-2 report (biorxiv 2025.06.14.659707)"),
    # The benchmark runs Biohub's ESMFold2 (ESMC-6B backbone), pulled from the
    # plain `biohub/ESMFold2` HF repo = the production "cutoff2025" checkpoint
    # (NOT the original ESM-2 ESMFold, nor the experimental cutoff2021 variant).
    # The preprint's recent held-out PPI eval set starts 2025-06-30, so the
    # production model trains on PDB released before that date.
    "ESMFold2":            ("2025-06-30", "Biohub ESMFold2 cutoff2025 production model (biohub preprint 2026)"),
}

# Colourblind-safe pair (Okabe-Ito): orange = potential leakage, blue = clean.
IN_TRAIN_COLOR = "#D55E00"   # vermillion — released on/before cutoff
CLEAN_COLOR    = "#0072B2"   # blue       — released after cutoff (clean test case)


def fetch_release_date(pdb):
    """Query the RCSB data API for a PDB's initial release date (YYYY-MM-DD) or None."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb}"
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            d = json.load(r)
        ts = d.get("rcsb_accession_info", {}).get("initial_release_date", "")
        return ts.split("T")[0] or None
    except Exception as e:
        print(f"  WARN: could not fetch {pdb}: {e}", file=sys.stderr)
        return None


def load_dates(dates_csv, pdbs):
    """Return {pdb: pd.Timestamp}. Fetch+append any pdb missing from the cache."""
    cache = {}
    if os.path.isfile(dates_csv):
        with open(dates_csv, newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("release_date"):
                    cache[row["pdb"]] = row["release_date"]

    missing = [p for p in pdbs if p not in cache]
    if missing:
        print(f"Fetching {len(missing)} missing release date(s) from RCSB: "
              f"{', '.join(missing)}")
        for p in missing:
            d = fetch_release_date(p)
            if d:
                cache[p] = d
        # rewrite cache (sorted), preserving any deposit_date column if present
        with open(dates_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["pdb", "release_date"])
            for p in sorted(cache):
                w.writerow([p, cache[p]])

    return {p: pd.Timestamp(cache[p]) for p in pdbs if p in cache}


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.normpath(os.path.join(here, ".."))
    default_bench = os.path.join(repo, "experiments", "benchmarks")
    default_dates = os.path.join(repo, "data", "metrics", "pdb_release_dates.csv")

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmarks-dir", default=default_bench,
                    help="dir of <PDB>/ benchmark folders (default: %(default)s)")
    ap.add_argument("--dates-csv", default=default_dates,
                    help="cached pdb,release_date CSV (default: %(default)s)")
    ap.add_argument("--assessed-only", default=None,
                    help="restrict to PDBs present in this combined_metrics.csv "
                         "(the structures actually scored), instead of every "
                         "benchmark folder")
    ap.add_argument("--outdir", default=None,
                    help="output dir for the plot (default: alongside dates-csv, "
                         "in plots/)")
    args = ap.parse_args()

    if not os.path.isdir(args.benchmarks_dir):   # pragma: no cover
        sys.exit(f"ERROR: benchmarks dir not found: {args.benchmarks_dir}")
    outdir = args.outdir or os.path.join(os.path.dirname(os.path.abspath(args.dates_csv)),
                                         "plots")
    os.makedirs(outdir, exist_ok=True)

    pdbs = sorted(d for d in os.listdir(args.benchmarks_dir)
                  if os.path.isdir(os.path.join(args.benchmarks_dir, d)))
    scope = "all benchmark folders"
    if args.assessed_only:   # pragma: no cover
        cm = pd.read_csv(args.assessed_only)
        assessed = set(cm["pdb"].astype(str))
        pdbs = [p for p in pdbs if p in assessed]
        scope = "assessed structures only"

    dates = load_dates(args.dates_csv, pdbs)
    missing = [p for p in pdbs if p not in dates]
    if missing:   # pragma: no cover
        print(f"NOTE: no release date for {len(missing)} PDB(s), excluded: "
              f"{', '.join(missing)}", file=sys.stderr)
    rel = pd.Series({p: dates[p] for p in pdbs if p in dates})
    N = len(rel)

    models = sorted(MODEL_CUTOFFS, key=lambda m: MODEL_CUTOFFS[m][0])
    in_train = [int((rel <= pd.Timestamp(MODEL_CUTOFFS[m][0])).sum()) for m in models]
    clean = [N - k for k in in_train]

    # ── plot ────────────────────────────────────────────────────────────────
    x = np.arange(len(models))
    fig, ax = plt.subplots(figsize=(max(9, 1.1 * len(models) + 2), 6))
    ax.bar(x, in_train, color=IN_TRAIN_COLOR, label="in training set")
    ax.bar(x, clean, bottom=in_train, color=CLEAN_COLOR,
           label="after cutoff")
    for xi, k in zip(x, in_train):
        ax.text(xi, k / 2 if k else 0.3, f"{k}", ha="center", va="center",
                color="white", fontweight="bold", fontsize=12)
        ax.text(xi, N + 0.3, f"{k}/{N}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{m}\n(≤ {MODEL_CUTOFFS[m][0]})" for m in models],
                       fontsize=9)
    ax.set_ylabel("Number of benchmark structures")
    ax.set_ylim(0, N + 5)
    # Legend inside, above the bars (the stacked bars all reach N, leaving the
    # band above N free).
    ax.legend(loc="upper center", ncol=2, frameon=False, fontsize=10)

    base = os.path.join(outdir, "training_set_membership")
    for ext in ("png", "svg"):
        fig.savefig(f"{base}.{ext}", dpi=150, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    # ── text summary ──────────────────────────────────────────────────────
    print(f"\nScope: {scope}  (N = {N} structures with release dates)")
    print(f"{'model':<14} {'cutoff':<12} {'in_train':>8}  source")
    for m, k in zip(models, in_train):
        cut, src = MODEL_CUTOFFS[m]
        print(f"{m:<14} {cut:<12} {k:>5}/{N}  {src}")
    print(f"\nWrote:\n  {base}.png\n  {base}.svg")


if __name__ == "__main__":
    main()

"""Shared loading, constants and helpers for the analyses.

Each analysis has its own folder with its own thesis-figures/ and
supplementary-figures/. Quarto renders a document with the working directory set
to its folder, so save_fig writes relative to the CWD and figures land in the
right place without the module knowing which analysis called it.

The metric CSVs are shared inputs rather than any one analysis's output, so they
live in data/metrics/.

This module sits in bin/ because several analyses use it, and carries the
leading underscore that marks a bin/ module as a library rather than a CLI
entrypoint.

Usage from a .qmd:
    import sys; sys.path.insert(0, "../../bin")
    from _analysis_common import *
    d = load()
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

HERE = Path(__file__).resolve().parent          # bin/
REPO_ROOT = HERE.parent
DATA_DIR = REPO_ROOT / "data" / "metrics"       # shared metric CSVs
MANIFEST = REPO_ROOT / "data" / "benchmark_complexes.tsv"

TIER = 1              # small HMA / integrated-domain : effector pairs
RA_COL = "rmsd_effector_receptor_aligned"
RA_EFF_THRESHOLD = 5.0
CONF_THRESHOLD = 0.8

NO_MSA_COLOR = "#4C72B0"
MSA_COLOR = "#C44E52"
SYS_COLORS = {"Pik": "#4C72B0", "RGA5": "#DD8452"}

RA_TICKS = [0.5, 1, 2, 5, 10, 20, 40]      # Angstrom
TIME_TICKS = [2, 3, 5, 10, 20, 50, 100]    # minutes

DISPLAY = {
    "af2m": "AlphaFold2-Multimer", "af3": "AlphaFold3",
    "boltz1": "Boltz-1", "boltz1_constrained": "Boltz-1 (constr)",
    "boltz2": "Boltz-2", "boltz2_constrained": "Boltz-2 (constr)",
    "chai1": "Chai-1", "colabfold": "ColabFold", "esmfold2": "ESMFold2",
}
MSA_LABEL = {"msa": "MSA", "no_msa": "No MSA"}

# Published PDB training cutoffs (release-date filter).
MODEL_CUTOFFS = {
    "AlphaFold2-Multimer": "2018-04-30",
    "ColabFold": "2018-04-30",
    "Chai-1": "2021-01-12",
    "AlphaFold3": "2021-09-30",
    "Boltz-1": "2021-09-30",
    "Boltz-2": "2023-06-01",
    "ESMFold2": "2025-06-30",
}
# Per base model, for within-model splits.
MODEL_CUT = {"af2m": "2018-04-30", "colabfold": "2018-04-30",
             "chai1": "2021-01-12", "af3": "2021-09-30",
             "boltz1": "2021-09-30", "boltz1_constrained": "2021-09-30",
             "boltz2": "2023-06-01", "boltz2_constrained": "2023-06-01",
             "esmfold2": "2025-06-30"}

# The constrained variants derive their restraints FROM THE REFERENCE COMPLEX,
# so they cannot be run on a novel target.  They are not predictors: they are a
# positive control showing how well a model can be driven to a known interface.
# Ranking them alongside genuine predictors would be misleading, so they are
# flagged here and called out wherever they appear.
#
# Their earlier results were also unusable for a separate reason: constraints
# were written in the reference PDB's author numbering, which Boltz cannot
# resolve (Boltz renumbers its input 1..N).  Boltz-2 discarded them outright,
# making its constrained run byte identical to the unconstrained one; Boltz-1's
# were equally mis-targeted.  bin/_constraint_geometry.read_ca_indexed now emits
# Boltz token indices and both variants have been re-run, so the current rows
# are valid.
USES_REFERENCE = {"boltz1_constrained", "boltz2_constrained"}
INVALID_COMBOS = set()

# Chai-1 v0.6.x emits no PAE matrix, so everything derived from PAE is
# uncomputable for it, but compute_metrics.py writes 0.0 rather than a blank.
# Detected data-driven below so any future predictor with the same gap is caught.
PAE_DERIVED = ["pae_mean", "ipae", "ipsae_ab", "ipsae_ba", "ipsae_min"]
CONF_COLS = ["avg_plddt", "iptm", "ptm", "ipsae_min", "ranking_score"]

# Relative to the CWD, which Quarto sets to the rendering document's folder.
THESIS_FIGDIR = Path("thesis-figures")
SUPP_FIGDIR = Path("supplementary-figures")


def save_fig(fig, name, thesis=False, subdir=None):
    """Save a PNG. thesis=True routes to thesis-figures/, else supplementary."""
    out = THESIS_FIGDIR if thesis else SUPP_FIGDIR
    if subdir:
        out = out / subdir
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=200, bbox_inches="tight",
                pad_inches=0.05)


def pretty(model):
    """Display name, wrapped before a parenthesised variant to fit axis ticks."""
    return DISPLAY.get(model, model).replace(" (", "\n(")


def combo_label(model, msa_flag):
    return f"{pretty(model)} / {MSA_LABEL[msa_flag]}".replace("\n", " ")


def is_oracle(model):
    """Report whether this combo builds its input from the reference structure."""
    return model in USES_REFERENCE


def readable_log(axis, ticks):
    """Plain-number ticks on a log axis; the decade defaults are unreadable."""
    axis.set_ticks(ticks)
    axis.set_ticklabels([f"{v:g}" for v in ticks])
    axis.set_minor_locator(plt.NullLocator())


def auc(score, y):
    """ROC AUC, equivalently the Mann-Whitney U statistic normalised."""
    score, y = np.asarray(score, float), np.asarray(y, bool)
    m = ~np.isnan(score)
    score, y = score[m], y[m]
    if not y.any() or y.all():
        return np.nan
    u = mannwhitneyu(score[y], score[~y], alternative="two-sided").statistic
    return float(u / (y.sum() * (~y).sum()))


def spearman(a, b):
    return float(spearmanr(a, b).statistic)


def _clean(frame):
    """Shared cleaning: numeric coercion, pLDDT scale, blank uncomputable PAE."""
    for c in CONF_COLS + PAE_DERIVED:
        if c in frame.columns:
            frame[c] = pd.to_numeric(frame[c], errors="coerce")
    # Boltz reports pLDDT on 0-1, the AlphaFold lineage on 0-100.  Both are as
    # the tools emit them, but they cannot be compared without rescaling.
    frame["plddt"] = np.where(frame["avg_plddt"] <= 1.0,
                              frame["avg_plddt"] * 100.0, frame["avg_plddt"])
    for c in PAE_DERIVED:
        if c not in frame.columns:
            continue
        for m in frame["model"].unique():
            mask = frame["model"] == m
            if mask.any() and (frame.loc[mask, c] == 0).all():
                frame.loc[mask, c] = np.nan
    return frame


def load(verbose=True):
    """Load every input, cleaned and Tier-filtered. Returns a dict."""
    manifest = pd.read_csv(MANIFEST, sep="\t")
    tier = manifest[manifest["tier"] == TIER]
    tier_pdbs = set(tier["pdb"])

    df = pd.read_csv(DATA_DIR / "combined_metrics.csv")
    df = df[df["msa"].isin(["msa", "no_msa"])].copy()
    df[RA_COL] = pd.to_numeric(df[RA_COL], errors="coerce")
    df = df.dropna(subset=[RA_COL])
    dropped = sorted(set(df["pdb"]) - tier_pdbs)
    df = df[df["pdb"].isin(tier_pdbs)].copy()
    df = df.merge(manifest[["pdb", "system", "description"]], on="pdb", how="left")
    df = _clean(df)

    pred = pd.read_csv(DATA_DIR / "per_prediction_metrics.csv")
    pred = pred[pred["pdb"].isin(tier_pdbs)].copy()
    pred[RA_COL] = pd.to_numeric(pred[RA_COL], errors="coerce")
    pred["rmsd_receptor"] = pd.to_numeric(pred["rmsd_receptor"], errors="coerce")
    pred = pred.dropna(subset=[RA_COL])
    pred = _clean(pred)
    pred["correct"] = pred[RA_COL] < RA_EFF_THRESHOLD

    rt = pd.read_csv(DATA_DIR / "runtime_stats.csv")
    rt = rt[rt["pdb"].isin(tier_pdbs)].copy()
    for c in ["standalone_elapsed_s", "peak_rss_gb"]:
        rt[c] = pd.to_numeric(rt[c], errors="coerce")
    rt = rt.dropna(subset=["standalone_elapsed_s", "peak_rss_gb"])
    # standalone_elapsed_s already folds in the shared COLABFOLD_SEARCH MSA
    # step for the models that depend on it (boltz1_msa, boltz2_msa,
    # colabfold); it equals elapsed_s for every other model.
    rt["elapsed_min"] = rt["standalone_elapsed_s"] / 60.0

    for m, f in INVALID_COMBOS:   # pragma: no cover
        df = df[~((df["model"] == m) & (df["msa"] == f))]
        pred = pred[~((pred["model"] == m) & (pred["msa"] == f))]
        rt = rt[~((rt["model"] == m) & (rt["msa"] == f))]

    sim = pd.read_csv(DATA_DIR / "reference_similarity.csv")
    near = sim[(sim["receptor_identity"] >= 95) &
               (sim["effector_identity"] >= 95)].copy()

    rel = pd.read_csv(DATA_DIR / "pdb_release_dates.csv", parse_dates=["release_date"])
    rel = rel[rel["pdb"].isin(tier_pdbs)]

    combos = ordered_combos(df)
    lookup = {(r["model"], r["msa"], r["pdb"]): r for _, r in df.iterrows()}

    if verbose:
        print(f"Tier {TIER}: {len(tier_pdbs)} complexes, "
              f"{df['pdb'].nunique()} with results, {len(combos)} model/MSA combos")
        if dropped:
            print(f"  excluded (not tier {TIER}): {dropped}")
        if INVALID_COMBOS:   # pragma: no cover
            print(f"  excluded: {sorted(m for m, _ in INVALID_COMBOS)}")
        oracle = sorted(set(df["model"]) & USES_REFERENCE)
        if oracle:
            print(f"  reference-derived (NOT predictors, see USES_REFERENCE): "
                  f"{oracle}")
        blank = sorted({m for m in df['model'].unique()
                        if df.loc[df['model'] == m, 'ipsae_min'].isna().all()})
        if blank:
            print(f"  PAE-derived metrics not computable for: {blank}")
        print(f"  per-prediction rows: {len(pred)};  runtime rows: {len(rt)};  "
              f"near-identical pairs: {len(near)}")

    return dict(manifest=manifest, tier=tier, tier_pdbs=tier_pdbs, df=df,
                pred=pred, rt=rt, sim=sim, near=near, rel=rel,
                combos=combos, lookup=lookup)


def ordered_combos(df):
    """No-MSA block then MSA block, each ranked by correct-pose count."""
    def ordered(flag):
        g = df[df["msa"] == flag].groupby("model")[RA_COL]
        stat = pd.DataFrame({
            "n_ok": g.apply(lambda v: int((v < RA_EFF_THRESHOLD).sum())),
            "median": g.median(),
        })
        return list(stat.sort_values(["n_ok", "median"],
                                     ascending=[False, True]).index)
    return ([("no_msa", m) for m in ordered("no_msa")] +
            [("msa", m) for m in ordered("msa")])


def block_positions(combos, gap=1.0):
    """X positions with a gap between the No-MSA and MSA blocks."""
    n_no = sum(1 for f, _ in combos if f == "no_msa")
    xs = list(range(n_no)) + [i + gap for i in range(n_no, len(combos))]
    return xs, n_no


def msa_legend(ax, **kw):
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec="black")
               for c in (NO_MSA_COLOR, MSA_COLOR)]
    defaults = dict(loc="center left", bbox_to_anchor=(1.01, 0.5),
                    frameon=True, fontsize=10)
    defaults.update(kw)
    return ax.legend(handles, ["No MSA", "MSA"], **defaults)

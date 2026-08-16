#!/usr/bin/env python3
"""Classify which receptor surface each prediction placed the effector on.

The Pik HMA has two effector-binding surfaces in this set, the AVR-Pik surface
and the AVR-Pia surface of 6Q76, about 27 Å apart. Each structure is reduced to
the set of receptor residues its effector contacts, expressed in a common frame,
and compared by Jaccard overlap against both. Contacts rather than centroid
distance, because a centroid test ignores orientation and passes an effector
sitting on the right surface the wrong way round.

Reads best_models/, which holds each model's confidence-selected structure
rather than its lowest-RMSD one, so the counts are what a model returns
unprompted rather than its achievable ceiling.

Usage:
    classify_predicted_interface.py --plot-dir DIR --output interfaces.csv
"""

import argparse
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import matplotlib

matplotlib.use("Agg")
import _structure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from combine_metrics import MODEL_MAP  # noqa: E402
from common_wrong_interface import (  # noqa: E402
    aligner,
    apply_tf,
    kabsch,
    matched_indices,
)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FRAME = "6Q76"

DISPLAY = {
    "af2m": "AlphaFold2-Multimer", "af3": "AlphaFold3",
    "boltz1": "Boltz-1", "boltz1_constrained": "Boltz-1 (constr)",
    "boltz2": "Boltz-2", "boltz2_constrained": "Boltz-2 (constr)",
    "chai1": "Chai-1", "colabfold": "ColabFold", "esmfold2": "ESMFold2",
}


ca_and_seq = _structure.read_chain


def contact_set(rec_xyz, eff_xyz, rec_index_to_frame, cutoff):
    """Receptor residues contacting the effector, in COMMON-FRAME numbering.

    rec_xyz / eff_xyz are (N,3) Ca arrays; rec_index_to_frame maps a receptor
    row index to its position in the frame receptor, so sets from different
    targets are directly comparable.
    """
    if len(rec_xyz) == 0 or len(eff_xyz) == 0:
        return set()
    d = np.linalg.norm(rec_xyz[:, None, :] - eff_xyz[None, :, :], axis=-1)
    hit = np.where((d <= cutoff).any(axis=1))[0]
    return {rec_index_to_frame[i] for i in hit if i in rec_index_to_frame}


def jaccard(a, b):
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def find_pred_flat(flat_dir, pdb, predictor):
    for ext in (".pdb", ".cif"):
        p = os.path.join(flat_dir, f"{pdb}__{predictor}{ext}")
        if os.path.isfile(p):
            return p
    return None


def find_pred(pdb, predictor, base_model):
    d = os.path.join(REPO, "experiments", "benchmarks", pdb,
                     f"{pdb}_benchmark_results", "best_models",
                     f"{predictor}_best")
    for ext in (".pdb", ".cif"):
        p = os.path.join(d, f"{base_model}_best{ext}")
        if os.path.isfile(p):
            return p
    if os.path.isdir(d):                       # fall back to whatever is there
        for f in sorted(os.listdir(d)):
            if f.endswith((".pdb", ".cif")):
                return os.path.join(d, f)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--contact-cutoff", type=float, default=8.0,
                    help="Ca-Ca distance defining an interface contact (Å)")
    ap.add_argument("--min-jaccard", type=float, default=0.10,
                    help="below this against both interfaces, label 'other'")
    ap.add_argument("--output", default="predicted_interfaces.csv")
    ap.add_argument("--plot-dir", default=None)
    ap.add_argument("--pred-dir", default=None,
                    help="flat dir of <PDB>__<tag>.{pdb,cif} files to use "
                         "instead of the published best_models/ trees")
    args = ap.parse_args()

    man = pd.read_csv(os.path.join(REPO, "data", "benchmark_complexes.tsv"), sep="\t")
    man = man[(man["tier"] == 1) & (man["system"] == "Pik")]
    targets = sorted(man["pdb"])

    def ref_path(p):
        return os.path.join(REPO, "data", "complexes_for_benchmarking", f"{p}.pdb")

    al = aligner()
    frame_rec, frame_rec_seq = ca_and_seq(ref_path(FRAME), "A")
    frame_eff, _ = ca_and_seq(ref_path(FRAME), "B")

    # For each target: the TRUE AVR-Pik interface as a set of receptor residues
    # in the frame's numbering, plus its own reference (for recomputing ra_eff).
    def frame_map(rec_seq):
        """Row index in this receptor -> row index in the frame receptor."""
        return {i: j for i, j in matched_indices(al, rec_seq, frame_rec_seq)}

    true_iface, own_ref = {}, {}
    for pdb in targets:
        rec, rec_seq = ca_and_seq(ref_path(pdb), "A")
        eff, eff_seq = ca_and_seq(ref_path(pdb), "B")
        true_iface[pdb] = contact_set(rec, eff, frame_map(rec_seq),
                                      args.contact_cutoff)
        own_ref[pdb] = (rec, rec_seq, eff, eff_seq)

    # The AVR-Pia interface, from the frame target itself.
    pia_iface = contact_set(frame_rec, frame_eff,
                            {i: i for i in range(len(frame_rec))},
                            args.contact_cutoff)

    # The two interfaces SHARE residues -- about a third of the AVR-Pia patch in
    # the Pik targets. A pose touching only the shared region would score partly
    # on both and the label would turn on noise. Compare against the
    # DISCRIMINATING residues only, so the assignment depends solely on what
    # distinguishes the two surfaces.
    shared = {r for p, v in true_iface.items() if p != FRAME for r in (v & pia_iface)}
    pia_only = pia_iface - shared
    print(f"shared between the two interfaces: {len(shared)} residues "
          f"{sorted(shared)} -- excluded from both sets")
    print(f"discriminating: AVR-Pia {len(pia_only)} residues; AVR-Pik "
          f"{min(len(true_iface[p] - shared) for p in targets if p != FRAME)}-"
          f"{max(len(true_iface[p] - shared) for p in targets if p != FRAME)}")
    print(f"AVR-Pia interface: {len(pia_iface)} receptor residues; "
          f"true AVR-Pik interfaces: "
          f"{min(len(v) for k, v in true_iface.items() if k != FRAME)}-"
          f"{max(len(v) for k, v in true_iface.items() if k != FRAME)} residues")
    ov = [jaccard(v, pia_iface) for k, v in true_iface.items() if k != FRAME]
    print(f"raw overlap between the two true interfaces: Jaccard "
          f"{min(ov):.2f}-{max(ov):.2f}")

    # AlphaFold2-Multimer does not always emit chains A/B (it wrote B/C for four
    # targets), so use the chain IDs compute_metrics recorded for the chosen
    # prediction rather than assuming.
    chain_map = {}
    if args.pred_dir:
        cpath = os.path.join(args.pred_dir, "_chains.csv")
        if os.path.isfile(cpath):
            cm = pd.read_csv(cpath)
            chain_map = {(r["pdb"], r["tag"]): (r["rec_chain"], r["eff_chain"])
                         for _, r in cm.iterrows()}

    rows = []
    for predictor, (model, msa) in MODEL_MAP.items():
        for pdb in targets:
            pp = (find_pred_flat(args.pred_dir, pdb, predictor)
                  if args.pred_dir else find_pred(pdb, predictor, model))
            if pp is None:
                continue
            rec_ch, eff_ch = chain_map.get((pdb, predictor), ("A", "B"))
            try:
                rec, rec_seq = ca_and_seq(pp, rec_ch)
                eff, eff_seq = ca_and_seq(pp, eff_ch)
            except Exception as e:
                print(f"  skip {predictor}/{pdb}: {e}")
                continue
            m = matched_indices(al, rec_seq, frame_rec_seq)
            if len(m) < 20:
                print(f"  skip {predictor}/{pdb}: only {len(m)} aligned residues")
                continue
            # ra_eff OF THIS STRUCTURE.  combined_metrics.csv holds the
            # lowest-ra_eff of the 25 predictions, whereas this classifies the
            # model's own confidence-selected pick -- a different structure.
            r_rec, r_rec_seq, r_eff, r_eff_seq = own_ref[pdb]
            mo = matched_indices(al, rec_seq, r_rec_seq)
            tf_own = kabsch(rec[[i for i, _ in mo]], r_rec[[j for _, j in mo]])
            eff_own = apply_tf(eff, tf_own)
            me = matched_indices(al, eff_seq, r_eff_seq)
            ra_self = float(np.sqrt(
                ((eff_own[[i for i, _ in me]] - r_eff[[j for _, j in me]]) ** 2)
                .sum(1).mean())) if me else float("nan")

            # Contact-based assignment.  No superposition is involved: the
            # contact set is a property of the predicted complex itself, so
            # orientation is captured and frame choice cannot bias it.
            pred_iface = contact_set(rec, eff, frame_map(rec_seq),
                                     args.contact_cutoff)
            # Residues common to both true interfaces carry no information about
            # which surface was used, so they are removed from all three sets.
            pred_d = pred_iface - shared
            j_own = jaccard(pred_d, true_iface[pdb] - shared)
            j_pia = jaccard(pred_d, pia_only)
            if max(j_own, j_pia) < args.min_jaccard:
                site = "other"
            elif j_own >= j_pia:
                site = "correct site"
            else:
                site = "AVR-Pia site"

            rows.append({"predictor": predictor, "model": model, "msa": msa,
                         "pdb": pdb,
                         "jaccard_true": round(j_own, 3),
                         "jaccard_pia": round(j_pia, 3),
                         "n_contacts": len(pred_iface),
                         "n_contacts_discriminating": len(pred_d),
                         "n_pik_only": len(pred_iface & (true_iface[pdb] - shared)),
                         "n_shared": len(pred_iface & shared),
                         "n_pia_only": len(pred_iface & pia_only),
                         # denominators, so coverage can be computed downstream
                         "tot_pik_only": len(true_iface[pdb] - shared),
                         "tot_shared": len(shared),
                         "tot_pia_only": len(pia_only),
                         "ra_eff_this_structure": round(ra_self, 2),
                         "site": site})

    out = pd.DataFrame(rows)
    out.to_csv(args.output, index=False)

    # 6Q76 defines the frame and its correct site IS the AVR-Pia site.
    tally = out[out["pdb"] != FRAME]
    n_targets = tally["pdb"].nunique()
    counts = (tally.groupby(["model", "msa", "site"]).size()
                   .unstack(fill_value=0)
                   .reindex(columns=["correct site", "AVR-Pia site", "other"],
                            fill_value=0))
    counts["combo"] = [f"{DISPLAY.get(m, m)} / {'MSA' if s == 'msa' else 'No MSA'}"
                       for m, s in counts.index]
    counts = counts.sort_values(["correct site", "AVR-Pia site"],
                                ascending=[False, False])
    print(f"\nInterface used, over {n_targets} Pik targets ({FRAME} excluded):")
    print(counts.set_index("combo").to_string())

    if args.plot_dir:
        os.makedirs(args.plot_dir, exist_ok=True)
        fig, ax = plt.subplots(figsize=(10.5, 6.2))
        x = np.arange(len(counts))
        bottom = np.zeros(len(counts))
        for col, colour in [("correct site", "#0072B2"),
                            ("AVR-Pia site", "#D55E00"),
                            ("other", "#BBBBBB")]:
            ax.bar(x, counts[col], bottom=bottom, color=colour,
                   edgecolor="black", linewidth=0.6, label=col)
            bottom += counts[col].values
        ax.set_xticks(x)
        ax.set_xticklabels(counts["combo"], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(f"Pik targets (of {n_targets})")
        ax.set_title("Which HMA surface did each model place the effector on?",
                     fontsize=11)
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=True,
                  fontsize=9)
        fig.savefig(os.path.join(args.plot_dir, "interface_choice_by_model.png"),
                    dpi=200, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        print(f"\nWrote plot to {args.plot_dir}")


if __name__ == "__main__":
    main()

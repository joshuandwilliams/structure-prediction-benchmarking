#!/usr/bin/env python3
"""
Classify which interface each model predicted, across every model.

The Pik-family HMA has (at least) two effector-binding surfaces represented in
this benchmark: the AVR-Pik surface used by every Pik target except 6Q76, and a
distinct surface used by AVR-Pia in 6Q76.  The two sit about 27 A apart, so a
prediction can be assigned to one or the other by where its effector lands once
all receptors are put in a common frame.

Assignment is by INTERFACE CONTACTS, not by centroid distance.  A centroid test
answers the wrong question: it ignores orientation entirely, so an effector
sitting on the correct surface but rotated the wrong way scores as correct.
Boltz-1 restrained exposed this — every target passed a 13 A centroid test while
RMSDs ran to 22 A.

Instead, each structure is reduced to the SET of receptor residues its effector
contacts (any Ca pair within --contact-cutoff A), expressed in the common frame's
numbering so sets are comparable across targets.  That set is compared by
Jaccard overlap against
  (a) the target's own true AVR-Pik interface, and
  (b) the AVR-Pia interface from 6Q76,
and the prediction is labelled:

    correct site   higher Jaccard with its own true interface, above --min-jaccard
    AVR-Pia site   higher Jaccard with the AVR-Pia interface, above --min-jaccard
    other          below --min-jaccard against both

6Q76 is EXCLUDED from the tally: it defines the frame, and its correct site is
the AVR-Pia site, so it cannot distinguish the two categories.

Caveat.  This reads best_models/, which holds each model's own confident pick
rather than its lowest-RMSD pose (compute_metrics.py selects the highest
avg_plddt; earlier runs selected the highest actifpTM).  That is the right basis
for asking what a model would hand you unprompted, but it is not the model's
best achievable pose, so the counts here are not an accuracy ceiling.

Usage:
    classify_predicted_interface.py --plot-dir <dir> --output interfaces.csv
"""

import argparse
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from combine_metrics import MODEL_MAP  # noqa: E402
from test_common_wrong_interface import (  # noqa: E402
    aligner, apply_tf, kabsch, matched_indices,
)

from Bio.PDB import MMCIFParser, PDBParser  # noqa: E402
from Bio.PDB.MMCIF2Dict import MMCIF2Dict  # noqa: E402
from Bio.PDB.Polypeptide import protein_letters_3to1  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FRAME = "6Q76"

DISPLAY = {
    "af2m": "AlphaFold2-Multimer", "af3": "AlphaFold3",
    "boltz1": "Boltz-1", "boltz1_constrained": "Boltz-1 (constr)",
    "boltz2": "Boltz-2", "boltz2_constrained": "Boltz-2 (constr)",
    "chai1": "Chai-1", "colabfold": "ColabFold", "esmfold2": "ESMFold2",
}


def _ca_from_cif_dict(path, chain_id):
    """Minimal mmCIF CA reader.

    Biopython's MMCIFParser requires _atom_site.occupancy, which ESMFold2 does
    not emit, so fall back to reading the atom_site loop directly.
    """
    d = MMCIF2Dict(path)
    chains = d["_atom_site.auth_asym_id"] if "_atom_site.auth_asym_id" in d \
        else d["_atom_site.label_asym_id"]
    names = d["_atom_site.label_atom_id"]
    comps = d["_atom_site.label_comp_id"]
    xs, ys, zs = (d["_atom_site.Cartn_x"], d["_atom_site.Cartn_y"],
                  d["_atom_site.Cartn_z"])
    coords, seq = [], []
    for ch, nm, comp, x, y, z in zip(chains, names, comps, xs, ys, zs):
        if ch != chain_id or nm != "CA":
            continue
        try:
            seq.append(protein_letters_3to1[comp])
        except KeyError:
            continue
        coords.append((float(x), float(y), float(z)))
    return np.asarray(coords, float), "".join(seq)


def ca_and_seq(path, chain_id):
    """Ca coordinates and sequence for one chain, from PDB or mmCIF."""
    if path.endswith(".cif"):
        try:
            st = MMCIFParser(QUIET=True).get_structure("x", path)
        except Exception:
            return _ca_from_cif_dict(path, chain_id)
    else:
        st = PDBParser(QUIET=True).get_structure("x", path)
    coords, seq = [], []
    for res in st[0][chain_id]:
        if "CA" not in res:
            continue
        try:
            seq.append(protein_letters_3to1[res.get_resname()])
        except KeyError:
            continue
        coords.append(res["CA"].get_coord())
    return np.asarray(coords, float), "".join(seq)


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
        """row index in this receptor -> row index in the frame receptor."""
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
    shared_all = set.intersection(*[true_iface[p] for p in targets if p != FRAME]) \
        if len(targets) > 1 else set()
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
            tf = kabsch(rec[[i for i, _ in m]], frame_rec[[j for _, j in m]])
            c = apply_tf(eff, tf).mean(0)

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

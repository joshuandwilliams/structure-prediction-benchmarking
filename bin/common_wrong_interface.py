#!/usr/bin/env python3
"""Test whether a predictor's failed poses land on the same wrong interface.

For Boltz-2 without an MSA most Pik targets fail at 19 to 27 Å while 6Q76
(Pikp1-HMA with AVR-Pia) is solved to 0.6 Å. AVR-Pia engages a different HMA
surface, so the hypothesis is that the model places AVR-Pik effectors there,
giving one preferred binding mode that is right for 6Q76 and wrong for the rest.

Every receptor is a close homologue, so all structures are put in one frame and
each predicted effector centroid is compared against its own true position and
against the AVR-Pia position.

Usage:
    common_wrong_interface.py --model boltz2 --msa no_msa --system Pik
"""

import argparse
import importlib.util
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

if importlib.util.find_spec("Bio") is None:
    sys.exit("ERROR: biopython required. pip install biopython")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


ca_and_seq = _structure.read_chain


def aligner():
    return _structure.aligner()


def matched_indices(al, seq_a, seq_b):
    """Index pairs of aligned, non-gap positions. ``al`` is accepted and ignored."""
    return list(zip(*_structure.matched_indices(seq_a, seq_b)))


def kabsch(mobile, target):
    """Rotation and translation mapping mobile onto target."""
    _, R, t, _ = _structure.kabsch(mobile, target)
    return R, t


def apply_tf(coords, tf):
    return _structure.apply_transform(coords, *tf)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="boltz2")
    ap.add_argument("--msa", default="no_msa", choices=["msa", "no_msa"])
    ap.add_argument("--frame", default="6Q76",
                    help="target whose reference receptor defines the frame")
    ap.add_argument("--system", default="Pik",
                    help="restrict to one system, or 'all'")
    ap.add_argument("--tier", type=int, default=1)
    ap.add_argument("--output", default="common_interface_test.csv")
    ap.add_argument("--plot-dir", default=None,
                    help="write the two diagnostic plots here")
    ap.add_argument("--pred-dir", default=None,
                    help="flat dir of <PDB>__<tag>.{pdb,cif} files to use "
                         "instead of the published best_models/ trees")
    args = ap.parse_args()

    man = pd.read_csv(os.path.join(REPO, "data", "benchmark_complexes.tsv"), sep="\t")
    man = man[man["tier"] == args.tier]
    if args.system != "all":
        man = man[man["system"] == args.system]
    targets = sorted(man["pdb"])
    if args.frame not in targets:
        sys.exit(f"ERROR: frame target {args.frame} not in the selected set")


    def ref_path(p):
        return os.path.join(REPO, "data", "complexes_for_benchmarking", f"{p}.pdb")

    # best_models/ is keyed by the PREDICTOR tag (boltz2, boltz2_msa, ...),
    # while the file inside is named after the base model, so both are needed.
    sys.path.insert(0, os.path.join(REPO, "bin"))
    from combine_metrics import MODEL_MAP
    tags = [k for k, v in MODEL_MAP.items() if v == (args.model, args.msa)]
    if not tags:
        sys.exit(f"ERROR: no predictor tag for {args.model}/{args.msa}")
    predictor_tag = tags[0]

    def pred_path(p):
        if args.pred_dir:
            for ext in (".pdb", ".cif"):
                flat = os.path.join(args.pred_dir, f"{p}__{predictor_tag}{ext}")
                if os.path.isfile(flat):
                    return flat
            return os.path.join(args.pred_dir, f"{p}__{predictor_tag}.pdb")
        return os.path.join(REPO, "experiments", "benchmarks", p,
                            f"{p}_benchmark_results", "best_models",
                            f"{predictor_tag}_best", f"{args.model}_best.pdb")

    al = aligner()
    frame_rec, frame_rec_seq = ca_and_seq(ref_path(args.frame), "A")
    frame_eff, _ = ca_and_seq(ref_path(args.frame), "B")
    frame_eff_centroid = frame_eff.mean(0)

    rows, clouds = [], {}
    for pdb in targets:
        rp, pp = ref_path(pdb), pred_path(pdb)
        if not (os.path.isfile(rp) and os.path.isfile(pp)):
            print(f"  skip {pdb}: missing structure")
            continue

        ref_rec, ref_rec_seq = ca_and_seq(rp, "A")
        ref_eff, ref_eff_seq = ca_and_seq(rp, "B")
        pred_rec, pred_rec_seq = ca_and_seq(pp, "A")
        pred_eff, pred_eff_seq = ca_and_seq(pp, "B")

        # Both receptors -> common frame, effectors carried along.
        pr = matched_indices(al, ref_rec_seq, frame_rec_seq)
        tf_ref = kabsch(ref_rec[[i for i, _ in pr]], frame_rec[[j for _, j in pr]])
        pp_ = matched_indices(al, pred_rec_seq, frame_rec_seq)
        tf_pred = kabsch(pred_rec[[i for i, _ in pp_]], frame_rec[[j for _, j in pp_]])

        true_eff_f = apply_tf(ref_eff, tf_ref)
        pred_eff_f = apply_tf(pred_eff, tf_pred)
        clouds[pdb] = pred_eff_f

        # ra_eff is recomputed here rather than joined from
        # combined_metrics.csv.  That file describes the lowest-ra_eff of the 25
        # predictions, whereas best_models/ holds the model's own
        # confidence-selected pick, which is a different structure (11th of 25
        # by RMSD for 6FUD, 4th of 25 for 6G11).  Joining the two would pair a
        # metric with a structure it did not come from, so the RMSD for the
        # structure actually analysed is computed directly.
        tf_self = kabsch(pred_rec[[i for i, _ in matched_indices(al, pred_rec_seq, ref_rec_seq)]],
                         ref_rec[[j for _, j in matched_indices(al, pred_rec_seq, ref_rec_seq)]])
        pe_self = apply_tf(pred_eff, tf_self)
        me = matched_indices(al, pred_eff_seq, ref_eff_seq)
        ra_self = float(np.sqrt(
            ((pe_self[[i for i, _ in me]] - ref_eff[[j for _, j in me]]) ** 2)
            .sum(1).mean()))

        c_pred, c_true = pred_eff_f.mean(0), true_eff_f.mean(0)
        rows.append({
            "pdb": pdb,
            "ra_eff": ra_self,
            "d_pred_to_own_true": float(np.linalg.norm(c_pred - c_true)),
            "d_pred_to_frame_eff": float(np.linalg.norm(c_pred - frame_eff_centroid)),
            "d_true_to_frame_eff": float(np.linalg.norm(c_true - frame_eff_centroid)),
        })

    out = pd.DataFrame(rows).sort_values("ra_eff")

    # Pairwise distance between predicted effector centroids.
    keys = list(out["pdb"])
    cents = np.array([clouds[k].mean(0) for k in keys])
    pair = pd.DataFrame(
        np.linalg.norm(cents[:, None, :] - cents[None, :, :], axis=-1),
        index=keys, columns=keys)

    if args.plot_dir:
        os.makedirs(args.plot_dir, exist_ok=True)
        tag = f"{args.model}_{args.msa}"
        ok = out["ra_eff"] < 5.0
        frame_row = out["pdb"] == args.frame
        typical_true = float(out.loc[~frame_row, "d_true_to_frame_eff"].median())

        # (1) how wrong is the pose, against how close it is to the frame site
        fig, ax = plt.subplots(figsize=(8.2, 6.4))
        ax.axhline(typical_true, ls="--", color="grey", lw=1)
        # Annotate on the left; the legend occupies the top right.
        ax.text(0.01, typical_true, f" true AVR-Pik site sits {typical_true:.0f} Å "
                f"from the {args.frame} site", transform=ax.get_yaxis_transform(),
                ha="left", va="bottom", fontsize=8, color="grey")
        for mask, col, lab in [(ok & ~frame_row, "#0072B2", "pose correct (<5 Å)"),
                               (~ok & ~frame_row, "#D55E00", "pose incorrect")]:
            s = out[mask]
            ax.scatter(s["d_pred_to_own_true"], s["d_pred_to_frame_eff"], s=70,
                       color=col, edgecolor="black", linewidth=0.5, label=lab,
                       zorder=3)
        s = out[frame_row]
        ax.scatter(s["d_pred_to_own_true"], s["d_pred_to_frame_eff"], s=110,
                   marker="*", color="black", label=f"{args.frame} (defines frame)",
                   zorder=4)
        for _, r in out.iterrows():
            ax.annotate(r["pdb"], (r["d_pred_to_own_true"], r["d_pred_to_frame_eff"]),
                        fontsize=7, xytext=(6, 0), textcoords="offset points",
                        va="center", zorder=5)
        ax.set_xlabel("Predicted effector to its OWN true position (Å)")
        ax.set_ylabel(f"Predicted effector to the {args.frame} effector position (Å)")
        ax.set_title(f"{args.model} / {args.msa}: do failed poses land on the "
                     f"{args.frame} interface?", fontsize=10)
        ax.legend(frameon=False, loc="upper right", fontsize=9)
        fig.savefig(os.path.join(args.plot_dir, f"wrong_interface_scatter_{tag}.png"),
                    dpi=200, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)

        # (2) do the predicted effectors cluster with each other
        fig, ax = plt.subplots(figsize=(8.4, 7.0))
        im = ax.imshow(pair.values, cmap="viridis_r")
        ax.set_xticks(range(len(keys)))
        ax.set_xticklabels(keys, rotation=90, fontsize=8)
        ax.set_yticks(range(len(keys)))
        ax.set_yticklabels([f"{k}  ({v:.0f} Å)" for k, v in
                            zip(keys, out["ra_eff"])], fontsize=8)
        for i in range(len(keys)):
            for j in range(len(keys)):
                v = pair.values[i, j]
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6,
                        color="white" if v > pair.values.max() * 0.55 else "black")
        ax.set_title("Distance between predicted effector centroids (Å)\n"
                     "rows ordered by pose error, shown in brackets", fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.8, label="centroid separation (Å)")
        fig.savefig(os.path.join(args.plot_dir, f"wrong_interface_clustering_{tag}.png"),
                    dpi=200, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        print(f"\nWrote 2 plots to {args.plot_dir}")

    out.to_csv(args.output, index=False)
    pair.to_csv(args.output.replace(".csv", "_pairwise.csv"))
    print(out.round(2).to_string(index=False))
    print(f"\nWrote {args.output} and its _pairwise.csv "
          f"({len(keys)} targets, frame = {args.frame})")


if __name__ == "__main__":
    main()

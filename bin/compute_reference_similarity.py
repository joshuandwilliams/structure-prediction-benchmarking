#!/usr/bin/env python3
"""
Pairwise sequence similarity between benchmark reference complexes.

Tier 1 contains allelic effector series (AVR-PikD/E/A/C/F) and engineered HMA
variants, so many references differ by only one or two residues.  That makes it
possible to ask whether near-identical inputs are predicted with near-identical
accuracy, or whether small sequence changes produce large accuracy swings.

For every pair of references this emits the receptor and effector percent
identity, the substitution count, and the number of alignment gap positions.
Gap positions are reported separately because the references are crystal
constructs: two entries can have the same underlying protein but different
modelled boundaries, which is a construct difference rather than a sequence
difference.

Identity is computed over the shorter sequence, so it is not deflated by one
construct simply being longer than the other.

Requires biopython.  Output is small (n*(n-1)/2 rows) and meant to be committed
so the analysis renders without re-running this.

Usage:
    compute_reference_similarity.py --tier 1 --output reference_similarity.csv
"""

import argparse
import csv
import itertools
import os
import sys
import warnings

warnings.filterwarnings("ignore")

try:
    from Bio import Align
    from Bio.Align import substitution_matrices
    from Bio.PDB import PDBParser, PPBuilder
except ImportError:
    sys.exit("ERROR: biopython required.  pip install biopython")

RECEPTOR_CHAIN = "A"    # plant protein, by the data/ two-chain convention
EFFECTOR_CHAIN = "B"    # pathogen effector


def build_aligner():
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    al.mode = "global"
    return al


def chain_sequences(pdb_path):
    parser = PDBParser(QUIET=True)
    ppb = PPBuilder()
    struct = parser.get_structure("x", pdb_path)
    out = {}
    for chain in struct[0]:
        out[chain.id] = "".join(str(pp.get_sequence())
                                for pp in ppb.build_peptides(chain))
    return out


def compare(aligner, a, b):
    """Return (pct_identity_over_shorter, n_substitutions, n_gap_positions)."""
    if not a or not b:
        return None, None, None
    aln = aligner.align(a, b)[0]
    ga, gb = aln[0], aln[1]
    matches = subs = gaps = 0
    for x, y in zip(ga, gb):
        if x == "-" or y == "-":
            gaps += 1
        elif x == y:
            matches += 1
        else:
            subs += 1
    return 100.0 * matches / min(len(a), len(b)), subs, gaps


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs-dir", default=os.path.normpath(
        os.path.join(here, "..", "data", "complexes_for_benchmarking")))
    ap.add_argument("--manifest", default=os.path.normpath(
        os.path.join(here, "..", "data", "benchmark_complexes.tsv")))
    ap.add_argument("--tier", type=int, default=None)
    ap.add_argument("--output", default="reference_similarity.csv")
    args = ap.parse_args()

    with open(args.manifest, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if args.tier is not None:
        rows = [r for r in rows if int(r["tier"]) == args.tier]
    pdbs = sorted(r["pdb"] for r in rows)
    system = {r["pdb"]: r["system"] for r in rows}
    if len(pdbs) < 2:
        sys.exit("ERROR: need at least two references to compare.")

    seqs = {}
    for pdb in pdbs:
        path = os.path.join(args.refs_dir, f"{pdb}.pdb")
        if not os.path.isfile(path):
            print(f"  WARNING: missing reference {path}")
            continue
        seqs[pdb] = chain_sequences(path)
    print(f"Read {len(seqs)} references")

    aligner = build_aligner()
    out = []
    for a, b in itertools.combinations(sorted(seqs), 2):
        rec_id, rec_sub, rec_gap = compare(
            aligner, seqs[a].get(RECEPTOR_CHAIN, ""), seqs[b].get(RECEPTOR_CHAIN, ""))
        eff_id, eff_sub, eff_gap = compare(
            aligner, seqs[a].get(EFFECTOR_CHAIN, ""), seqs[b].get(EFFECTOR_CHAIN, ""))
        if rec_id is None or eff_id is None:
            continue
        out.append({
            "pdb_a": a, "pdb_b": b,
            "system_a": system.get(a, ""), "system_b": system.get(b, ""),
            "same_system": system.get(a) == system.get(b),
            "receptor_identity": round(rec_id, 2),
            "receptor_substitutions": rec_sub,
            "receptor_gap_positions": rec_gap,
            "effector_identity": round(eff_id, 2),
            "effector_substitutions": eff_sub,
            "effector_gap_positions": eff_gap,
            "min_identity": round(min(rec_id, eff_id), 2),
            "total_substitutions": rec_sub + eff_sub,
        })

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"Wrote {len(out)} pairs -> {args.output}")


if __name__ == "__main__":
    main()

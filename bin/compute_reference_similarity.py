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

Alignment is done by EMBOSS ``needle`` (Needleman-Wunsch, global), called as a
subprocess with ``-aformat3 fasta`` so the gapped sequences come back directly.
Using the standard tool rather than a hand-rolled aligner means the parameters
are the published EMBOSS defaults (EBLOSUM62, gap open 10.0, gap extend 0.5) and
the alignment is reproducible by anyone with EMBOSS.

Percent identity has no single agreed definition, so **both** conventional
denominators are emitted and neither is left implicit:

* ``*_identity`` divides by the shorter ungapped sequence.  A construct is then
  not penalised for residues the other entry simply did not model, which is what
  we want here, since several references are the same protein resolved to
  different boundaries.
* ``*_identity_over_alignment`` divides by the full alignment length, including
  gap columns.  This is the number ``needle`` itself prints.

Requires biopython (sequence extraction from coordinates) and EMBOSS on PATH.
Output is small (n*(n-1)/2 rows) and meant to be committed so the analysis
renders without re-running this.

Usage:
    compute_reference_similarity.py --tier 1 --output reference_similarity.csv
"""

import argparse
import csv
import itertools
import os
import shutil
import subprocess
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")

try:
    from Bio.PDB import PDBParser, PPBuilder
except ImportError:
    sys.exit("ERROR: biopython required.  pip install biopython")

RECEPTOR_CHAIN = "A"    # plant protein, by the data/ two-chain convention
EFFECTOR_CHAIN = "B"    # pathogen effector

# EMBOSS needle defaults, stated explicitly so the run is reproducible.
GAP_OPEN = "10.0"
GAP_EXTEND = "0.5"
MATRIX = "EBLOSUM62"


def needle_align(a, b, workdir):
    """Global-align two sequences with EMBOSS needle. Returns the gapped pair."""
    fa = os.path.join(workdir, "a.fasta")
    fb = os.path.join(workdir, "b.fasta")
    out = os.path.join(workdir, "aln.fasta")
    with open(fa, "w") as fh:
        fh.write(f">a\n{a}\n")
    with open(fb, "w") as fh:
        fh.write(f">b\n{b}\n")
    subprocess.run(
        ["needle", "-asequence", fa, "-bsequence", fb,
         "-gapopen", GAP_OPEN, "-gapextend", GAP_EXTEND,
         "-datafile", MATRIX, "-aformat3", "fasta",
         "-outfile", out, "-auto"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

    seqs, cur = [], []
    with open(out) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    if len(seqs) != 2 or len(seqs[0]) != len(seqs[1]):
        raise RuntimeError("unexpected needle output")
    return seqs[0].upper(), seqs[1].upper()


def chain_sequences(pdb_path):
    parser = PDBParser(QUIET=True)
    ppb = PPBuilder()
    struct = parser.get_structure("x", pdb_path)
    out = {}
    for chain in struct[0]:
        out[chain.id] = "".join(str(pp.get_sequence())
                                for pp in ppb.build_peptides(chain))
    return out


def compare(workdir, a, b):
    """Identity by both denominators, plus substitution and gap counts."""
    if not a or not b:
        return None, None, None, None
    ga, gb = needle_align(a, b, workdir)
    matches = subs = gaps = 0
    for x, y in zip(ga, gb):
        if x == "-" or y == "-":
            gaps += 1
        elif x == y:
            matches += 1
        else:
            subs += 1
    return (100.0 * matches / min(len(a), len(b)),
            100.0 * matches / len(ga),
            subs, gaps)


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

    if shutil.which("needle") is None:
        sys.exit("ERROR: EMBOSS needle not on PATH.  "
                 "conda install -c bioconda emboss")

    out = []
    with tempfile.TemporaryDirectory() as workdir:
        for a, b in itertools.combinations(sorted(seqs), 2):
            rec = compare(workdir, seqs[a].get(RECEPTOR_CHAIN, ""),
                          seqs[b].get(RECEPTOR_CHAIN, ""))
            eff = compare(workdir, seqs[a].get(EFFECTOR_CHAIN, ""),
                          seqs[b].get(EFFECTOR_CHAIN, ""))
            if rec[0] is None or eff[0] is None:
                continue
            rec_id, rec_id_aln, rec_sub, rec_gap = rec
            eff_id, eff_id_aln, eff_sub, eff_gap = eff
            out.append({
                "pdb_a": a, "pdb_b": b,
                "system_a": system.get(a, ""), "system_b": system.get(b, ""),
                "same_system": system.get(a) == system.get(b),
                "receptor_identity": round(rec_id, 2),
                "receptor_identity_over_alignment": round(rec_id_aln, 2),
                "receptor_substitutions": rec_sub,
                "receptor_gap_positions": rec_gap,
                "effector_identity": round(eff_id, 2),
                "effector_identity_over_alignment": round(eff_id_aln, 2),
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

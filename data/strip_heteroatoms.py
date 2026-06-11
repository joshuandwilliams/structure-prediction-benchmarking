#!/usr/bin/env python3
"""
strip_heteroatoms.py
--------------------
Remove non-polymer records (HETATM) and supplementary records (ANISOU, CONECT)
from the two-chain reference PDBs in complexes_for_benchmarking/.

Why
===
The references are extracted from experimental structures that carry
crystallisation additives, ions, and cofactors (HOH, SO4, EDO, ADP, ATP, ZN,
CL, MN, CA, MPD, …) on the receptor/target chains. Those heteroatoms break
downstream tools — notably the negative-steering engine's effector-template
extraction, which fails with `_atom_site.label_seq_id is '.'` when a HETATM
(e.g. SO4) sits in the effector chain. Only the protein matters for the
benchmark (sequence, fold, interface, RMSD), so we strip the rest.

Safety
======
None of the benchmark PDBs contain modified amino acids (no MSE etc.) — every
HETATM is a water/additive/ion/cofactor — so removing all HETATM never deletes
a polymer residue. The script asserts the chain-A and chain-B Cα counts are
unchanged before overwriting each file.

This is invoked automatically by extract_benchmark_complexes.py after it writes
each complex, and can also be run standalone:

    python strip_heteroatoms.py                      # clean complexes_for_benchmarking/
    python strip_heteroatoms.py path/to/file.pdb …   # clean specific files
"""

from __future__ import annotations

import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DIR = os.path.join(HERE, "complexes_for_benchmarking")

_DROP_PREFIXES = ("HETATM", "ANISOU", "CONECT")


def _ca_counts(lines):
    """{chain_id: n_CA} for ATOM CA records, by author chain column (22)."""
    counts = {}
    for ln in lines:
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            counts[ln[21]] = counts.get(ln[21], 0) + 1
    return counts


def strip_file(path: str) -> bool:
    """Strip HETATM/ANISOU/CONECT from one PDB in place. Returns True if changed."""
    lines = open(path).read().splitlines()
    if not any(ln.startswith("HETATM") for ln in lines):
        return False
    before = _ca_counts(lines)
    kept = [ln for ln in lines if not ln.startswith(_DROP_PREFIXES)]
    after = _ca_counts(kept)
    assert before == after, f"{path}: Cα counts changed {before} -> {after} (aborting)"
    with open(path, "w") as fh:
        fh.write("\n".join(kept) + "\n")
    return True


def main(argv):
    targets = argv or sorted(glob.glob(os.path.join(DEFAULT_DIR, "*.pdb")))
    if not targets:
        print(f"No PDBs found (looked in {DEFAULT_DIR}).", file=sys.stderr)
        return 1
    changed = 0
    for path in targets:
        if strip_file(path):
            changed += 1
            print(f"stripped {os.path.basename(path)}")
    print(f"\nDone: {changed} of {len(targets)} file(s) had heteroatoms removed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

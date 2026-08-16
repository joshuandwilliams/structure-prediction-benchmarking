#!/usr/bin/env python3
"""Remove HETATM, ANISOU and CONECT records from the reference PDBs.

Waters, ions and cofactors are not part of the complex being scored, and the
structure-negative-steering engine fails on heteroatoms in the effector chain.
No modified residues are present in this set, so no polymer is lost.

Called automatically by extract_benchmark_complexes.py, and runnable standalone
to re-clean existing files.
"""

from __future__ import annotations

import glob
import os
import sys

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
DEFAULT_DIR = os.path.join(DATA, "complexes_for_benchmarking")

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

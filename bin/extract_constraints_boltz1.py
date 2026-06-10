#!/usr/bin/env python3
"""
extract_constraints_boltz1.py
-----------------------------
Extract a pocket-only constraint from a two-chain reference PDB for Boltz-1.

Boltz-1 schema.py enforces max_distance == 6.0 for ALL constraint types,
so contact constraints (which have variable distances) are NOT generated.
Only the pocket constraint is written.

Usage:
    python extract_constraints_boltz1.py <pdb> <rec_chain> <eff_chain> \
        <pocket_cutoff> <pocket_max_distance>

Writes YAML-formatted constraint block to stdout; prints stats to stderr.
"""
import sys

from _constraint_geometry import (
    format_pocket_block,
    pocket_residues,
    read_ca_by_chain,
)


def main():
    if len(sys.argv) != 6:
        print("Usage: extract_constraints_boltz1.py PDB REC EFF POCKET_CUTOFF POCKET_MAX_D",
              file=sys.stderr)
        sys.exit(2)

    pdb_path      = sys.argv[1]
    rec_chain     = sys.argv[2]
    eff_chain     = sys.argv[3]
    pocket_cutoff = float(sys.argv[4])
    pocket_max_d  = float(sys.argv[5])

    coords = read_ca_by_chain(pdb_path, chains={rec_chain, eff_chain})
    rec_ca = coords.get(rec_chain, {})
    eff_ca = coords.get(eff_chain, {})

    if not rec_ca:
        print(f"ERROR: No Cα atoms for chain {rec_chain}", file=sys.stderr)
        sys.exit(1)
    if not eff_ca:
        print(f"ERROR: No Cα atoms for chain {eff_chain}", file=sys.stderr)
        sys.exit(1)

    print(f"Receptor Cα atoms: {len(rec_ca)}", file=sys.stderr)
    print(f"Effector Cα atoms: {len(eff_ca)}", file=sys.stderr)

    pocket = pocket_residues(rec_ca, eff_ca, pocket_cutoff)
    print(f"Pocket residues:   {len(pocket)} (within {pocket_cutoff} Å)",
          file=sys.stderr)
    if pocket:
        print(f"  Range: {min(pocket)}-{max(pocket)}", file=sys.stderr)

    if not pocket:
        print("ERROR: No pocket residues found. Check chain IDs.", file=sys.stderr)
        sys.exit(1)

    print("constraints:")
    print(format_pocket_block(rec_chain, eff_chain, pocket, pocket_max_d))


if __name__ == "__main__":
    main()

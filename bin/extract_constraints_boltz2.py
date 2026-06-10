#!/usr/bin/env python3
"""
extract_constraints_boltz2.py
-----------------------------
Extract pocket + contact constraints from a two-chain reference PDB for Boltz-2.

Boltz-2 (unlike Boltz-1) supports arbitrary max_distance values for both
pocket and contact constraints, enabling dense residue-pair restraints.

Usage:
    python extract_constraints_boltz2.py <pdb> <rec_chain> <eff_chain> \
        <contact_cutoff> <contact_max> <contact_tolerance> \
        <pocket_cutoff> <pocket_max_distance>

Writes YAML-formatted constraint block to stdout; prints stats to stderr.
"""
import sys

from _constraint_geometry import (
    contact_pairs,
    format_contact_block,
    format_pocket_block,
    pocket_residues,
    read_ca_by_chain,
)


def main():
    if len(sys.argv) != 9:
        print("Usage: extract_constraints_boltz2.py PDB REC EFF "
              "CONTACT_CUTOFF CONTACT_MAX CONTACT_TOL POCKET_CUTOFF POCKET_MAX_D",
              file=sys.stderr)
        sys.exit(2)

    pdb_path       = sys.argv[1]
    rec_chain      = sys.argv[2]
    eff_chain      = sys.argv[3]
    contact_cutoff = float(sys.argv[4])
    contact_max    = int(sys.argv[5])
    contact_tol    = float(sys.argv[6])
    pocket_cutoff  = float(sys.argv[7])
    pocket_max_d   = float(sys.argv[8])

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

    contacts = contact_pairs(rec_ca, eff_ca, contact_cutoff, contact_max)
    print(f"Contact pairs:     {len(contacts)} (within {contact_cutoff} Å, capped at {contact_max})",
          file=sys.stderr)
    if contacts:
        print(f"  Closest: {rec_chain}{contacts[0][0]}-{eff_chain}{contacts[0][1]}: {contacts[0][2]:.1f} A",
              file=sys.stderr)
        print(f"  Widest:  {rec_chain}{contacts[-1][0]}-{eff_chain}{contacts[-1][1]}: {contacts[-1][2]:.1f} A",
              file=sys.stderr)

    if not pocket and not contacts:
        print("ERROR: No constraints generated. Check chain IDs.", file=sys.stderr)
        sys.exit(1)

    print("constraints:")

    if pocket:
        print(format_pocket_block(rec_chain, eff_chain, pocket, pocket_max_d))

    for r_rn, e_rn, d in contacts:
        max_d = round(d + contact_tol, 1)
        print(format_contact_block(rec_chain, eff_chain, r_rn, e_rn, max_d))


if __name__ == "__main__":
    main()

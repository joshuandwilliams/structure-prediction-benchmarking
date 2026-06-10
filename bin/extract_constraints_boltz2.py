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

    rec_ca = {}
    eff_ca = {}

    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                ch     = line[21]
                resnum = int(line[22:26].strip())
                xyz    = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                if ch == rec_chain:
                    rec_ca[resnum] = xyz
                elif ch == eff_chain:
                    eff_ca[resnum] = xyz

    if not rec_ca:
        print(f"ERROR: No Cα atoms for chain {rec_chain}", file=sys.stderr)
        sys.exit(1)
    if not eff_ca:
        print(f"ERROR: No Cα atoms for chain {eff_chain}", file=sys.stderr)
        sys.exit(1)

    print(f"Receptor Cα atoms: {len(rec_ca)}", file=sys.stderr)
    print(f"Effector Cα atoms: {len(eff_ca)}", file=sys.stderr)

    pocket_residues = set()
    for r_rn, r_xyz in rec_ca.items():
        for e_xyz in eff_ca.values():
            d = sum((a - b) ** 2 for a, b in zip(r_xyz, e_xyz)) ** 0.5
            if d <= pocket_cutoff:
                pocket_residues.add(r_rn)
                break

    pocket_residues = sorted(pocket_residues)
    print(f"Pocket residues:   {len(pocket_residues)} (within {pocket_cutoff} Å)",
          file=sys.stderr)
    if pocket_residues:
        print(f"  Range: {min(pocket_residues)}-{max(pocket_residues)}", file=sys.stderr)

    contacts = []
    for r_rn, r_xyz in rec_ca.items():
        for e_rn, e_xyz in eff_ca.items():
            d = sum((a - b) ** 2 for a, b in zip(r_xyz, e_xyz)) ** 0.5
            if d <= contact_cutoff:
                contacts.append((r_rn, e_rn, d))

    contacts.sort(key=lambda x: x[2])
    contacts = contacts[:contact_max]

    print(f"Contact pairs:     {len(contacts)} (within {contact_cutoff} Å, capped at {contact_max})",
          file=sys.stderr)
    if contacts:
        print(f"  Closest: {rec_chain}{contacts[0][0]}-{eff_chain}{contacts[0][1]}: {contacts[0][2]:.1f} A",
              file=sys.stderr)
        print(f"  Widest:  {rec_chain}{contacts[-1][0]}-{eff_chain}{contacts[-1][1]}: {contacts[-1][2]:.1f} A",
              file=sys.stderr)

    if not pocket_residues and not contacts:
        print("ERROR: No constraints generated. Check chain IDs.", file=sys.stderr)
        sys.exit(1)

    print("constraints:")

    if pocket_residues:
        contacts_str = ", ".join(f"[{rec_chain}, {r}]" for r in pocket_residues)
        print(f"  - pocket:")
        print(f"      binder: {eff_chain}")
        print(f"      contacts: [{contacts_str}]")
        print(f"      max_distance: {pocket_max_d}")
        print(f"      force: true")

    for r_rn, e_rn, d in contacts:
        max_d = round(d + contact_tol, 1)
        print(f"  - contact:")
        print(f"      token1: [{rec_chain}, {r_rn}]")
        print(f"      token2: [{eff_chain}, {e_rn}]")
        print(f"      max_distance: {max_d}")
        print(f"      force: true")


if __name__ == "__main__":
    main()

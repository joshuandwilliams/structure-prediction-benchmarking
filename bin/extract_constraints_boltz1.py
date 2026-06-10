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

    if not pocket_residues:
        print("ERROR: No pocket residues found. Check chain IDs.", file=sys.stderr)
        sys.exit(1)

    contacts_str = ", ".join(f"[{rec_chain}, {r}]" for r in pocket_residues)
    print("constraints:")
    print(f"  - pocket:")
    print(f"      binder: {eff_chain}")
    print(f"      contacts: [{contacts_str}]")
    print(f"      max_distance: {pocket_max_d}")
    print(f"      force: true")


if __name__ == "__main__":
    main()

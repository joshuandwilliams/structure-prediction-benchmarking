#!/usr/bin/env python3
"""
validate_boltz_yaml.py
----------------------
Validate a Boltz input YAML before running `boltz predict`.

Checks:
  - sequences and constraints parse
  - For --model boltz1: max_distance must be exactly 6.0 on every constraint
    block, and no `contact` constraints are allowed (schema rejects them)
  - For --model boltz2: token1 and token2 on contact blocks are well-formed
    [chain, resnum] pairs

Exits non-zero on any structural problem.

Usage:
    python validate_boltz_yaml.py <input.yaml> --model {boltz1,boltz2}
"""
import argparse
import sys

import yaml


def main():
    p = argparse.ArgumentParser()
    p.add_argument("yaml_path")
    p.add_argument("--model", choices=["boltz1", "boltz2"], required=True)
    args = p.parse_args()

    with open(args.yaml_path) as f:
        data = yaml.safe_load(f)

    seqs        = data.get("sequences", []) or []
    constraints = data.get("constraints", []) or []
    n_pocket    = sum(1 for c in constraints if "pocket" in c)
    n_contact   = sum(1 for c in constraints if "contact" in c)
    print(f"  Sequences: {len(seqs)}, Constraints: {n_pocket} pocket + {n_contact} contact")

    if args.model == "boltz1":
        # Boltz-1 schema rejects contact constraints outright
        if n_contact > 0:
            print("  ERROR: Contact constraints present — Boltz-1 will reject them!")
            sys.exit(1)
        for c in constraints:
            key = "pocket" if "pocket" in c else "contact" if "contact" in c else None
            if key:
                md = c[key].get("max_distance")
                if md != 6.0:
                    print(f"  ERROR: {key} max_distance={md} — Boltz-1 requires exactly 6.0!")
                    sys.exit(1)
    else:  # boltz2
        for c in constraints:
            if "contact" in c:
                ct = c["contact"]
                if not isinstance(ct.get("token1"), list) or len(ct["token1"]) != 2:
                    print(f"  ERROR: malformed token1: {ct.get('token1')}")
                    sys.exit(1)
                if not isinstance(ct.get("token2"), list) or len(ct["token2"]) != 2:
                    print(f"  ERROR: malformed token2: {ct.get('token2')}")
                    sys.exit(1)

    print("  YAML structure OK")


if __name__ == "__main__":
    main()

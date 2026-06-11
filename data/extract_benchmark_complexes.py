#!/usr/bin/env python3
"""
extract_benchmark_complexes.py
------------------------------
Build two-chain binder-target reference complexes for the benchmark.

Reads data/benchmark_complexes.tsv. For each entry it loads the downloaded
structure from solved_NLR_structures/ (.cif preferred, .pdb fallback), chooses
the receptor-chain / target-chain pair in closest contact (so a multi-copy
crystal still yields a genuinely bound pair), relabels:

    receptor -> chain A      (always the PLANT protein)
    target   -> chain B      (always the PATHOGEN effector)

and writes a clean 2-chain PDB to complexes_for_benchmarking/<PDB>.pdb.

The A = plant / B = pathogen convention is uniform across all tiers. For Tier-3
effector/host-target complexes the plant host protein is still chain A even
though biologically the effector is the "binder" there — uniform labelling
keeps every benchmark case trackable and lets one params.yml
(receptor_chain: A, effector_chain: B) drive them all.

Requires gemmi:  pip install gemmi   (or run inside the benchmark container).

Usage:
    python extract_benchmark_complexes.py            # all entries in the manifest
    python extract_benchmark_complexes.py 6G10 9QLV  # only these PDB IDs
"""

from __future__ import annotations

import csv
import os
import sys

import gemmi

import strip_heteroatoms

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "solved_NLR_structures")
OUT = os.path.join(HERE, "complexes_for_benchmarking")
MANIFEST = os.path.join(HERE, "benchmark_complexes.tsv")
CONTACT_CUTOFF = 8.0  # Å, Cα–Cα


def load_structure(pdb_id):
    for ext in (".cif", ".pdb"):
        path = os.path.join(SRC, pdb_id + ext)
        if os.path.exists(path):
            st = gemmi.read_structure(path)
            st.setup_entities()
            return st
    return None


def ca_positions(chain):
    out = []
    for res in chain:
        atom = res.find_atom("CA", "*")
        if atom is not None:
            out.append(atom.pos)
    return out


def count_contacts(ca_a, ca_b, cutoff=CONTACT_CUTOFF):
    n = 0
    for p in ca_a:
        for q in ca_b:
            if p.dist(q) <= cutoff:
                n += 1
    return n


def pick_bound_pair(model, rec_chains, tgt_chains):
    """Return ((rec_name, tgt_name), n_contacts) for the receptor/target chain
    pair with the most Cα–Cα contacts, or (None, 0) if no chains are found."""
    wanted = set(rec_chains) | set(tgt_chains)
    ca = {ch.name: ca_positions(ch) for ch in model if ch.name in wanted}
    best, best_n = None, -1
    for r in rec_chains:
        for t in tgt_chains:
            if r in ca and t in ca:
                n = count_contacts(ca[r], ca[t])
                if n > best_n:
                    best, best_n = (r, t), n
    return best, max(best_n, 0)


def build_two_chain(st, rec_name, tgt_name):
    """New Structure containing only rec_name -> A and tgt_name -> B."""
    model = st[0]
    new = gemmi.Structure()
    new.name = st.name
    new.cell = st.cell
    new.spacegroup_hm = st.spacegroup_hm
    nm = gemmi.Model("1")
    for src, dst in ((rec_name, "A"), (tgt_name, "B")):
        nc = gemmi.Chain(dst)
        for res in model[src]:
            nc.add_residue(res)
        nm.add_chain(nc)
    new.add_model(nm)
    new.setup_entities()
    return new


def main():
    only = {a.strip().upper() for a in sys.argv[1:]}
    os.makedirs(OUT, exist_ok=True)
    with open(MANIFEST) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    written = warned = missing = 0
    for row in rows:
        pid = row["pdb"].strip().upper()
        if only and pid not in only:
            continue
        rec_chains = [c.strip() for c in row["receptor_chains"].split(",") if c.strip()]
        tgt_chains = [c.strip() for c in row["target_chains"].split(",") if c.strip()]

        st = load_structure(pid)
        if st is None:
            print(f"[{pid}] MISSING source — run download_solved_structures.sh first",
                  file=sys.stderr)
            missing += 1
            continue

        (pair, n) = pick_bound_pair(st[0], rec_chains, tgt_chains)
        if pair is None:
            print(f"[{pid}] ERROR: none of receptor {rec_chains} / target "
                  f"{tgt_chains} chains found in structure", file=sys.stderr)
            missing += 1
            continue

        rec, tgt = pair
        if n == 0:
            print(f"[{pid}] WARNING: chosen {rec}/{tgt} have no Cα contacts "
                  f"(<{CONTACT_CUTOFF} Å) — check the chain assignment", file=sys.stderr)
            warned += 1

        new = build_two_chain(st, rec, tgt)
        out_path = os.path.join(OUT, f"{pid}.pdb")
        new.write_pdb(out_path)
        # Drop heteroatoms (waters/additives/ions/cofactors copied through from
        # the source chains) so the reference is protein-only — see
        # strip_heteroatoms.py for the why.
        strip_heteroatoms.strip_file(out_path)
        print(f"[{pid}] tier{row['tier']:>2} {row['system']:<9} "
              f"{rec}->A {tgt}->B  ({n} contacts)  -> {os.path.basename(out_path)}")
        written += 1

    print(f"\nDone: {written} written, {warned} contact-warnings, {missing} missing/failed.",
          file=sys.stderr)
    if warned or missing:
        sys.exit(1)


if __name__ == "__main__":
    main()

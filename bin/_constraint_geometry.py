"""Cα geometry and YAML formatting shared by the two Boltz constraint extractors.

Boltz-1 emits a pocket block only, Boltz-2 adds dense per-pair contacts.
Everything else is common.
"""

from __future__ import annotations

import numpy as np
from _structure import contact_pairs as _contact_pairs
from _structure import read_ca_by_chain  # noqa: F401  (re-exported)


def read_ca_indexed(pdb_path, chains=None):
    """Cα coordinates keyed by Boltz token index rather than author numbering.

    Boltz renumbers whatever sequence it receives as 1..N per chain, so a
    constraint in author numbering targets the wrong residue or none at all.
    9IMU's receptor is authored 996-1070 but reaches Boltz as 1-75.

    Indices follow file order over the polymer, matching the sequence
    extract_sequences.py hands to Boltz. Residues are enumerated from every
    ATOM record, so a residue with no Cα still consumes an index. Keying on
    (resseq, icode) keeps insertion codes from colliding with a bare integer.

    Returns (coords, labels) where coords is {chain: {index: (x, y, z)}} and
    labels is {chain: {index: "1065"}} for logging.
    """
    wanted = set(chains) if chains is not None else None
    order: dict[str, list[tuple[int, str]]] = {}
    ca: dict[str, dict[tuple[int, str], tuple[float, float, float]]] = {}

    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            ch = line[21]
            if wanted is not None and ch not in wanted:
                continue
            key = (int(line[22:26].strip()), line[26].strip())
            seen = order.setdefault(ch, [])
            if not seen or seen[-1] != key:
                if key not in seen:
                    seen.append(key)
            if line[12:16].strip() == "CA":
                ca.setdefault(ch, {})[key] = (
                    float(line[30:38]), float(line[38:46]), float(line[46:54]),
                )

    coords, labels = {}, {}
    for ch, keys in order.items():
        idx_of = {k: i + 1 for i, k in enumerate(keys)}
        labels[ch] = {i + 1: f"{k[0]}{k[1]}" for i, k in enumerate(keys)}
        coords[ch] = {idx_of[k]: xyz for k, xyz in ca.get(ch, {}).items()}
    return coords, labels


def distance(a, b):
    return float(np.linalg.norm(np.asarray(a, float) - np.asarray(b, float)))


def pocket_residues(rec_ca, eff_ca, cutoff):
    """Receptor residue numbers with a Cα within cutoff of any effector Cα."""
    if not rec_ca or not eff_ca:
        return []
    rec_ids, eff_ids = list(rec_ca), list(eff_ca)
    hits = _contact_pairs([rec_ca[i] for i in rec_ids],
                          [eff_ca[j] for j in eff_ids], cutoff)
    return sorted({rec_ids[i] for i, _, _ in hits})


def contact_pairs(rec_ca, eff_ca, cutoff, max_pairs):
    """Closest inter-chain Cα residue pairs, as (rec, eff, dist), nearest first."""
    if not rec_ca or not eff_ca:
        return []
    rec_ids, eff_ids = list(rec_ca), list(eff_ca)
    hits = _contact_pairs([rec_ca[i] for i in rec_ids],
                          [eff_ca[j] for j in eff_ids], cutoff)
    return [(rec_ids[i], eff_ids[j], d) for i, j, d in hits[:max_pairs]]


def format_pocket_block(rec_chain, eff_chain, residues, max_distance):
    contacts_str = ", ".join(f"[{rec_chain}, {r}]" for r in residues)
    return "\n".join([
        "  - pocket:",
        f"      binder: {eff_chain}",
        f"      contacts: [{contacts_str}]",
        f"      max_distance: {max_distance}",
        "      force: true",
    ])


def format_contact_block(rec_chain, eff_chain, rec_resnum, eff_resnum, max_distance):
    """max_distance is emitted verbatim, so callers round it as needed."""
    return "\n".join([
        "  - contact:",
        f"      token1: [{rec_chain}, {rec_resnum}]",
        f"      token2: [{eff_chain}, {eff_resnum}]",
        f"      max_distance: {max_distance}",
        "      force: true",
    ])

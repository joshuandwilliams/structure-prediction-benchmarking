"""Shared Cα-geometry helpers for the Boltz constraint extractors.

Both extract_constraints_boltz1.py and extract_constraints_boltz2.py derive
restraints from the Cα coordinates of a two-chain reference complex. The PDB
parsing, the pairwise-distance scan, and the YAML block formatting are
identical between them — only *which* blocks each emits differs (Boltz-1 is
pocket-only with max_distance pinned to 6.0; Boltz-2 adds dense per-pair
contact blocks). This module holds the shared logic so the two CLI scripts
stay thin and the geometry can be unit-tested directly.

This file is a helper, not a CLI entrypoint — Nextflow never invokes it by
name. It is imported by its sibling scripts; because Python puts a script's
own directory on sys.path[0], `import _constraint_geometry` resolves whether
the scripts run from the repo or from inside the Singularity container.
"""

from __future__ import annotations

import math


def read_ca_by_chain(pdb_path, chains=None):
    """Parse Cα coordinates from a PDB file, grouped by chain.

    Returns ``{chain_id: {resnum: (x, y, z)}}``. Only ATOM records whose atom
    name is exactly ``CA`` are read. If ``chains`` is given (a set/iterable of
    chain IDs), atoms on other chains are skipped.
    """
    wanted = set(chains) if chains is not None else None
    coords: dict[str, dict[int, tuple[float, float, float]]] = {}
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                ch = line[21]
                if wanted is not None and ch not in wanted:
                    continue
                resnum = int(line[22:26].strip())
                xyz = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
                coords.setdefault(ch, {})[resnum] = xyz
    return coords


def distance(a, b):
    """Euclidean distance between two 3-tuples."""
    return math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b)))


def pocket_residues(rec_ca, eff_ca, cutoff):
    """Receptor residue numbers with any Cα within ``cutoff`` Å of an effector
    Cα. Returned sorted ascending.

    ``rec_ca`` / ``eff_ca`` are ``{resnum: (x, y, z)}`` maps as produced by
    :func:`read_ca_by_chain`.
    """
    selected = set()
    for resnum, r_xyz in rec_ca.items():
        for e_xyz in eff_ca.values():
            if distance(r_xyz, e_xyz) <= cutoff:
                selected.add(resnum)
                break
    return sorted(selected)


def contact_pairs(rec_ca, eff_ca, cutoff, max_pairs):
    """Closest inter-chain Cα residue pairs within ``cutoff`` Å.

    Returns a list of ``(rec_resnum, eff_resnum, dist)`` tuples sorted by
    distance ascending and truncated to ``max_pairs``.
    """
    pairs = []
    for r_rn, r_xyz in rec_ca.items():
        for e_rn, e_xyz in eff_ca.items():
            d = distance(r_xyz, e_xyz)
            if d <= cutoff:
                pairs.append((r_rn, e_rn, d))
    pairs.sort(key=lambda x: x[2])
    return pairs[:max_pairs]


def format_pocket_block(rec_chain, eff_chain, residues, max_distance):
    """Render a Boltz ``pocket`` constraint block (without the leading
    ``constraints:`` header). Lines are returned joined by newlines, with no
    trailing newline.
    """
    contacts_str = ", ".join(f"[{rec_chain}, {r}]" for r in residues)
    return "\n".join([
        "  - pocket:",
        f"      binder: {eff_chain}",
        f"      contacts: [{contacts_str}]",
        f"      max_distance: {max_distance}",
        "      force: true",
    ])


def format_contact_block(rec_chain, eff_chain, rec_resnum, eff_resnum, max_distance):
    """Render a single Boltz ``contact`` constraint block (without the leading
    ``constraints:`` header). ``max_distance`` is emitted verbatim — callers
    round it as needed.
    """
    return "\n".join([
        "  - contact:",
        f"      token1: [{rec_chain}, {rec_resnum}]",
        f"      token2: [{eff_chain}, {eff_resnum}]",
        f"      max_distance: {max_distance}",
        "      force: true",
    ])

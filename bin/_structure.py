"""Shared structure I/O and superposition helpers.

Every script that reads coordinates goes through here, so PDB and mmCIF are
handled identically everywhere. gemmi does the parsing and Biopython the
sequence alignment, rather than each script carrying its own reader.
"""

from __future__ import annotations

import numpy as np


def one_letter(resname):
    """Three-letter residue name to one letter, X if unknown.

    gemmi's table already maps the modified residues these predictors emit
    (MSE, SEP, TPO, PTR) onto their parents.
    """
    import gemmi
    info = gemmi.find_tabulated_residue(resname)
    code = info.one_letter_code.upper() if info else ""
    return code if code.isalpha() else "X"


def read_chains(path, chains=None):
    """Parse Cα coordinates and sequence per chain.

    Returns ``{chain_id: (coords (N,3) float array, sequence str, resnums list)}``
    for the first model.

    Selection is any amino-acid residue carrying a Cα, rather than gemmi's
    polymer view, which declines to classify very short chains and would
    return nothing for them. Waters, ligands and ions are dropped because they
    are not amino acids. Only the first altloc of a residue is kept, so a
    residue with alternate conformations is not counted twice.
    """
    import gemmi
    st = gemmi.read_structure(str(path))

    # A file gemmi cannot make sense of yields zero models rather than raising,
    # and indexing st[0] then throws IndexError. Callers already handle "no
    # chains found", so degrade to that instead of taking the whole run down.
    # gemmi 0.6.5 hits this on the minimal mmCIF the esm library writes, which
    # 0.7.5 reads without complaint.
    if len(st) == 0:
        return {}

    st.remove_waters()

    wanted = set(chains) if chains is not None else None
    out = {}
    for chain in st[0]:
        if wanted is not None and chain.name not in wanted:
            continue
        xyz, seq, nums, seen = [], [], [], set()
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid():
                continue
            key = (res.seqid.num, res.seqid.icode)
            if key in seen:
                continue
            atom = res.find_atom("CA", "*")
            if atom is None:
                continue
            seen.add(key)
            xyz.append([atom.pos.x, atom.pos.y, atom.pos.z])
            seq.append(one_letter(res.name))
            nums.append(res.seqid.num)
        if xyz:
            out[chain.name] = (np.asarray(xyz, dtype=float), "".join(seq), nums)
    return out


def read_chain(path, chain_id):
    """Cα coordinates and sequence for one chain. Raises if absent."""
    chains = read_chains(path, [chain_id])
    if chain_id not in chains:
        raise KeyError(f"chain {chain_id!r} not found in {path}")
    coords, seq, _ = chains[chain_id]
    return coords, seq


def read_ca_by_chain(path, chains=None):
    """``{chain_id: {resnum: (x, y, z)}}``, the form the constraint code wants."""
    return {cid: {n: tuple(xyz) for n, xyz in zip(nums, coords)}
            for cid, (coords, _, nums) in read_chains(path, chains).items()}


# ── Sequence-based residue pairing ────────────────────────────────────────

_ALIGNER = None


def aligner():
    """Global Needleman-Wunsch aligner, ChimeraX matchmaker defaults.

    BLOSUM62 with gap open -10 and extend -0.5, so residue pairing matches
    matchmaker's for homologous chains with missing termini or loops.
    """
    global _ALIGNER
    if _ALIGNER is None:
        from Bio.Align import PairwiseAligner, substitution_matrices
        a = PairwiseAligner()
        a.mode = "global"
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score = -10.0
        a.extend_gap_score = -0.5
        _ALIGNER = a
    return _ALIGNER


def matched_indices(seq_a, seq_b):
    """Index pairs of ungapped aligned positions, as two parallel lists.

    Subset coordinate arrays with these to get row-matched Cα pairs.
    """
    if not seq_a or not seq_b:
        return [], []
    try:
        alns = aligner().align(seq_a, seq_b)
        if len(alns) == 0:   # pragma: no cover
            return [], []
        blocks_a, blocks_b = alns[0].aligned
    except Exception:
        return [], []

    idx_a, idx_b = [], []
    for (a0, a1), (b0, b1) in zip(blocks_a, blocks_b):
        for k in range(min(a1 - a0, b1 - b0)):
            idx_a.append(int(a0 + k))
            idx_b.append(int(b0 + k))
    return idx_a, idx_b


# ── Superposition ─────────────────────────────────────────────────────────

def kabsch(P, Q):
    """Superpose P onto Q. Returns (rmsd, R, t, n) with ``P @ R.T + t ~ Q``.

    Wraps Biopython's SVDSuperimposer so the SVD and the reflection guard are
    not reimplemented here.
    """
    from Bio.SVDSuperimposer import SVDSuperimposer
    n = min(len(P), len(Q))
    if n == 0:
        return float("nan"), np.eye(3), np.zeros(3), 0
    P, Q = np.asarray(P[:n], dtype=float), np.asarray(Q[:n], dtype=float)

    sup = SVDSuperimposer()
    sup.set(Q, P)               # (reference, mobile)
    sup.run()
    R, t = sup.get_rotran()
    return float(sup.get_rms()), R, t, n


def apply_transform(coords, R, t):
    return np.asarray(coords, dtype=float) @ R + t


def kabsch_pruned(P, Q, cutoff=2.0, max_iter=10, min_pairs=20):
    """Fit iteratively with outlier pruning, as ChimeraX matchmaker does.

    Fits, drops pairs further apart than ``cutoff``, refits, and repeats until
    the surviving set stops changing or would fall below ``min_pairs``. Returns
    the same tuple as ``kabsch`` plus the surviving mask and iteration count.
    """
    P, Q = np.asarray(P, dtype=float), np.asarray(Q, dtype=float)
    n = min(len(P), len(Q))
    P, Q = P[:n], Q[:n]
    mask = np.ones(n, dtype=bool)

    best = (float("nan"), np.eye(3), np.zeros(3), 0)
    best_mask, iters = mask.copy(), 0

    for iters in range(1, max_iter + 1):
        if mask.sum() < min_pairs:
            break
        rmsd, R, t, k = kabsch(P[mask], Q[mask])
        best, best_mask = (rmsd, R, t, k), mask.copy()
        d = np.linalg.norm(apply_transform(P, R, t) - Q, axis=1)
        new_mask = d <= cutoff
        if new_mask.sum() < min_pairs or np.array_equal(new_mask, mask):
            break
        mask = new_mask
    return (*best, best_mask, iters)


# ── Contacts ──────────────────────────────────────────────────────────────

def contact_pairs(a_coords, b_coords, cutoff):
    """Index pairs (i, j, distance) closer than ``cutoff``, nearest first.

    Uses a KD-tree, so this stays linear-ish rather than scanning every pair.
    """
    from scipy.spatial import cKDTree
    a, b = np.asarray(a_coords, dtype=float), np.asarray(b_coords, dtype=float)
    if not len(a) or not len(b):
        return []
    pairs = cKDTree(a).query_ball_tree(cKDTree(b), r=cutoff)
    out = [(i, j, float(np.linalg.norm(a[i] - b[j])))
           for i, js in enumerate(pairs) for j in js]
    out.sort(key=lambda p: p[2])
    return out

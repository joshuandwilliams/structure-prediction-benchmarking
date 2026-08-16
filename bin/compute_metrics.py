#!/usr/bin/env python3
"""Compute confidence and structural-accuracy metrics from predictor output.

Parsers cover Boltz-1/2, Chai-1, AlphaFold2-Multimer, AlphaFold 3, ColabFold and
ESMFold2. Confidence metrics are pLDDT, pTM, ipTM, mean and interface PAE, ipSAE
and the AlphaFold ranking score. With a reference complex it also emits Cα RMSDs,
of which rmsd_effector_receptor_aligned is the one the benchmark scores on.

Usage:
    compute_metrics.py --model boltz2 --prediction-dir DIR \
        --chain-lengths 245 112 --output-csv metrics.csv
"""

import argparse
import csv
import glob
import json
import os
import shutil
import sys
import traceback

import _structure
import numpy as np

# ═══════════════════════════════════════════════════════════════════════════
# ipSAE (Dunbrack 2025)
# ═══════════════════════════════════════════════════════════════════════════

def calc_d0(n_res):
    """TM-score distance normaliser d0, as ipSAE defines it.

    This is the standard TM-score d0(L) = 1.24*(L-15)^(1/3) - 1.8, with the
    two clamps ipSAE applies: L is floored at 27 (below that the cube root
    term collapses and d0 goes negative), and the result is floored at 1.0 Å.

    Reference: Dunbrack (2025), "Res ipSAE loquuntur", doi:10/nq4x.
    """
    L = max(float(n_res), 27.0)
    return max(1.0, 1.24 * (L - 15.0) ** (1.0 / 3.0) - 1.8)


def compute_ipsae_one_direction(pae_ij, cutoff=10.0):
    """IpSAE for direction i->j. Returns (score, best_residue_idx).

    For each residue i of the first chain, take the residues j of the second
    chain whose interchain PAE is below `cutoff`. That subset defines both the
    residues scored and the length that sets d0, which is what makes ipSAE
    insensitive to the size of the non-interacting remainder — the property
    ipTM lacks. The directional score is the best single residue's mean.
    """
    n_i, n_j = pae_ij.shape
    scores = []
    for i in range(n_i):
        row = pae_ij[i, :]
        below = row[row < cutoff]
        n_below = len(below)
        if n_below == 0:
            scores.append(0.0)
            continue
        d0 = calc_d0(n_below)
        s = 1.0 / (1.0 + (below / d0) ** 2)
        scores.append(float(s.mean()))
    if not scores:
        return 0.0, -1
    return max(scores), int(np.argmax(scores))


def compute_ipsae(pae_matrix, chain_lengths, cutoff=10.0):
    """Compute ipSAE from full pAE matrix and chain length list."""
    if pae_matrix is None or len(chain_lengths) < 2:
        return {"ipsae_ab": 0.0, "ipsae_ba": 0.0, "ipsae_min": 0.0}

    len_a = chain_lengths[0]
    len_b = sum(chain_lengths[1:])
    total = pae_matrix.shape[0]

    if total != len_a + len_b:
        print(f"  WARNING: pAE matrix size {total} != expected {len_a + len_b}")
        if total > 0 and len_a < total:
            len_b = total - len_a
        else:
            return {"ipsae_ab": 0.0, "ipsae_ba": 0.0, "ipsae_min": 0.0}

    pae_ab = pae_matrix[:len_a, len_a:len_a + len_b]
    pae_ba = pae_matrix[len_a:len_a + len_b, :len_a]

    ipsae_ab, _ = compute_ipsae_one_direction(pae_ab, cutoff)
    ipsae_ba, _ = compute_ipsae_one_direction(pae_ba, cutoff)

    return {
        "ipsae_ab": round(ipsae_ab, 4),
        "ipsae_ba": round(ipsae_ba, 4),
        "ipsae_min": round(min(ipsae_ab, ipsae_ba), 4),
    }


def compute_ipae(pae_matrix, chain_lengths):
    """Mean interchain PAE (both off-diagonal blocks)."""
    if pae_matrix is None or len(chain_lengths) < 2:
        return 0.0
    len_a = chain_lengths[0]
    total = pae_matrix.shape[0]
    if total != sum(chain_lengths):
        return 0.0
    block_ab = pae_matrix[:len_a, len_a:]
    block_ba = pae_matrix[len_a:, :len_a]
    all_interchain = np.concatenate([block_ab.flatten(), block_ba.flatten()])
    return round(float(all_interchain.mean()), 2)


def _pae_derived_metrics(pae_matrix, chain_lengths):
    """Compute all PAE-derived metrics in one call."""
    pae_mean = round(float(pae_matrix.mean()), 2) if pae_matrix is not None else 0.0
    ipae = compute_ipae(pae_matrix, chain_lengths)
    ipsae = compute_ipsae(pae_matrix, chain_lengths)
    return {"pae_mean": pae_mean, "ipae": ipae, **ipsae}


# ═══════════════════════════════════════════════════════════════════════════
# pLDDT from PDB B-factor column
# ═══════════════════════════════════════════════════════════════════════════

def plddt_from_pdb(pdb_path):
    """Extract average pLDDT from B-factor column of CA atoms."""
    bfactors = []
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    try:
                        bfactors.append(float(line[60:66]))
                    except (ValueError, IndexError):
                        pass
    except Exception:
        return 0.0
    if not bfactors:
        return 0.0
    return round(float(np.mean(bfactors)), 2)


def plddt_from_cif(cif_path):
    """Average pLDDT from the B-factor column of Cα atoms in an mmCIF."""
    import gemmi
    try:
        st = gemmi.read_structure(str(cif_path))
        b = [res.find_atom("CA", "*").b_iso
             for chain in st[0] for res in chain
             if res.find_atom("CA", "*") is not None]
    except Exception:
        return 0.0
    return round(float(np.mean(b)), 2) if b else 0.0


def _finalize_entry(entry, pae_matrix, chain_lengths):
    """Apply defaults, compute PAE-derived metrics, compute ranking_score."""
    entry.setdefault("ptm", 0.0)
    entry.setdefault("iptm", 0.0)
    entry.setdefault("avg_plddt", 0.0)
    entry.update(_pae_derived_metrics(pae_matrix, chain_lengths))
    entry["ranking_score"] = round(0.8 * entry["iptm"] + 0.2 * entry["ptm"], 4)
    return entry


# ═══════════════════════════════════════════════════════════════════════════
# Structural RMSD against reference PDB
# ═══════════════════════════════════════════════════════════════════════════

def _read_ca_seq_by_chain(path):
    """{chain: (coords (N,3), sequence)} for PDB or mmCIF alike."""
    return {cid: (xyz, seq) for cid, (xyz, seq, _) in _structure.read_chains(path).items()}


_read_ca_seq_by_chain_cif = _read_ca_seq_by_chain


def _seqalign_pair_indices(seq_pred, seq_ref):
    """Ungapped aligned index pairs, ChimeraX matchmaker alignment parameters."""
    return _structure.matched_indices(seq_pred, seq_ref)


def _kabsch_align(P, Q):
    """Superpose P onto Q. Returns (rmsd, R, t, n) with P @ R.T + t ~ Q."""
    rmsd, R, t, n = _structure.kabsch(P, Q)
    return rmsd, R.T, t, n


def _kabsch_align_iterative(P, Q, cutoff=2.0, max_iter=10, min_pairs=20):
    """Matchmaker-style iterative fit with outlier pruning.

    Returns (rmsd, R, t, n_used, n_start, n_iter).
    """
    rmsd, R, t, n, mask, iters = _structure.kabsch_pruned(
        P, Q, cutoff=cutoff, max_iter=max_iter, min_pairs=min_pairs)
    return rmsd, R.T, t, int(mask.sum()), min(len(P), len(Q)), iters


def _resolve_pred_chains(pred_chains_seq, ref_chains_seq, rec_chain, eff_chain, pred_pdb):
    """Match predicted chains to the reference by sequence, not by chain ID.

    Predicted chain IDs are not trustworthy. Boltz writes A and B regardless of
    what was passed in, and AF2M wrote B and C on four benchmark targets. Each
    predicted chain is scored against both reference chains by matched residues
    over the longer length, then assigned greedily, receptor first.
    """
    # Sanity-check the reference chain IDs.
    if rec_chain not in ref_chains_seq or len(ref_chains_seq[rec_chain][0]) == 0:   # pragma: no cover
        raise ValueError(
            f"Reference chain '{rec_chain}' not found in reference PDB. "
            f"Available reference chains: {list(ref_chains_seq.keys())}"
        )
    if eff_chain not in ref_chains_seq or len(ref_chains_seq[eff_chain][0]) == 0:   # pragma: no cover
        raise ValueError(
            f"Reference chain '{eff_chain}' not found in reference PDB. "
            f"Available reference chains: {list(ref_chains_seq.keys())}"
        )

    ref_rec_seq = ref_chains_seq[rec_chain][1]
    ref_eff_seq = ref_chains_seq[eff_chain][1]

    pred_chain_ids = list(pred_chains_seq.keys())
    if len(pred_chain_ids) < 2:
        raise ValueError(
            f"Predicted file has fewer than 2 chains: {pred_chain_ids}\n"
            f"  File: {pred_pdb}"
        )

    # Score every predicted chain against both reference targets.
    # Score = (#aligned ungapped pairs with matching residue) / max(len_pred, len_ref).
    # This gives a fraction in [0, 1] that's 1.0 for identical chains.
    def identity(pred_seq, ref_seq):
        if not pred_seq or not ref_seq:   # pragma: no cover
            return 0.0
        pi, ri = _seqalign_pair_indices(pred_seq, ref_seq)
        if not pi:   # pragma: no cover
            return 0.0
        matches = sum(1 for p, r in zip(pi, ri) if pred_seq[p] == ref_seq[r])
        return matches / max(len(pred_seq), len(ref_seq))

    scores = {}  # cid -> (id_vs_receptor, id_vs_effector)
    for cid in pred_chain_ids:
        pred_seq = pred_chains_seq[cid][1]
        scores[cid] = (
            identity(pred_seq, ref_rec_seq),
            identity(pred_seq, ref_eff_seq),
        )

    # Greedy assignment: receptor first.
    best_rec_id = max(pred_chain_ids, key=lambda c: scores[c][0])
    best_rec_score = scores[best_rec_id][0]

    # Effector from remaining chains.
    remaining = [c for c in pred_chain_ids if c != best_rec_id]
    best_eff_id = max(remaining, key=lambda c: scores[c][1])
    best_eff_score = scores[best_eff_id][1]

    MIN_IDENTITY = 0.30  # below this, the chain is almost certainly the wrong protein
    if best_rec_score < MIN_IDENTITY:   # pragma: no cover
        raise ValueError(
            f"No predicted chain has >{MIN_IDENTITY*100:.0f}% identity to "
            f"reference receptor (chain '{rec_chain}', "
            f"length {len(ref_rec_seq)}).\n"
            f"  Predicted chain identities to receptor: "
            f"{ {c: round(scores[c][0], 2) for c in pred_chain_ids} }\n"
            f"  File: {pred_pdb}\n"
            f"  Either the wrong reference PDB was supplied or the chain "
            f"IDs in params.yml are wrong."
        )
    if best_eff_score < MIN_IDENTITY:   # pragma: no cover
        raise ValueError(
            f"No predicted chain has >{MIN_IDENTITY*100:.0f}% identity to "
            f"reference effector (chain '{eff_chain}', "
            f"length {len(ref_eff_seq)}).\n"
            f"  Predicted chain identities to effector: "
            f"{ {c: round(scores[c][1], 2) for c in pred_chain_ids} }\n"
            f"  File: {pred_pdb}"
        )

    print(f"    chain mapping (sequence-identity): "
          f"pred {best_rec_id} → receptor "
          f"({best_rec_score*100:.0f}% id to ref {rec_chain}), "
          f"pred {best_eff_id} → effector "
          f"({best_eff_score*100:.0f}% id to ref {eff_chain})")

    pred_rec_coords = pred_chains_seq[best_rec_id][0]
    pred_eff_coords = pred_chains_seq[best_eff_id][0]
    return best_rec_id, best_eff_id, pred_rec_coords, pred_eff_coords

def compute_structural_rmsds(pred_pdb, ref_pdb, rec_chain="A", eff_chain="B",
                             matchmaker_cutoff=2.0, min_pairs=20, max_iter=10):
    """
    Compute chain-specific RMSD metrics against a reference PDB.

    Pairing strategy
    ----------------
    Predicted and reference Cα atoms are paired by **global sequence
    alignment** (Needleman-Wunsch, BLOSUM62, gap-open -10, gap-extend
    -0.5 — ChimeraX matchmaker defaults).  Only positions that align
    without a gap on either side are kept.  This handles the common
    case where a prediction has residues at a terminus that the
    reference's crystal structure did not model: file-order pairing
    would then shift every Cα by one position from the start of the
    chain, producing a uniform large RMSD on what is actually a fine
    prediction.  Sequence alignment correctly skips the unmatched
    residue.

    After pairing, two RMSDs are computed for every metric:

    PRUNED  (matchmaker-equivalent, primary reporting columns)
        Iterative Kabsch with per-iteration pruning of pairs whose
        post-fit distance exceeds `matchmaker_cutoff` Å (default 2.0).
        Converges when the surviving set stops changing, up to
        `max_iter` iterations, refusing to prune below `min_pairs`
        survivors.  Distance-blind: it does not distinguish interface
        residues from mobile loops.  Inspect n_*_pruned to see how
        much was cut before trusting the number.

    RAW  (single-pass, no pruning, no iteration — _raw columns)
        Plain Kabsch on all sequence-aligned pairs.  Sensitive to
        single large outliers that drag the fit.  Kept for comparison
        so you can see how much pruning helped.

    Note that BOTH variants now use sequence-aligned pairing — the only
    difference is whether iterative pruning is applied after the fit.
    The previous file-order pairing has been removed because it was
    silently producing wrong numbers whenever termini differed between
    prediction and reference.

    Three metric families
    ---------------------
    Each family has a "whole-fit" primary variant and a "core" pruned
    variant.  The primary variant (unsuffixed name) uses all sequence-
    aligned pairs in a single-pass Kabsch fit.  The core variant
    (`_core` suffix, or `_corefit` for the docking metric) applies
    ChimeraX matchmaker-equivalent iterative outlier pruning —
    drops pairs > matchmaker_cutoff Å after the fit, refits, repeats
    — giving a tighter local RMSD on the well-fitting subset.

    1. rmsd_receptor  /  rmsd_receptor_core
       Kabsch of predicted receptor Cα onto reference receptor Cα.
       Whole-fit vs pruned-core.  Receptor atoms only.

    2. rmsd_effector_independent  /  rmsd_effector_independent_core
       Kabsch of predicted effector Cα onto reference effector Cα,
       fitted independently.  Measures intrinsic effector fold
       accuracy, ignoring placement.

    3. rmsd_effector_receptor_aligned  /  rmsd_effector_receptor_aligned_corefit
       Take the receptor Kabsch transform from metric 1 (whole-fit
       for the primary, core-pruned for _corefit) and apply it
       rigidly to the predicted effector Cα.  Compute RMSD against
       the reference effector across ALL aligned effector pairs —
       NO effector pruning, because pruning would hide misplacement
       of the very atoms we want to evaluate.

       The primary `rmsd_effector_receptor_aligned` is the number
       to rank designs by: it measures "once I line up the whole
       receptor, where is the effector?" — which matches what a
       structural biologist visually evaluates.  The `_corefit`
       variant uses the pruned receptor transform and is typically
       slightly worse on well-folded predictions, because the
       pruned transform is optimised against the rigid core, not
       the whole chain, and that small rotation difference
       propagates to the effector position.

    Returns
    -------
    dict with keys:
        pred_receptor_chain                -- str, chain ID used in predicted file
        pred_effector_chain                -- str, chain ID used in predicted file

        # Whole-fit (primary) — use these for ranking
        rmsd_receptor                      -- float or None
        rmsd_effector_independent          -- float or None
        rmsd_effector_receptor_aligned     -- float or None  (primary docking metric)
        n_receptor_ca                      -- int, sequence-aligned pairs used
        n_effector_ca                      -- int

        # Core (matchmaker-pruned) — tighter local fit on rigid core
        rmsd_receptor_core                      -- float or None
        rmsd_effector_independent_core          -- float or None
        rmsd_effector_receptor_aligned_corefit  -- float or None
        n_receptor_ca_core_start                -- int, aligned pairs before pruning
        n_receptor_ca_core_used                 -- int, pairs after pruning
        n_receptor_ca_core_pruned               -- int, pairs dropped by pruning
        n_receptor_core_iter                    -- int, refit iterations
        n_effector_ca_core_start                -- int
        n_effector_ca_core_used                 -- int
        n_effector_ca_core_pruned               -- int
        n_effector_core_iter                    -- int
    """
    nan_result = {
        "pred_receptor_chain":                   None,
        "pred_effector_chain":                   None,
        # Whole-fit (primary)
        "rmsd_receptor":                         None,
        "rmsd_effector_independent":             None,
        "rmsd_effector_receptor_aligned":        None,
        "n_receptor_ca":                         0,
        "n_effector_ca":                         0,
        # Core (matchmaker-pruned)
        "rmsd_receptor_core":                    None,
        "rmsd_effector_independent_core":        None,
        "rmsd_effector_receptor_aligned_corefit": None,
        "n_receptor_ca_core_start":              0,
        "n_receptor_ca_core_used":               0,
        "n_receptor_ca_core_pruned":             0,
        "n_receptor_core_iter":                  0,
        "n_effector_ca_core_start":              0,
        "n_effector_ca_core_used":               0,
        "n_effector_ca_core_pruned":             0,
        "n_effector_core_iter":                  0,
    }

    if not ref_pdb or not os.path.exists(ref_pdb):   # pragma: no cover
        print(f"  ERROR: Reference PDB not found: {ref_pdb}")
        return nan_result
    if not pred_pdb or not os.path.exists(pred_pdb):
        print(f"  ERROR: Predicted file not found: {pred_pdb}")
        return nan_result

    # ── Read predicted + reference (coords AND sequences) ──────────────
    if pred_pdb.endswith(".cif"):   # pragma: no cover
        pred_chains_seq = _read_ca_seq_by_chain_cif(pred_pdb)
    else:
        pred_chains_seq = _read_ca_seq_by_chain(pred_pdb)
    ref_chains_seq = _read_ca_seq_by_chain(ref_pdb)

    if not pred_chains_seq:
        print(f"  ERROR: No Cα atoms read from predicted file: {pred_pdb}")
        return nan_result
    if not ref_chains_seq:
        print(f"  ERROR: No Cα atoms read from reference PDB: {ref_pdb}")
        return nan_result

    # Sequence-similarity-based chain resolution.  This handles cases
    # where the prediction has chain IDs that don't match the reference
    # (Boltz writes A/B regardless of params), or where the reference
    # has multiple chains with the same sequence (homodimers in the
    # crystal asymmetric unit).
    try:
        pred_rec_id, pred_eff_id, _, _ = _resolve_pred_chains(
            pred_chains_seq, ref_chains_seq, rec_chain, eff_chain, pred_pdb
        )
    except ValueError as e:
        print(f"  ERROR: Chain resolution failed for {os.path.basename(pred_pdb)}:\n"
              f"    {e}\n"
              f"  Predicted chains available: "
              f"{ {c: len(tup[0]) for c, tup in pred_chains_seq.items()} }\n"
              f"  RMSD will not be computed for this structure.")
        return nan_result

    result = dict(nan_result)
    result["pred_receptor_chain"] = pred_rec_id
    result["pred_effector_chain"] = pred_eff_id

    # Pull coords + sequences for the resolved chains.
    pred_rec_coords, pred_rec_seq = pred_chains_seq[pred_rec_id]
    pred_eff_coords, pred_eff_seq = pred_chains_seq[pred_eff_id]
    ref_rec_coords,  ref_rec_seq  = ref_chains_seq [rec_chain]
    ref_eff_coords,  ref_eff_seq  = ref_chains_seq [eff_chain]

    if len(pred_rec_coords) == 0 or len(ref_rec_coords) == 0:   # pragma: no cover
        print(f"  ERROR: Receptor Cα count zero "
              f"(pred {pred_rec_id}: {len(pred_rec_coords)}, "
              f"ref {rec_chain}: {len(ref_rec_coords)}). RMSD not computed.")
        return result
    if len(pred_eff_coords) == 0 or len(ref_eff_coords) == 0:   # pragma: no cover
        print(f"  ERROR: Effector Cα count zero "
              f"(pred {pred_eff_id}: {len(pred_eff_coords)}, "
              f"ref {eff_chain}: {len(ref_eff_coords)}). RMSD not computed.")
        return result

    # ── Sequence-align receptor and effector independently ─────────────
    # Each metric pair (receptor / effector) gets its own alignment, so a
    # missing residue in the effector chain doesn't affect receptor pairing
    # and vice versa.
    rec_pred_idx, rec_ref_idx = _seqalign_pair_indices(pred_rec_seq, ref_rec_seq)
    eff_pred_idx, eff_ref_idx = _seqalign_pair_indices(pred_eff_seq, ref_eff_seq)

    if len(rec_pred_idx) == 0:   # pragma: no cover
        print(f"  ERROR: Receptor sequence alignment produced zero pairs. "
              f"pred_seq_len={len(pred_rec_seq)}, ref_seq_len={len(ref_rec_seq)}. "
              f"RMSD not computed.")
        return result
    if len(eff_pred_idx) == 0:   # pragma: no cover
        print(f"  ERROR: Effector sequence alignment produced zero pairs. "
              f"pred_seq_len={len(pred_eff_seq)}, ref_seq_len={len(ref_eff_seq)}. "
              f"RMSD not computed.")
        return result

    # Subset coords to the aligned positions.
    pred_rec = pred_rec_coords[rec_pred_idx]
    ref_rec  = ref_rec_coords [rec_ref_idx]
    pred_eff = pred_eff_coords[eff_pred_idx]
    ref_eff  = ref_eff_coords [eff_ref_idx]

    # ══════════════════════════════════════════════════════════════════
    # RAW (single-pass) metrics — unpruned Kabsch on seq-aligned pairs
    # ══════════════════════════════════════════════════════════════════
    rmsd_rec_raw, R_raw, t_raw, n_rec_raw = _kabsch_align(pred_rec, ref_rec)
    result["rmsd_receptor"] = round(rmsd_rec_raw, 2)
    result["n_receptor_ca"] = n_rec_raw

    rmsd_eff_ind_raw, _, _, n_eff_raw = _kabsch_align(pred_eff, ref_eff)
    result["rmsd_effector_independent"] = round(rmsd_eff_ind_raw, 2)
    result["n_effector_ca"] = n_eff_raw

    Pe_transformed_raw = (R_raw @ pred_eff.T).T + t_raw
    rmsd_eff_rec_aligned_raw = float(
        np.sqrt(np.mean(np.sum((Pe_transformed_raw - ref_eff) ** 2, axis=1)))
    )
    result["rmsd_effector_receptor_aligned"] = round(rmsd_eff_rec_aligned_raw, 2)

    # ══════════════════════════════════════════════════════════════════
    # PRUNED (matchmaker-equivalent) metrics — iterative Kabsch
    # ══════════════════════════════════════════════════════════════════
    (rmsd_rec_p, R_p, t_p,
     n_rec_used, n_rec_start, n_rec_iter) = _kabsch_align_iterative(
        pred_rec, ref_rec,
        cutoff=matchmaker_cutoff, max_iter=max_iter, min_pairs=min_pairs
    )
    result["rmsd_receptor_core"]         = None if rmsd_rec_p != rmsd_rec_p else round(rmsd_rec_p, 2)
    result["n_receptor_ca_core_start"]   = n_rec_start
    result["n_receptor_ca_core_used"]    = n_rec_used
    result["n_receptor_ca_core_pruned"]  = max(0, n_rec_start - n_rec_used)
    result["n_receptor_core_iter"]       = n_rec_iter

    (rmsd_eff_p, _, _,
     n_eff_used, n_eff_start, n_eff_iter) = _kabsch_align_iterative(
        pred_eff, ref_eff,
        cutoff=matchmaker_cutoff, max_iter=max_iter, min_pairs=min_pairs
    )
    result["rmsd_effector_independent_core"] = None if rmsd_eff_p != rmsd_eff_p else round(rmsd_eff_p, 2)
    result["n_effector_ca_core_start"]       = n_eff_start
    result["n_effector_ca_core_used"]        = n_eff_used
    result["n_effector_ca_core_pruned"]      = max(0, n_eff_start - n_eff_used)
    result["n_effector_core_iter"]           = n_eff_iter

    # Receptor-aligned effector RMSD, pruned variant: pruned receptor R,t
    # applied to ALL aligned effector pairs (no effector pruning).
    if rmsd_rec_p == rmsd_rec_p:  # not NaN
        Pe_transformed_p = (R_p @ pred_eff.T).T + t_p
        rmsd_eff_rec_aligned_p = float(
            np.sqrt(np.mean(np.sum((Pe_transformed_p - ref_eff) ** 2, axis=1)))
        )
        result["rmsd_effector_receptor_aligned_corefit"] = round(rmsd_eff_rec_aligned_p, 2)
    else:
        print(f"  WARNING: Pruned receptor fit failed for "
              f"{os.path.basename(pred_pdb)}; receptor-aligned effector "
              f"RMSD (pruned) left as None.")

    return result

def _dedup_paths(paths):
    """Deduplicate file paths preserving order."""
    seen = set()
    out = []
    for p in paths:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def _is_reference(path):
    """Check if file is our copied reference PDB."""
    return os.path.basename(path) == "reference.pdb"


# ═══════════════════════════════════════════════════════════════════════════
# Model-specific parsers
# ═══════════════════════════════════════════════════════════════════════════

def parse_boltz2(pred_dir, chain_lengths):
    """
    Parse Boltz-2 outputs.

    Output structure:
        output/predictions/<job>/pdb/<n>_model_<M>.pdb
        output/predictions/<job>/pae/<n>_model_<M>.npz
        output/predictions/<job>/confidence/<n>_model_<M>.json
        output/predictions/<job>/plddt/<n>_model_<M>.npz

    Or flatter layouts depending on version.
    """
    results = []

    pdb_files = sorted(glob.glob(os.path.join(pred_dir, "**", "*.pdb"), recursive=True))
    pdb_files = [p for p in pdb_files if not _is_reference(p)]

    if not pdb_files:
        print(f"  WARNING: No PDB files found in {pred_dir}")
        return results

    for pdb_path in pdb_files:
        name = os.path.basename(pdb_path).replace(".pdb", "")
        pdb_dir = os.path.dirname(pdb_path)
        entry = {"pdb_path": pdb_path, "model_name": name}

        # --- Confidence JSON ---
        conf_files = _dedup_paths(
            glob.glob(os.path.join(pdb_dir, f"confidence_{name}.json")) +
            glob.glob(os.path.join(pdb_dir, "..", "confidence", f"confidence_{name}.json")) +
            glob.glob(os.path.join(pred_dir, "**", f"confidence_{name}.json"), recursive=True)
        )
        if conf_files:
            try:
                with open(conf_files[0]) as f:
                    conf = json.load(f)
                entry["ptm"] = round(float(conf.get("ptm", 0.0)), 4)
                entry["iptm"] = round(float(
                    conf.get("iptm", conf.get("protein_iptm", 0.0))
                ), 4)
            except Exception as e:
                print(f"  WARNING: conf parse failed for {name}: {e}")

        # --- pLDDT ---
        plddt_files = _dedup_paths(
            glob.glob(os.path.join(pdb_dir, f"plddt_{name}.npz")) +
            glob.glob(os.path.join(pdb_dir, "..", "plddt", f"plddt_{name}.npz")) +
            glob.glob(os.path.join(pred_dir, "**", f"plddt_{name}.npz"), recursive=True)
        )
        if plddt_files:
            try:
                data = np.load(plddt_files[0])
                arr = data[list(data.keys())[0]]
                if arr.ndim > 1:
                    arr = arr.flatten()
                entry["avg_plddt"] = round(float(arr.mean()), 2)
            except Exception as e:
                print(f"  WARNING: plddt parse failed for {name}: {e}")
        if "avg_plddt" not in entry:
            entry["avg_plddt"] = plddt_from_pdb(pdb_path)

        # --- PAE matrix ---
        pae_files = _dedup_paths(
            glob.glob(os.path.join(pdb_dir, f"pae_{name}.npz")) +
            glob.glob(os.path.join(pdb_dir, "..", "pae", f"pae_{name}.npz")) +
            glob.glob(os.path.join(pred_dir, "**", f"pae_{name}.npz"), recursive=True)
        )
        pae_matrix = None
        if pae_files:
            try:
                data = np.load(pae_files[0])
                for key in data.keys():
                    arr = np.array(data[key])
                    if arr.ndim == 2:
                        pae_matrix = arr
                        break
                    elif arr.ndim == 3:
                        pae_matrix = arr[0]
                        break
            except Exception as e:
                print(f"  WARNING: pae parse failed for {name}: {e}")

        _finalize_entry(entry, pae_matrix, chain_lengths)
        results.append(entry)

    return results


def parse_boltz1(pred_dir, chain_lengths):
    """Parse Boltz-1 outputs. Same structure as Boltz-2."""
    return parse_boltz2(pred_dir, chain_lengths)


def parse_chai1(pred_dir, chain_lengths):
    """
    Parse Chai-1 v0.6.1 outputs.

    Output structure:
        output/pred.model_idx_0.cif
        output/scores.model_idx_0.npz
    """
    results = []

    struct_files = sorted(
        glob.glob(os.path.join(pred_dir, "**", "pred.model_idx_*.cif"), recursive=True) +
        glob.glob(os.path.join(pred_dir, "**", "pred.model_idx_*.pdb"), recursive=True)
    )
    if not struct_files:   # pragma: no cover
        struct_files = sorted(
            glob.glob(os.path.join(pred_dir, "**", "*.cif"), recursive=True) +
            glob.glob(os.path.join(pred_dir, "**", "*.pdb"), recursive=True)
        )
        struct_files = [p for p in struct_files if not _is_reference(p)]

    for struct_path in _dedup_paths(struct_files):
        name = os.path.basename(struct_path).rsplit(".", 1)[0]
        base = os.path.dirname(struct_path)
        entry = {"pdb_path": struct_path, "model_name": name}

        # Match scores file: pred.model_idx_0 -> scores.model_idx_0
        score_name = name.replace("pred.", "scores.")
        score_files = _dedup_paths(
            glob.glob(os.path.join(base, f"{score_name}.npz")) +
            glob.glob(os.path.join(base, f"scores.{name.replace('pred.', '')}.npz")) +
            glob.glob(os.path.join(base, "scores*.npz"))
        )

        pae_matrix = None
        for sf in score_files:
            try:
                data = np.load(sf)
                keys = list(data.keys())
                # Chai-1 v0.6.x NPZ contains only scalar confidence scores,
                # no pLDDT or PAE arrays.
                for metric in ("ptm", "iptm"):
                    if metric in keys:
                        v = np.array(data[metric])
                        entry[metric] = round(float(v.item() if v.ndim == 0 else v.flat[0]), 4)
                # PAE not present in v0.6.x — leave pae_matrix as None
                break
            except Exception as e:   # pragma: no cover
                print(f"  WARNING: Score parse failed for {sf}: {e}")

        # pLDDT is stored in the CIF B-factor column in Chai-1 v0.6.x
        if struct_path.endswith(".cif"):   # pragma: no cover
            entry["avg_plddt"] = plddt_from_cif(struct_path)
        elif struct_path.endswith(".pdb"):
            entry["avg_plddt"] = plddt_from_pdb(struct_path)

        _finalize_entry(entry, pae_matrix, chain_lengths)
        results.append(entry)

    return results


def parse_af2m(pred_dir, chain_lengths):
    """
    Parse AlphaFold2-Multimer outputs.

    Output structure:
        output/<fasta_name>/
            ranked_0.pdb ... ranked_4.pdb
            result_model_1_multimer_v3_pred_0.pkl
            ranking_debug.json

    pkl files contain: predicted_aligned_error (NxN), plddt (N,), ptm, iptm.
    ranking_debug.json: order list, iptm+ptm dict (COMBINED score, not separate).
    """
    import pickle
    results = []

    subdirs = [d for d in glob.glob(os.path.join(pred_dir, "*")) if os.path.isdir(d)]
    search_dir = subdirs[0] if subdirs else pred_dir

    pdb_files = sorted(glob.glob(os.path.join(search_dir, "ranked_*.pdb")))
    if not pdb_files:
        pdb_files = sorted(
            glob.glob(os.path.join(search_dir, "relaxed_model_*.pdb")) +
            glob.glob(os.path.join(search_dir, "unrelaxed_model_*.pdb"))
        )

    if not pdb_files:
        print(f"  WARNING: No PDB files found in {search_dir}")
        return results

    # Load ranking_debug.json
    ranking = {}
    ranking_file = os.path.join(search_dir, "ranking_debug.json")
    if os.path.exists(ranking_file):
        try:
            with open(ranking_file) as f:
                ranking = json.load(f)
        except Exception:
            pass

    # Load all pkl files
    pkl_metrics = {}
    for pkl_path in sorted(glob.glob(os.path.join(search_dir, "result_model_*.pkl"))):
        try:
            with open(pkl_path, "rb") as f:
                pkl_data = pickle.load(f)
            model_key = os.path.basename(pkl_path).replace("result_", "").replace(".pkl", "")
            m = {}
            if "predicted_aligned_error" in pkl_data:
                m["pae"] = np.array(pkl_data["predicted_aligned_error"])
            if "plddt" in pkl_data:
                m["plddt"] = np.array(pkl_data["plddt"])
            if "ptm" in pkl_data:
                m["ptm"] = float(pkl_data["ptm"])
            if "iptm" in pkl_data:
                m["iptm"] = float(pkl_data["iptm"])
            pkl_metrics[model_key] = m
        except Exception as e:
            print(f"  WARNING: pkl parse failed for {pkl_path}: {e}")

    order = ranking.get("order", [])

    for pdb_path in pdb_files:
        name = os.path.basename(pdb_path).replace(".pdb", "")
        entry = {"pdb_path": pdb_path, "model_name": name}

        pkl_data = None
        if name.startswith("ranked_"):
            try:
                rank_idx = int(name.replace("ranked_", ""))
                if rank_idx < len(order):
                    pkl_data = pkl_metrics.get(order[rank_idx])
            except ValueError:   # pragma: no cover
                pass
        else:   # pragma: no cover
            model_key = name.replace("relaxed_", "").replace("unrelaxed_", "")
            pkl_data = pkl_metrics.get(model_key)

        pae_matrix = None
        if pkl_data:
            if "ptm" in pkl_data:
                entry["ptm"] = round(pkl_data["ptm"], 4)
            if "iptm" in pkl_data:
                entry["iptm"] = round(pkl_data["iptm"], 4)
            if "plddt" in pkl_data:
                entry["avg_plddt"] = round(float(pkl_data["plddt"].mean()), 2)
            if "pae" in pkl_data:
                pae_matrix = pkl_data["pae"]

        if "avg_plddt" not in entry:
            entry["avg_plddt"] = plddt_from_pdb(pdb_path)

        _finalize_entry(entry, pae_matrix, chain_lengths)
        results.append(entry)

    return results


def _parse_af3_confidences(conf_path, chain_lengths):
    """
    Parse an AF3 confidences.json (atom-level pLDDT + full PAE matrix).

    Returns (avg_plddt, pae_matrix).  pae_matrix may be None.

    AF3 confidences.json fields:
        atom_plddts:  list of per-atom pLDDT values (0-100)
        pae:          NxN nested list of token-level PAE (Angstroms)
                      N = total residues across all chains (not atoms)

    The PAE matrix may be larger than sum(chain_lengths) when the PDB input
    had more chains than the two we are benchmarking (e.g. a 3-chain crystal
    structure used as the reference).  We trim to the expected size rather
    than failing.
    """
    avg_plddt = 0.0
    pae_matrix = None
    try:
        with open(conf_path) as f:
            cdata = json.load(f)

        if "atom_plddts" in cdata:
            avg_plddt = round(float(np.mean(cdata["atom_plddts"])), 2)
        elif "plddt" in cdata:
            val = cdata["plddt"]
            avg_plddt = round(float(np.mean(val) if isinstance(val, list) else val), 2)

        if "pae" in cdata:
            pae_raw = cdata["pae"]
            if isinstance(pae_raw, list) and pae_raw:
                if isinstance(pae_raw[0], list):
                    pae_matrix = np.array(pae_raw, dtype=np.float32)
                else:   # pragma: no cover
                    pae_flat = np.array(pae_raw, dtype=np.float32)
                    n = int(np.sqrt(len(pae_flat)))
                    if n * n == len(pae_flat):   # pragma: no cover
                        pae_matrix = pae_flat.reshape(n, n)

        # Trim PAE to expected size (handles 3-chain reference → 2-chain prediction)
        if pae_matrix is not None and chain_lengths:
            expected = sum(chain_lengths)
            if pae_matrix.shape[0] > expected:   # pragma: no cover
                pae_matrix = pae_matrix[:expected, :expected]
            elif pae_matrix.shape[0] < expected:   # pragma: no cover
                pae_matrix = None  # too small — discard rather than silently wrong

    except Exception as e:
        print(f"  WARNING: AF3 confidences parse failed for {conf_path}: {e}")
        traceback.print_exc()

    return avg_plddt, pae_matrix


def parse_af3(pred_dir, chain_lengths):
    """
    Parse AlphaFold 3 outputs.

    AF3 (run_alphafold.py, Dec 2024 NBI build) produces:

      output/<job_name>/                              ← top-level job dir
          <job_name>_model.cif                        ← aggregate best model
          <job_name>_confidences.json                 ← aggregate confidences
          <job_name>_summary_confidences.json         ← aggregate ptm/iptm
          ranking_scores.csv
          seed-<S>_sample-<N>/                        ← per-seed/sample dirs
              model.cif
              confidences.json
              summary_confidences.json

    We collect ALL per-seed/sample structures (the full 25-entry ensemble),
    plus the top-level aggregate model if present.  The aggregate is ranked
    by ranking_score alongside the per-seed entries.

    pred_dir is expected to be the 'output/' directory (one level above the
    job subdir).  We walk two levels deep to handle both layouts.
    """
    results = []

    # ── Find job-level subdirs (e.g. output/benchmark_test/) ──────────────
    job_dirs = [d for d in glob.glob(os.path.join(pred_dir, "*")) if os.path.isdir(d)]
    search_dirs = job_dirs if job_dirs else [pred_dir]

    for job_dir in search_dirs:
        # ── 1. Per-seed/sample subdirs: seed-N_sample-M/ ──────────────────
        seed_dirs = sorted([
            d for d in glob.glob(os.path.join(job_dir, "seed-*_sample-*"))
            if os.path.isdir(d)
        ])

        for seed_dir in seed_dirs:
            cif_path = os.path.join(seed_dir, "model.cif")
            if not os.path.exists(cif_path):   # pragma: no cover
                continue

            seed_sample = os.path.basename(seed_dir)  # e.g. "seed-42_sample-3"
            entry = {"pdb_path": cif_path, "model_name": seed_sample}

            # confidences.json lives alongside model.cif in the same seed dir
            conf_path = os.path.join(seed_dir, "confidences.json")
            avg_plddt, pae_matrix = (0.0, None)
            if os.path.exists(conf_path):
                avg_plddt, pae_matrix = _parse_af3_confidences(conf_path, chain_lengths)
            else:
                avg_plddt = plddt_from_cif(cif_path)

            entry["avg_plddt"] = avg_plddt

            # summary_confidences.json for ptm/iptm
            summary_path = os.path.join(seed_dir, "summary_confidences.json")
            if os.path.exists(summary_path):
                try:
                    with open(summary_path) as f:
                        sdata = json.load(f)
                    if "ptm" in sdata:
                        entry["ptm"] = round(float(sdata["ptm"]), 4)
                    if "iptm" in sdata:
                        entry["iptm"] = round(float(sdata["iptm"]), 4)
                except Exception as e:   # pragma: no cover
                    print(f"  WARNING: AF3 summary conf parse failed for {summary_path}: {e}")

            _finalize_entry(entry, pae_matrix, chain_lengths)
            results.append(entry)

        # ── 2. Top-level aggregate model (optional, already ranked by AF3) ─
        # Named <job_name>_model.cif at the job_dir level.
        # We include it so it participates in ranking but mark it clearly.
        agg_cifs = _dedup_paths(
            glob.glob(os.path.join(job_dir, "*_model.cif"))
        )
        for agg_cif in agg_cifs:
            name = os.path.basename(agg_cif).rsplit(".", 1)[0]  # e.g. benchmark_test_model
            entry = {"pdb_path": agg_cif, "model_name": f"{name}__aggregate"}

            conf_base = name.replace("_model", "")
            conf_path = os.path.join(job_dir, f"{conf_base}_confidences.json")
            avg_plddt, pae_matrix = (0.0, None)
            if os.path.exists(conf_path):
                avg_plddt, pae_matrix = _parse_af3_confidences(conf_path, chain_lengths)
            else:   # pragma: no cover
                avg_plddt = plddt_from_cif(agg_cif)

            entry["avg_plddt"] = avg_plddt

            summary_path = os.path.join(job_dir, f"{conf_base}_summary_confidences.json")
            if os.path.exists(summary_path):
                try:
                    with open(summary_path) as f:
                        sdata = json.load(f)
                    if "ptm" in sdata:
                        entry["ptm"] = round(float(sdata["ptm"]), 4)
                    if "iptm" in sdata:
                        entry["iptm"] = round(float(sdata["iptm"]), 4)
                except Exception as e:   # pragma: no cover
                    print(f"  WARNING: AF3 aggregate summary conf parse failed: {e}")

            _finalize_entry(entry, pae_matrix, chain_lengths)
            results.append(entry)

        # ── 3. Fallback: flat layout (old AF3 format, no seed subdirs) ─────
        if not seed_dirs:
            struct_files = _dedup_paths(
                glob.glob(os.path.join(job_dir, "*_model.cif")) +
                glob.glob(os.path.join(job_dir, "*.cif"))
            )
            if not struct_files:
                struct_files = [p for p in sorted(glob.glob(os.path.join(job_dir, "*.pdb")))
                                if not _is_reference(p)]
            for struct_path in struct_files:
                name = os.path.basename(struct_path).rsplit(".", 1)[0]
                base = os.path.dirname(struct_path)
                entry = {"pdb_path": struct_path, "model_name": name}
                conf_base = name.replace("_model", "")
                conf_path = os.path.join(base, f"{conf_base}_confidences.json")
                avg_plddt, pae_matrix = (0.0, None)
                if os.path.exists(conf_path):
                    avg_plddt, pae_matrix = _parse_af3_confidences(conf_path, chain_lengths)
                elif struct_path.endswith(".cif"):   # pragma: no cover
                    avg_plddt = plddt_from_cif(struct_path)
                else:   # pragma: no cover
                    avg_plddt = plddt_from_pdb(struct_path)
                entry["avg_plddt"] = avg_plddt
                summary_path = os.path.join(base, f"{conf_base}_summary_confidences.json")
                if os.path.exists(summary_path):
                    try:
                        with open(summary_path) as f:
                            sdata = json.load(f)
                        if "ptm" in sdata:
                            entry["ptm"] = round(float(sdata["ptm"]), 4)
                        if "iptm" in sdata:
                            entry["iptm"] = round(float(sdata["iptm"]), 4)
                    except Exception as e:   # pragma: no cover
                        print(f"  WARNING: AF3 flat summary conf parse failed: {e}")
                _finalize_entry(entry, pae_matrix, chain_lengths)
                results.append(entry)

    if not results:
        print(f"  WARNING: No AF3 predictions found in {pred_dir}")

    return results


def parse_colabfold(pred_dir, chain_lengths):
    """
    Parse ColabFold outputs.

    Output structure:
        output/<name>_unrelaxed_rank_001_..._model_1_seed_000.pdb
        output/<name>_scores_rank_001_..._model_1_seed_000.json

    Scores JSON: plddt (per-residue list 0-100), ptm, iptm, pae (NxN nested list)
    """
    results = []

    pdb_files = sorted(glob.glob(os.path.join(pred_dir, "**", "*.pdb"), recursive=True))
    pdb_files = [p for p in pdb_files if not _is_reference(p)]

    for pdb_path in pdb_files:
        name = os.path.basename(pdb_path).replace(".pdb", "")
        base = os.path.dirname(pdb_path)
        entry = {"pdb_path": pdb_path, "model_name": name}

        # ColabFold: _unrelaxed_ or _relaxed_ in PDB name → _scores_ in JSON
        scores_name = name.replace("_unrelaxed_", "_scores_").replace("_relaxed_", "_scores_")
        score_files = _dedup_paths(
            glob.glob(os.path.join(base, f"{scores_name}.json")) +
            glob.glob(os.path.join(base, name.replace("unrelaxed_", "scores_") + ".json"))
        )
        if not score_files:
            score_files = sorted(glob.glob(os.path.join(base, "*scores*.json")))

        pae_matrix = None
        for jf in _dedup_paths(score_files):
            try:
                with open(jf) as f:
                    sdata = json.load(f)

                if "ptm" in sdata:
                    entry["ptm"] = round(float(sdata["ptm"]), 4)
                if "iptm" in sdata:
                    entry["iptm"] = round(float(sdata["iptm"]), 4)
                if "plddt" in sdata and isinstance(sdata["plddt"], list):
                    entry["avg_plddt"] = round(float(np.mean(sdata["plddt"])), 2)
                if "pae" in sdata:
                    pae_matrix = np.array(sdata["pae"], dtype=np.float32)

                break
            except Exception as e:
                print(f"  WARNING: ColabFold score parse for {jf}: {e}")

        if "avg_plddt" not in entry:
            entry["avg_plddt"] = plddt_from_pdb(pdb_path)

        _finalize_entry(entry, pae_matrix, chain_lengths)
        results.append(entry)

    return results


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def parse_esmfold2(pred_dir, chain_lengths):
    """
    Parse ESMFold2 outputs.

    modules/esmfold2.nf writes, per seed, into the aggregated staging tree:

        all_outputs/<seed_tag>/esmfold2_pred.pdb      predicted complex (PDB)
        all_outputs/<seed_tag>/esmfold2_pred.cif      the same complex as mmCIF
        all_outputs/<seed_tag>/confidences.json       {plddt:[0-100], ptm, iptm?, pae?}

    The PDB is preferred where present. The mmCIF the esm library writes omits
    columns gemmi 0.6.5 requires, and that version returns zero models for it
    rather than raising, so the metrics container cannot read the mmCIF at all.
    Older runs predate the PDB, so the mmCIF remains the fallback.

    pLDDT is already rescaled to 0-100 by bin/esmfold2_fold.py. ipTM / PAE are
    present only if the model exposes them; when absent, the interface metrics
    fall back to 0 exactly as for any other model without a PAE matrix. The
    confidences.json layout mirrors AF3's, so this parser is a close cousin of
    parse_af3.
    """
    results = []

    cif_files = sorted(glob.glob(os.path.join(pred_dir, "**", "*.cif"), recursive=True))
    cif_files = [p for p in cif_files if not _is_reference(p)]

    # Swap in the sibling PDB wherever one exists.
    cif_files = [
        (p[:-4] + ".pdb") if os.path.isfile(p[:-4] + ".pdb") else p
        for p in cif_files
    ]

    for cif_path in cif_files:
        base = os.path.dirname(cif_path)
        name = os.path.relpath(cif_path, pred_dir).replace(os.sep, "_").rsplit(".", 1)[0]
        entry = {"pdb_path": cif_path, "model_name": name}

        pae_matrix = None
        conf_path = os.path.join(base, "confidences.json")
        if os.path.exists(conf_path):
            try:
                with open(conf_path) as f:
                    cdata = json.load(f)
                if "plddt" in cdata:
                    val = cdata["plddt"]
                    entry["avg_plddt"] = round(
                        float(np.mean(val) if isinstance(val, list) else val), 2)
                if "ptm" in cdata:
                    entry["ptm"] = round(float(cdata["ptm"]), 4)
                if "iptm" in cdata:
                    entry["iptm"] = round(float(cdata["iptm"]), 4)
                pae_raw = cdata.get("pae")
                if isinstance(pae_raw, list) and pae_raw and isinstance(pae_raw[0], list):
                    pae_matrix = np.array(pae_raw, dtype=np.float32)
                    if chain_lengths:
                        expected = sum(chain_lengths)
                        if pae_matrix.shape[0] > expected:
                            pae_matrix = pae_matrix[:expected, :expected]
                        elif pae_matrix.shape[0] < expected:
                            pae_matrix = None  # too small — discard rather than guess
            except Exception as e:
                print(f"  WARNING: ESMFold2 confidences parse failed for {conf_path}: {e}")

        if "avg_plddt" not in entry:
            entry["avg_plddt"] = plddt_from_cif(cif_path)

        _finalize_entry(entry, pae_matrix, chain_lengths)
        results.append(entry)

    if not results:
        print(f"  WARNING: No ESMFold2 predictions found in {pred_dir}")

    return results


PARSERS = {
    "boltz2": parse_boltz2,
    "boltz2_constrained": parse_boltz2,  # same output format; pocket + contact constraints
    "boltz1": parse_boltz1,
    "boltz1_constrained": parse_boltz1,  # same output format; pocket constraint only
    "chai1": parse_chai1,
    "af2m": parse_af2m,
    "colabfold_nomsa": parse_colabfold,  # same output format, single-sequence MSA input
    "af3": parse_af3,
    "af3_nomsa": parse_af3,              # same output format, no MSA/template search
    "colabfold": parse_colabfold,
    "esmfold2": parse_esmfold2,          # single-sequence diffusion complex predictor
}

CSV_FIELDS = [
    "model", "model_name", "pdb_path",
    "avg_plddt", "ptm", "iptm", "pae_mean", "ipae",
    "ipsae_ab", "ipsae_ba", "ipsae_min", "ranking_score",
    # Structural RMSD vs reference (None when --reference-pdb not provided).
    # pred_receptor_chain / pred_effector_chain record which chain IDs were
    # actually used in the predicted file — models relabel chains internally
    # so these may differ from the --receptor-chain / --effector-chain args.
    "pred_receptor_chain",
    "pred_effector_chain",
    # ── Whole-fit metrics (primary) ────────────────────────────────────
    # Sequence-aligned single-pass Kabsch across all paired residues.
    # These are the default reporting values and what the ranking CSV
    # sorts on.  rmsd_effector_receptor_aligned is the primary
    # docking-quality metric.
    "rmsd_receptor",
    "rmsd_effector_independent",
    "rmsd_effector_receptor_aligned",
    "n_receptor_ca",
    "n_effector_ca",
    # ── Core (matchmaker-pruned) metrics ───────────────────────────────
    # ChimeraX-matchmaker iterative Kabsch: fit, drop pairs > cutoff Å,
    # refit, repeat until stable.  Tighter local fit on a rigid core
    # subset.  The _corefit effector RMSD uses the pruned-receptor
    # transform; on well-folded predictions it is typically slightly
    # worse than the whole-fit version because the pruned transform
    # is optimised against the core, not the full chain.
    "rmsd_receptor_core",
    "rmsd_effector_independent_core",
    "rmsd_effector_receptor_aligned_corefit",
    "n_receptor_ca_core_start",
    "n_receptor_ca_core_used",
    "n_receptor_ca_core_pruned",
    "n_receptor_core_iter",
    "n_effector_ca_core_start",
    "n_effector_ca_core_used",
    "n_effector_ca_core_pruned",
    "n_effector_core_iter",
]


def main():
    parser = argparse.ArgumentParser(
        description="Compute confidence metrics from structure prediction outputs",
    )
    parser.add_argument("--model", required=True, choices=list(PARSERS.keys()))
    parser.add_argument("--prediction-dir", required=True)
    parser.add_argument("--chain-lengths", nargs="+", type=int, required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--best-model-dir", default=None,
                        help="Copy the best model (by ranking_score) here")
    parser.add_argument("--pae-cutoff", type=float, default=10.0)
    parser.add_argument(
        "--reference-pdb", default=None,
        help="Reference (crystal) PDB for structural RMSD calculation. "
             "When provided, computes rmsd_receptor, rmsd_effector_independent, "
             "and rmsd_effector_receptor_aligned (whole-fit, primary) plus the "
             "matchmaker-pruned _core/_corefit variants for every predicted model.",
    )
    parser.add_argument(
        "--receptor-chain", default="A",
        help="Chain ID of the receptor in predicted and reference PDBs (default: A)",
    )
    parser.add_argument(
        "--effector-chain", default="B",
        help="Chain ID of the effector in predicted and reference PDBs (default: B)",
    )
    parser.add_argument(
        "--matchmaker-cutoff", type=float, default=2.0,
        help="Distance cutoff (Angstroms) for iterative pruning in the "
             "matchmaker-equivalent RMSD calculation.  Pairs whose post-fit "
             "distance exceeds this value are dropped and the fit is refitted. "
             "ChimeraX's matchmaker default is 2.0 Å.  Lower = stricter.",
    )
    parser.add_argument(
        "--matchmaker-min-pairs", type=int, default=20,
        help="Minimum number of Cα pairs to keep in the iterative pruning. "
             "Pruning stops if it would drop survivors below this count, "
             "preventing meaningless fits on tiny residue sets.  Default: 20.",
    )
    parser.add_argument(
        "--matchmaker-max-iter", type=int, default=10,
        help="Maximum number of prune-and-refit iterations for the "
             "matchmaker-equivalent RMSD.  Default: 10 (usually converges "
             "in 2-4).",
    )
    args = parser.parse_args()

    print(f"=== Computing metrics for {args.model} ===")
    print(f"Prediction dir: {args.prediction_dir}")
    print(f"Chain lengths:  {args.chain_lengths}")

    if not os.path.exists(args.prediction_dir):   # pragma: no cover
        print(f"ERROR: Prediction directory not found: {args.prediction_dir}")
        sys.exit(1)

    parse_fn = PARSERS[args.model]

    try:
        results = parse_fn(args.prediction_dir, args.chain_lengths)
    except Exception as e:   # pragma: no cover
        print(f"ERROR: Parser failed for {args.model}: {e}")
        traceback.print_exc()
        results = []

    if not results:
        print(f"WARNING: No predictions found for {args.model}")
        os.makedirs(os.path.dirname(args.output_csv) or ".", exist_ok=True)
        with open(args.output_csv, "w") as f:
            f.write(",".join(CSV_FIELDS) + "\n")
        return

    for r in results:
        r["model"] = args.model

    # ── Structural RMSD vs reference ──────────────────────────────────────
    if args.reference_pdb:
        if not os.path.exists(args.reference_pdb):
            print(f"WARNING: --reference-pdb not found: {args.reference_pdb}")
        else:
            print(f"\nComputing structural RMSDs vs {args.reference_pdb}")
            print(f"  Reference receptor chain: {args.receptor_chain}, "
                  f"reference effector chain: {args.effector_chain}")
            print("  Method: whole-fit single-pass Kabsch (primary) "
                  "+ matchmaker-pruned core Kabsch")
            for r in results:
                pred_pdb = r.get("pdb_path", "")
                rmsds = compute_structural_rmsds(
                    pred_pdb, args.reference_pdb,
                    rec_chain=args.receptor_chain,
                    eff_chain=args.effector_chain,
                    matchmaker_cutoff=args.matchmaker_cutoff,
                    min_pairs=args.matchmaker_min_pairs,
                    max_iter=args.matchmaker_max_iter,
                )
                r.update(rmsds)
                print(f"  {os.path.basename(pred_pdb)}: "
                      f"pred_rec={rmsds['pred_receptor_chain']} "
                      f"pred_eff={rmsds['pred_effector_chain']}")
                print(f"    whole:  "
                      f"rec={rmsds['rmsd_receptor']} A  "
                      f"eff_ind={rmsds['rmsd_effector_independent']} A  "
                      f"eff_rec_aligned={rmsds['rmsd_effector_receptor_aligned']} A  "
                      f"(n_rec={rmsds['n_receptor_ca']}, "
                      f"n_eff={rmsds['n_effector_ca']})")
                print(f"    core:   "
                      f"rec={rmsds['rmsd_receptor_core']} A  "
                      f"eff_ind={rmsds['rmsd_effector_independent_core']} A  "
                      f"eff_rec_aligned={rmsds['rmsd_effector_receptor_aligned_corefit']} A")
                print(f"            rec pairs: "
                      f"{rmsds['n_receptor_ca_core_used']}/{rmsds['n_receptor_ca_core_start']} "
                      f"used ({rmsds['n_receptor_ca_core_pruned']} pruned, "
                      f"{rmsds['n_receptor_core_iter']} iter)  |  "
                      f"eff pairs: "
                      f"{rmsds['n_effector_ca_core_used']}/{rmsds['n_effector_ca_core_start']} "
                      f"used ({rmsds['n_effector_ca_core_pruned']} pruned, "
                      f"{rmsds['n_effector_core_iter']} iter)")
    else:
        for r in results:
            r.setdefault("pred_receptor_chain", None)
            r.setdefault("pred_effector_chain", None)
            r.setdefault("rmsd_receptor_core", None)
            r.setdefault("rmsd_effector_independent_core", None)
            r.setdefault("rmsd_effector_receptor_aligned_corefit", None)
            r.setdefault("n_receptor_ca_core_start", 0)
            r.setdefault("n_receptor_ca_core_used", 0)
            r.setdefault("n_receptor_ca_core_pruned", 0)
            r.setdefault("n_receptor_core_iter", 0)
            r.setdefault("n_effector_ca_core_start", 0)
            r.setdefault("n_effector_ca_core_used", 0)
            r.setdefault("n_effector_ca_core_pruned", 0)
            r.setdefault("n_effector_core_iter", 0)
            r.setdefault("rmsd_receptor", None)
            r.setdefault("rmsd_effector_independent", None)
            r.setdefault("rmsd_effector_receptor_aligned", None)
            r.setdefault("n_receptor_ca", 0)
            r.setdefault("n_effector_ca", 0)

    results.sort(key=lambda x: x.get("ranking_score", 0.0), reverse=True)

    os.makedirs(os.path.dirname(args.output_csv) or ".", exist_ok=True)
    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print(f"\nWrote {len(results)} entries to {args.output_csv}")

    # Which prediction to publish as this model's representative.
    #
    # The choice must be reference-free: selecting the lowest-RMSD prediction
    # would need the answer the benchmark is trying to measure, and could not be
    # reproduced on a novel target where no reference exists.  results is sorted
    # by AF ranking score above, so results[0] was previously that pick.  On the Tier-1
    # set, average pLDDT recovers more of the achievable ceiling than any
    # interface score: it poses 49.5% of targets correctly against 47.2% for
    # AF ranking score and 46.8% for ipTM, where selecting on RMSD would give 56.5%.
    # Rank-averaged combinations of the confidence scores did no better than
    # pLDDT alone.
    #
    # pLDDT is only compared within one model's own predictions here, so the
    # fact that Boltz reports 0-1 while the AlphaFold lineage reports 0-100 does
    # not affect the comparison.
    def _plddt_key(entry):
        v = entry.get("avg_plddt")
        return v if isinstance(v, (int, float)) and v > 0 else None

    scored = [r for r in results if _plddt_key(r) is not None]
    if scored:
        best = max(scored, key=_plddt_key)
        basis = "highest avg_plddt"
    else:
        best = results[0]
        basis = "highest ranking_score (no usable pLDDT)"

    print(f"\nSelected model ({args.model}): "
          f"{best.get('model_name', 'unknown')} [{basis}]")
    confidence_fields = [f for f in CSV_FIELDS[3:]
                         if not f.startswith("rmsd_") and not f.startswith("n_")]
    rmsd_fields = [f for f in CSV_FIELDS if f.startswith("rmsd_") or f.startswith("n_")]
    for k in confidence_fields:
        print(f"  {k:12s} = {best.get(k, 'N/A')}")
    if args.reference_pdb:
        print("  --- Structural RMSD vs reference ---")
        for k in rmsd_fields:
            v = best.get(k)
            print(f"  {k:35s} = {v if v is not None else 'N/A'}")

    if args.best_model_dir:
        os.makedirs(args.best_model_dir, exist_ok=True)
        src = best["pdb_path"]
        ext = os.path.splitext(src)[1]
        dst = os.path.join(args.best_model_dir, f"{args.model}_best{ext}")
        if os.path.exists(src):
            shutil.copy2(src, dst)
            print(f"Copied best model to {dst}")
        else:   # pragma: no cover
            print(f"WARNING: Best model file not found: {src}")


if __name__ == "__main__":
    main()

"""Unit tests for the PAE-derived confidence metrics in bin/compute_metrics.py
(ipSAE, ipAE). Pure numpy — runs on a laptop.
"""

import compute_metrics as cm
import numpy as np
import pytest

pytestmark = pytest.mark.local_unit


# ── compute_ipsae_one_direction ───────────────────────────────────────────

def test_one_direction_perfect_pae_scores_one():
    # All PAE = 0 (< cutoff): every term 1/(1+0) = 1 → mean 1 per row → 1.0.
    pae = np.zeros((4, 5))
    score, idx = cm.compute_ipsae_one_direction(pae, cutoff=10.0)
    assert score == pytest.approx(1.0)
    assert idx == 0


def test_one_direction_all_above_cutoff_scores_zero():
    # No entries below cutoff anywhere → every row contributes 0.
    pae = np.full((3, 3), 20.0)
    score, _ = cm.compute_ipsae_one_direction(pae, cutoff=10.0)
    assert score == 0.0


def test_one_direction_picks_best_row():
    # Row 1 is all-zero (score 1.0); other rows all above cutoff (score 0).
    pae = np.full((3, 4), 20.0)
    pae[1, :] = 0.0
    score, idx = cm.compute_ipsae_one_direction(pae, cutoff=10.0)
    assert score == pytest.approx(1.0)
    assert idx == 1


# ── compute_ipsae ─────────────────────────────────────────────────────────

def test_ipsae_perfect_interface():
    pae = np.zeros((4, 4))            # len_a=2, len_b=2
    out = cm.compute_ipsae(pae, [2, 2])
    assert out == {"ipsae_ab": 1.0, "ipsae_ba": 1.0, "ipsae_min": 1.0}


def test_ipsae_none_matrix_returns_zeros():
    assert cm.compute_ipsae(None, [2, 2]) == {
        "ipsae_ab": 0.0, "ipsae_ba": 0.0, "ipsae_min": 0.0
    }


def test_ipsae_single_chain_returns_zeros():
    assert cm.compute_ipsae(np.zeros((3, 3)), [3]) == {
        "ipsae_ab": 0.0, "ipsae_ba": 0.0, "ipsae_min": 0.0
    }


def test_ipsae_min_is_minimum_of_directions():
    # Asymmetric interface: A->B block confident (0), B->A block not (20).
    pae = np.full((4, 4), 20.0)
    pae[:2, 2:] = 0.0                 # A->B block confident
    out = cm.compute_ipsae(pae, [2, 2])
    assert out["ipsae_ab"] == pytest.approx(1.0)
    assert out["ipsae_ba"] == 0.0
    assert out["ipsae_min"] == 0.0


# ── compute_ipae ──────────────────────────────────────────────────────────

def test_ipae_mean_of_interchain_blocks():
    pae = np.zeros((4, 4))
    pae[:2, 2:] = 4.0                 # A->B block
    pae[2:, :2] = 6.0                 # B->A block
    # mean over both off-diagonal blocks = (4*4 + 6*4)/8 = 5.0
    assert cm.compute_ipae(pae, [2, 2]) == pytest.approx(5.0)


def test_ipae_size_mismatch_returns_zero():
    assert cm.compute_ipae(np.zeros((5, 5)), [2, 2]) == 0.0


def test_ipae_none_returns_zero():
    assert cm.compute_ipae(None, [2, 2]) == 0.0


# ── _pae_derived_metrics roll-up ──────────────────────────────────────────

def test_pae_derived_metrics_keys():
    pae = np.zeros((4, 4))
    out = cm._pae_derived_metrics(pae, [2, 2])
    assert set(out) == {"pae_mean", "ipae", "ipsae_ab", "ipsae_ba", "ipsae_min"}
    assert out["pae_mean"] == 0.0


# ── calc_d0: the TM-score normaliser ipSAE uses ───────────────────────────
# These pin the exact formula. The previous implementation used
# d0 = 0.5*sqrt(n), which is not ipSAE and gives materially different scores
# at the interface sizes this benchmark sees.

def test_calc_d0_matches_the_tm_score_formula():
    # L = 100 -> 1.24*(85)^(1/3) - 1.8
    expected = 1.24 * (85.0 ** (1.0 / 3.0)) - 1.8
    assert cm.calc_d0(100) == pytest.approx(expected, abs=1e-9)


def test_calc_d0_floors_short_interfaces_at_L_27():
    # Below 27 residues the cube-root term collapses, so L is clamped.
    assert cm.calc_d0(5) == pytest.approx(cm.calc_d0(27), abs=1e-12)
    assert cm.calc_d0(26) == pytest.approx(cm.calc_d0(27), abs=1e-12)


def test_calc_d0_never_drops_below_one_angstrom():
    assert cm.calc_d0(1) >= 1.0
    assert cm.calc_d0(27) >= 1.0


def test_calc_d0_increases_with_interface_size():
    vals = [cm.calc_d0(n) for n in (27, 50, 100, 400)]
    assert vals == sorted(vals)
    assert len(set(vals)) == len(vals)


def test_calc_d0_is_not_the_old_sqrt_heuristic():
    # Guards against a silent revert to d0 = 0.5*sqrt(n).
    assert cm.calc_d0(100) != pytest.approx(0.5 * np.sqrt(100), abs=1e-6)


# ── ipSAE with non-degenerate PAE ─────────────────────────────────────────
# An all-zero PAE gives 1.0 under any d0, so these use real error values.

def test_ipsae_uses_calc_d0_for_the_subset_size():
    # 30 partners, all at PAE 5.0, all below the 10 A cutoff.
    pae = np.full((1, 30), 5.0)
    d0 = cm.calc_d0(30)
    expected = 1.0 / (1.0 + (5.0 / d0) ** 2)
    score, idx = cm.compute_ipsae_one_direction(pae, cutoff=10.0)
    assert score == pytest.approx(expected, abs=1e-9)
    assert idx == 0


def test_ipsae_excludes_partners_above_the_cutoff():
    # 30 good partners at PAE 2, plus 70 at PAE 20 that must be ignored:
    # both the mean AND d0 must be computed on the surviving 30 only.
    row = np.concatenate([np.full(30, 2.0), np.full(70, 20.0)])
    d0 = cm.calc_d0(30)
    expected = 1.0 / (1.0 + (2.0 / d0) ** 2)
    score, _ = cm.compute_ipsae_one_direction(row.reshape(1, -1), cutoff=10.0)
    assert score == pytest.approx(expected, abs=1e-9)


def test_ipsae_takes_the_best_residue_not_the_mean():
    pae = np.vstack([np.full(30, 9.0), np.full(30, 1.0)])   # one poor, one good
    score, idx = cm.compute_ipsae_one_direction(pae, cutoff=10.0)
    best_alone, _ = cm.compute_ipsae_one_direction(pae[1:2, :], cutoff=10.0)
    assert idx == 1
    assert score == pytest.approx(best_alone, abs=1e-12)


def test_ipsae_falls_to_zero_when_no_partner_clears_the_cutoff():
    # Every residue scores 0, so the returned index is arbitrary (argmax of a
    # flat array) and carries no meaning — only the score does.
    score, _ = cm.compute_ipsae_one_direction(np.full((4, 6), 25.0), cutoff=10.0)
    assert score == 0.0


def test_ipsae_returns_sentinel_index_for_an_empty_matrix():
    score, idx = cm.compute_ipsae_one_direction(np.zeros((0, 6)), cutoff=10.0)
    assert score == 0.0
    assert idx == -1


def test_ipsae_is_insensitive_to_non_interacting_padding():
    # The property ipSAE exists to provide, and the one ipTM lacks: adding
    # residues that interact with nothing must not move the score.
    core = np.full((2, 25), 3.0)
    padded = np.hstack([core, np.full((2, 200), 30.0)])
    a, _ = cm.compute_ipsae_one_direction(core, cutoff=10.0)
    b, _ = cm.compute_ipsae_one_direction(padded, cutoff=10.0)
    assert a == pytest.approx(b, abs=1e-12)

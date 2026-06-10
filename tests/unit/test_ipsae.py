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

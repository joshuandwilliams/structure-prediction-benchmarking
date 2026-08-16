"""Unit tests for the shared analysis library.

Every analysis loads its data through load(), so a silent change in the
cleaning or filtering here moves numbers in all ten documents at once.
"""


import _analysis_common as ac
import numpy as np
import pandas as pd
import pytest

pytestmark = pytest.mark.local_unit


# ── Display helpers ───────────────────────────────────────────────────────

def test_pretty_wraps_before_a_parenthesised_variant():
    assert ac.pretty("boltz2") == "Boltz-2"
    assert ac.pretty("boltz2_constrained") == "Boltz-2\n(constr)"


def test_pretty_passes_unknown_names_through():
    assert ac.pretty("some_new_model") == "some_new_model"


def test_combo_label_joins_model_and_msa_on_one_line():
    assert ac.combo_label("boltz2", "no_msa") == "Boltz-2 / No MSA"
    assert "\n" not in ac.combo_label("boltz2_constrained", "msa")


def test_is_oracle_flags_only_the_reference_derived_variants():
    assert ac.is_oracle("boltz1_constrained")
    assert ac.is_oracle("boltz2_constrained")
    assert not ac.is_oracle("boltz2")
    assert not ac.is_oracle("af3")


# ── Statistics ────────────────────────────────────────────────────────────

def test_auc_is_one_for_a_perfect_separator():
    assert ac.auc([1.0, 2.0, 3.0, 4.0], [False, False, True, True]) == 1.0


def test_auc_is_half_for_an_uninformative_score():
    assert ac.auc([1.0, 1.0, 1.0, 1.0], [True, False, True, False]) == 0.5


def test_auc_inverts_for_a_reversed_separator():
    assert ac.auc([4.0, 3.0, 2.0, 1.0], [False, False, True, True]) == 0.0


def test_auc_is_nan_when_one_class_is_absent():
    assert np.isnan(ac.auc([1.0, 2.0], [True, True]))
    assert np.isnan(ac.auc([1.0, 2.0], [False, False]))


def test_auc_ignores_nan_scores():
    assert ac.auc([1.0, np.nan, 3.0], [False, True, True]) == 1.0


def test_spearman_is_one_for_a_monotonic_relationship():
    assert ac.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert ac.spearman([1, 2, 3, 4], [40, 30, 20, 10]) == pytest.approx(-1.0)


def test_spearman_ranks_rather_than_fitting_a_line():
    """Non-linear but monotonic still gives 1, which is why it is used."""
    assert ac.spearman([1, 2, 3, 4], [1, 4, 9, 16]) == pytest.approx(1.0)


# ── Cleaning ──────────────────────────────────────────────────────────────

def test_plddt_is_rescaled_onto_a_common_axis():
    """Boltz reports 0-1 and the AlphaFold lineage 0-100. Plotting them
    together requires one scale."""
    df = pd.DataFrame({"model": ["boltz2", "af3"], "avg_plddt": [0.85, 85.0]})
    out = ac._clean(df.copy())
    assert out["plddt"].tolist() == [85.0, 85.0]


def test_all_zero_pae_columns_are_blanked_per_model():
    """Chai-1 emits no PAE, and compute_metrics writes 0.0 rather than blank.
    Reading those zeros as real would drag its ipSAE to the floor."""
    df = pd.DataFrame({
        "model": ["chai1", "chai1", "boltz2"],
        "avg_plddt": [80.0, 80.0, 80.0],
        "ipsae_min": [0.0, 0.0, 0.4],
    })
    out = ac._clean(df.copy())
    assert out.loc[out["model"] == "chai1", "ipsae_min"].isna().all()
    assert out.loc[out["model"] == "boltz2", "ipsae_min"].tolist() == [0.4]


def test_a_genuine_zero_survives_when_the_model_has_other_values():
    df = pd.DataFrame({
        "model": ["boltz2", "boltz2"],
        "avg_plddt": [80.0, 80.0],
        "ipsae_min": [0.0, 0.4],
    })
    out = ac._clean(df.copy())
    assert out["ipsae_min"].tolist() == [0.0, 0.4]


# ── Ordering and layout ───────────────────────────────────────────────────

def _combo_frame():
    return pd.DataFrame({
        "model": ["boltz2", "boltz2", "af3", "af3"],
        "msa": ["no_msa", "msa", "no_msa", "msa"],
        ac.RA_COL: [20.0, 1.0, 2.0, 1.5],
    })


def test_combos_are_ordered_no_msa_block_first():
    combos = ac.ordered_combos(_combo_frame())
    flags = [f for f, _ in combos]
    assert flags == sorted(flags, key=lambda f: f != "no_msa")


def test_combos_rank_by_correct_count_then_median():
    combos = ac.ordered_combos(_combo_frame())
    no_msa = [m for f, m in combos if f == "no_msa"]
    assert no_msa[0] == "af3"          # 2.0 Å passes, boltz2's 20.0 Å does not


def test_block_positions_leave_a_gap_between_the_two_blocks():
    combos = [("no_msa", "a"), ("no_msa", "b"), ("msa", "a")]
    xs, n_no = ac.block_positions(combos, gap=1.0)
    assert n_no == 2
    assert xs[1] - xs[0] == 1.0        # within a block
    assert xs[2] - xs[1] == 2.0        # across the gap


# ── save_fig ──────────────────────────────────────────────────────────────

def test_save_fig_writes_into_the_calling_analysis_folder(tmp_path, monkeypatch):
    """Paths are relative to the CWD, which Quarto sets to the document's
    folder, so figures land beside the analysis that made them."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    monkeypatch.chdir(tmp_path)
    fig = plt.figure()
    ac.save_fig(fig, "thing", thesis=True)
    ac.save_fig(fig, "other")
    ac.save_fig(fig, "nested", subdir="sub")
    plt.close(fig)

    assert (tmp_path / "thesis-figures" / "thing.png").is_file()
    assert (tmp_path / "supplementary-figures" / "other.png").is_file()
    assert (tmp_path / "supplementary-figures" / "sub" / "nested.png").is_file()


# ── load ──────────────────────────────────────────────────────────────────

@pytest.mark.local_integration
def test_load_reads_the_committed_metrics(capsys):
    """Guards the wiring between data/metrics/ and DATA_DIR."""
    d = ac.load(verbose=False)
    for key in ("df", "pred", "rt", "sim", "near", "rel", "combos", "lookup"):
        assert key in d
    assert len(d["df"]) > 0
    assert d["df"]["pdb"].nunique() <= 18


@pytest.mark.local_integration
def test_load_filters_to_the_manifest_targets():
    d = ac.load(verbose=False)
    manifest = set(pd.read_csv(ac.MANIFEST, sep="\t")["pdb"])
    assert set(d["df"]["pdb"]) <= manifest


@pytest.mark.local_integration
def test_load_marks_correctness_on_every_prediction():
    d = ac.load(verbose=False)
    pred = d["pred"]
    assert pred["correct"].dtype == bool
    assert (pred["correct"] == (pred[ac.RA_COL] < ac.RA_EFF_THRESHOLD)).all()


@pytest.mark.local_integration
def test_load_uses_standalone_runtime_so_msa_time_is_charged():
    d = ac.load(verbose=False)
    rt = d["rt"]
    assert (rt["elapsed_min"] == rt["standalone_elapsed_s"] / 60.0).all()


@pytest.mark.local_integration
def test_verbose_load_reports_what_it_dropped(capsys):
    ac.load(verbose=True)
    out = capsys.readouterr().out
    assert "complexes" in out
    assert "combos" in out

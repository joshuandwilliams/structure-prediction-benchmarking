"""Interface classification against poses with known ground truth.

The classifier decides which of two receptor surfaces a prediction placed the
effector on. Testing that against a benchmark run is circular, because the run
is what the classifier is meant to describe, and it cannot run before the run
exists. Instead each pose here is built from a committed reference by a rigid
transform, so the expected label follows from the construction.

    boltz2   the reference complex, so the true AVR-Pik surface
    boltz1   the frame effector superposed on, so the AVR-Pia surface
    chai1    the true effector displaced 40 A, so neither surface

The predictor tags carry no meaning. Both scripts enumerate MODEL_MAP to find
files, so a pose has to borrow a real tag to be seen at all.

One hpc-tier test still runs the classifier over a real benchmark when one is
present, which is what catches a published tree the fixtures cannot model.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest
from _analysis_common import REPO_ROOT

sys.path.insert(0, str(Path(__file__).parent / "fixtures"))
from build_interface_fixtures import POSE_TAGS, TARGETS, build  # noqa: E402

pytestmark = pytest.mark.local_integration

FRAME = "6Q76"
NON_FRAME = [p for p in TARGETS if p != FRAME]

CORRECT, PIA, OTHER = (POSE_TAGS["correct"], POSE_TAGS["pia"],
                       POSE_TAGS["other"])


@pytest.fixture(scope="module")
def preds(tmp_path_factory):
    """The synthetic prediction set, built once for the module."""
    d = tmp_path_factory.mktemp("interface_preds")
    build(str(d))
    return d


def _run(script, *args, cwd):
    return subprocess.run(
        [sys.executable, str(REPO_ROOT / "bin" / script), *map(str, args)],
        capture_output=True, text=True, cwd=str(cwd))


@pytest.fixture(scope="module")
def classified(preds, tmp_path_factory):
    import pandas as pd
    out = tmp_path_factory.mktemp("classify") / "interfaces.csv"
    r = _run("classify_predicted_interface.py", "--pred-dir", preds,
             "--output", out, cwd=out.parent)
    assert r.returncode == 0, r.stderr
    return pd.read_csv(out)


def _site(df, pdb, model):
    rows = df[(df["pdb"] == pdb) & (df["model"] == model)]
    assert len(rows) == 1, f"expected one row for {pdb}/{model}, got {len(rows)}"
    return rows.iloc[0]


# ── classify_predicted_interface ──────────────────────────────────────────

def test_it_emits_a_row_per_prediction(classified):
    assert len(classified) == len(TARGETS) * len(POSE_TAGS)
    assert {"pdb", "model", "msa", "site", "jaccard_true",
            "jaccard_pia"} <= set(classified)


def test_a_native_pose_is_called_the_correct_site(classified):
    """The reference complex is on the true surface by definition, so anything
    other than 'correct site' here means the classifier is broken."""
    for pdb in TARGETS:
        row = _site(classified, pdb, CORRECT)
        assert row["site"] == "correct site", f"{pdb}: got {row['site']}"
        assert row["jaccard_true"] == pytest.approx(1.0)


def test_a_pose_on_the_avr_pia_surface_is_called_that(classified):
    """This is the discrimination the whole analysis rests on. The effector is
    on the receptor, contacting a real patch, just the wrong one."""
    for pdb in NON_FRAME:
        row = _site(classified, pdb, PIA)
        assert row["site"] == "AVR-Pia site", f"{pdb}: got {row['site']}"
        assert row["jaccard_pia"] > row["jaccard_true"]


def test_a_pose_on_neither_surface_is_called_other(classified):
    for pdb in TARGETS:
        row = _site(classified, pdb, OTHER)
        assert row["site"] == "other", f"{pdb}: got {row['site']}"
        assert row["jaccard_true"] == 0.0
        assert row["jaccard_pia"] == 0.0


def test_the_frame_target_cannot_distinguish_the_two_sites(classified):
    """6Q76's true site IS the AVR-Pia site, so a pose there scores 1.0 against
    both. This is why the frame is excluded from the published tally."""
    row = _site(classified, FRAME, CORRECT)
    assert row["jaccard_true"] == pytest.approx(1.0)
    assert row["jaccard_pia"] == pytest.approx(1.0)


def test_every_label_is_one_of_the_three(classified):
    assert set(classified["site"]) <= {"correct site", "AVR-Pia site", "other"}


def test_rmsd_is_recomputed_on_the_classified_structure(classified):
    """The site call and the RMSD must describe one molecule. A native pose has
    to read 0 A, and the displaced pose has to read back the offset it was
    built with, which no cached per-run value could produce."""
    assert "ra_eff_this_structure" in classified
    for pdb in TARGETS:
        assert _site(classified, pdb, CORRECT)["ra_eff_this_structure"] == \
            pytest.approx(0.0, abs=1e-6)
        assert _site(classified, pdb, OTHER)["ra_eff_this_structure"] == \
            pytest.approx(40.0, abs=0.5)


# ── common_wrong_interface ────────────────────────────────────────────────

@pytest.fixture(scope="module")
def common(preds, tmp_path_factory):
    import pandas as pd
    out = tmp_path_factory.mktemp("common") / "common.csv"
    r = _run("common_wrong_interface.py", "--model", "boltz2", "--msa",
             "no_msa", "--frame", FRAME, "--system", "Pik",
             "--pred-dir", preds, "--output", out, cwd=out.parent)
    assert r.returncode == 0, r.stderr
    return pd.read_csv(out)


def test_wrong_interface_reports_distances_to_both_sites(common):
    assert {"pdb", "ra_eff", "d_pred_to_own_true",
            "d_pred_to_frame_eff"} <= set(common)
    assert len(common) == len(TARGETS)


def test_the_two_surfaces_are_far_apart_in_every_target(common):
    """The argument rests on the AVR-Pik and AVR-Pia surfaces being genuinely
    distinct rather than a marginal shift."""
    non_frame = common[common["pdb"] != FRAME]
    assert len(non_frame) > 0
    assert (non_frame["d_true_to_frame_eff"] > 20.0).all()


def test_the_frame_sits_on_top_of_itself(common):
    """6Q76 defines the AVR-Pia position, so its distance to it must be zero.
    A non-zero value means the frame superposition drifted."""
    row = common[common["pdb"] == FRAME].iloc[0]
    assert row["d_true_to_frame_eff"] == pytest.approx(0.0, abs=1e-6)


# ── hpc tier ──────────────────────────────────────────────────────────────

BENCH = Path(os.environ.get("SPB_BENCHMARKS_DIR",
                            REPO_ROOT / "experiments" / "benchmarks"))


def _has_results():
    if not BENCH.is_dir():
        return False
    return any((BENCH / p / f"{p}_benchmark_results" / "best_models").is_dir()
               for p in os.listdir(BENCH))


@pytest.mark.hpc
@pytest.mark.skipif(not _has_results(),
                    reason=f"no published best_models/ under {BENCH}")
def test_it_runs_over_a_real_benchmark(tmp_path):
    """The fixtures pin the logic. This pins that the published tree layout is
    still what the script expects to walk."""
    import pandas as pd
    out = tmp_path / "interfaces.csv"
    r = _run("classify_predicted_interface.py", "--output", out, cwd=tmp_path)
    assert r.returncode == 0, r.stderr
    df = pd.read_csv(out)
    assert len(df) > 0
    assert set(df["site"]) <= {"correct site", "AVR-Pia site", "other"}

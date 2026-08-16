"""Unit tests for the CSV aggregation scripts.

combine_metrics picks the representative prediction per (model, msa, pdb),
which is the selection every accuracy figure rests on. build_per_prediction_csv
keeps all 25 instead, and build_runtime_csv collects the timings.
"""

import csv
import sys

import build_per_prediction_csv as bpp
import build_runtime_csv as brc
import combine_metrics as cmb
import pytest

pytestmark = pytest.mark.local_unit

COLUMNS = ["model", "model_name", "avg_plddt", "iptm", "ptm",
           "rmsd_effector_receptor_aligned", "pdb_path"]


def _all_metrics(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLUMNS})


def _bench(tmp_path, pdb, rows, suffix="_benchmark_results"):
    _all_metrics(tmp_path / pdb / f"{pdb}{suffix}" / "all_metrics.csv", rows)


def _run(module, tmp_path, out, *extra):
    sys.argv = [module.__name__, "--benchmarks-dir", str(tmp_path),
                "--output", str(out), *extra]
    module.main()
    with open(out, newline="") as fh:
        return list(csv.DictReader(fh))


# ── combine_metrics: the selection rule ───────────────────────────────────

def test_the_most_confident_prediction_is_kept(tmp_path, capsys):
    """Highest avg_plddt wins, not the lowest RMSD. Selecting on RMSD needs
    the answer and reports an oracle no pipeline could reach."""
    _bench(tmp_path, "6G10", [
        {"model": "boltz2", "model_name": "a", "avg_plddt": "70",
         "rmsd_effector_receptor_aligned": "1.0"},
        {"model": "boltz2", "model_name": "b", "avg_plddt": "90",
         "rmsd_effector_receptor_aligned": "20.0"},
    ])
    rows = _run(cmb, tmp_path, tmp_path / "out.csv")
    assert len(rows) == 1
    assert rows[0]["model_name"] == "b"                      # confident, not close
    assert rows[0]["rmsd_effector_receptor_aligned"] == "20.0"


def test_each_model_msa_pdb_combination_gets_one_row(tmp_path, capsys):
    _bench(tmp_path, "6G10", [
        {"model": "boltz2", "model_name": "a", "avg_plddt": "70",
         "rmsd_effector_receptor_aligned": "3.0"},
        {"model": "boltz2_msa", "model_name": "b", "avg_plddt": "80",
         "rmsd_effector_receptor_aligned": "3.0"},
        {"model": "af3_nomsa", "model_name": "c", "avg_plddt": "60",
         "rmsd_effector_receptor_aligned": "3.0"},
    ])
    _bench(tmp_path, "5A6W", [
        {"model": "boltz2", "model_name": "d", "avg_plddt": "75",
         "rmsd_effector_receptor_aligned": "3.0"},
    ])
    rows = _run(cmb, tmp_path, tmp_path / "out.csv")
    assert len(rows) == 4
    assert {(r["model"], r["msa"], r["pdb"]) for r in rows} == {
        ("boltz2", "no_msa", "6G10"), ("boltz2", "msa", "6G10"),
        ("af3", "no_msa", "6G10"), ("boltz2", "no_msa", "5A6W")}


def test_template_variants_stay_separate_from_their_twins(tmp_path, capsys):
    _bench(tmp_path, "6G10", [
        {"model": "boltz2", "model_name": "free", "avg_plddt": "70",
         "rmsd_effector_receptor_aligned": "3.0"},
        {"model": "boltz2_template", "model_name": "tmpl", "avg_plddt": "80",
         "rmsd_effector_receptor_aligned": "3.0"},
        {"model": "boltz2_msa_template", "model_name": "tmpl_msa",
         "avg_plddt": "85", "rmsd_effector_receptor_aligned": "3.0"},
    ])
    rows = _run(cmb, tmp_path, tmp_path / "out.csv")
    by = {(r["model"], r["msa"]): r["model_name"] for r in rows}
    assert by[("boltz2", "no_msa")] == "free"
    assert by[("boltz2_template", "no_msa")] == "tmpl"
    assert by[("boltz2_template", "msa")] == "tmpl_msa"


def test_rows_without_a_usable_score_are_dropped(tmp_path, capsys):
    _bench(tmp_path, "6G10", [
        {"model": "boltz2", "model_name": "blank", "avg_plddt": "",
         "rmsd_effector_receptor_aligned": "3.0"},
        {"model": "boltz2", "model_name": "zero", "avg_plddt": "0",
         "rmsd_effector_receptor_aligned": "3.0"},
        {"model": "boltz2", "model_name": "real", "avg_plddt": "55",
         "rmsd_effector_receptor_aligned": "3.0"},
    ])
    rows = _run(cmb, tmp_path, tmp_path / "out.csv")
    assert len(rows) == 1 and rows[0]["model_name"] == "real"


def test_a_combo_with_no_valid_rmsd_is_dropped(tmp_path, capsys):
    """A row needs both a confidence score and an ra_eff. Without the RMSD
    there is nothing to score the prediction against."""
    _bench(tmp_path, "6G10", [
        {"model": "boltz2", "model_name": "no_rmsd", "avg_plddt": "70"}])
    rows = _run(cmb, tmp_path, tmp_path / "out.csv")
    assert rows == []
    assert "no valid ra_eff" in capsys.readouterr().out


def test_a_target_with_no_metrics_file_is_reported_not_fatal(tmp_path, capsys):
    _bench(tmp_path, "6G10", [{"model": "boltz2", "model_name": "a",
                               "avg_plddt": "70",
                               "rmsd_effector_receptor_aligned": "3.0"}])
    (tmp_path / "5A6W").mkdir()
    rows = _run(cmb, tmp_path, tmp_path / "out.csv")
    assert len(rows) == 1
    assert "5A6W" in capsys.readouterr().out


def test_missing_benchmarks_dir_exits(tmp_path):
    sys.argv = ["combine_metrics.py", "--benchmarks-dir",
                str(tmp_path / "nope"), "--output", str(tmp_path / "o.csv")]
    with pytest.raises(SystemExit):
        cmb.main()


# ── build_per_prediction_csv: keeps everything ────────────────────────────

def test_every_prediction_is_kept_not_just_the_best(tmp_path, capsys):
    _bench(tmp_path, "6G10", [
        {"model": "boltz2", "model_name": f"s{i}", "avg_plddt": str(50 + i),
         "rmsd_effector_receptor_aligned": "3.0"} for i in range(5)])
    rows = _run(bpp, tmp_path, tmp_path / "pp.csv")
    assert len(rows) == 5
    assert {r["model_name"] for r in rows} == {f"s{i}" for i in range(5)}


def test_rows_without_an_rmsd_are_dropped(tmp_path, capsys):
    _bench(tmp_path, "6G10", [
        {"model": "boltz2", "model_name": "ok", "avg_plddt": "70",
         "rmsd_effector_receptor_aligned": "3.0"},
        {"model": "boltz2", "model_name": "no_rmsd", "avg_plddt": "70"},
    ])
    rows = _run(bpp, tmp_path, tmp_path / "pp.csv")
    assert [r["model_name"] for r in rows] == ["ok"]


# ── build_runtime_csv ─────────────────────────────────────────────────────

def test_runtime_rows_are_collected_per_target(tmp_path, capsys):
    for pdb in ("6G10", "5A6W"):
        d = tmp_path / pdb / f"{pdb}_benchmark_results"
        d.mkdir(parents=True)
        (d / "predictor_runtime_stats.csv").write_text(
            "model,status,elapsed_s,standalone_elapsed_s,peak_rss_gb\n"
            "boltz2,COMPLETED,300,300,5.00\n"
            "boltz2_msa,COMPLETED,300,2364,5.00\n")
    rows = _run(brc, tmp_path, tmp_path / "rt.csv")
    assert len(rows) == 4
    assert {r["pdb"] for r in rows} == {"6G10", "5A6W"}
    msa = [r for r in rows if r["model"] == "boltz2" and r["msa"] == "msa"]
    assert msa and msa[0]["standalone_elapsed_s"] == "2364"

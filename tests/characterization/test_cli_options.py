"""Command-line options and remaining branches across bin/.

Covers the tier filters, the network-backed date cache, and the structural
guards in the template extractor.
"""

import csv
import json
import sys
from urllib.error import URLError

import build_per_prediction_csv as bpp
import build_runtime_csv as brc
import extract_effector_template as eet
import plot_training_set_membership as ptsm
import pytest
from _analysis_common import REPO_ROOT

pytestmark = pytest.mark.local_unit

REFS = REPO_ROOT / "data" / "complexes_for_benchmarking"
COLUMNS = ["model", "model_name", "avg_plddt", "iptm", "ptm",
           "rmsd_effector_receptor_aligned", "pdb_path"]


def _bench(root, pdb, rows):
    d = root / pdb / f"{pdb}_benchmark_results"
    d.mkdir(parents=True)
    with open(d / "all_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLUMNS})


def _manifest(path, entries):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pdb", "tier", "system"])
        w.writerows(entries)


ROW = {"model": "boltz2", "model_name": "m", "avg_plddt": "70",
       "rmsd_effector_receptor_aligned": "3.0"}


# ── Tier filtering ────────────────────────────────────────────────────────

def test_per_prediction_tier_filter_keeps_only_matching_targets(tmp_path, capsys):
    _bench(tmp_path, "6G10", [ROW])
    _bench(tmp_path, "7XVG", [ROW])
    man = tmp_path / "manifest.tsv"
    _manifest(man, [("6G10", 1, "Pik"), ("7XVG", 2, "Sr35")])

    out = tmp_path / "pp.csv"
    sys.argv = ["build_per_prediction_csv.py", "--benchmarks-dir", str(tmp_path),
                "--manifest", str(man), "--tier", "1", "--output", str(out)]
    bpp.main()
    rows = list(csv.DictReader(open(out)))
    assert {r["pdb"] for r in rows} == {"6G10"}


def test_per_prediction_tier_filter_needs_a_manifest(tmp_path):
    _bench(tmp_path, "6G10", [ROW])
    sys.argv = ["build_per_prediction_csv.py", "--benchmarks-dir", str(tmp_path),
                "--manifest", str(tmp_path / "nope.tsv"), "--tier", "1",
                "--output", str(tmp_path / "pp.csv")]
    with pytest.raises(SystemExit):
        bpp.main()


def test_runtime_tier_filter_keeps_only_matching_targets(tmp_path, capsys):
    for pdb in ("6G10", "7XVG"):
        d = tmp_path / pdb / f"{pdb}_benchmark_results"
        d.mkdir(parents=True)
        (d / "predictor_runtime_stats.csv").write_text(
            "model,status,elapsed_s,standalone_elapsed_s,peak_rss_gb\n"
            "boltz2,COMPLETED,300,300,5.00\n")
    man = tmp_path / "manifest.tsv"
    _manifest(man, [("6G10", 1, "Pik"), ("7XVG", 2, "Sr35")])

    out = tmp_path / "rt.csv"
    sys.argv = ["build_runtime_csv.py", "--benchmarks-dir", str(tmp_path),
                "--manifest", str(man), "--tier", "1", "--output", str(out)]
    brc.main()
    assert {r["pdb"] for r in csv.DictReader(open(out))} == {"6G10"}


# ── The release-date cache ────────────────────────────────────────────────

def test_a_missing_pdb_is_fetched_and_written_back(tmp_path, monkeypatch, capsys):
    """The cache is self-maintaining, so a new target needs network once and
    never again."""
    dates = tmp_path / "dates.csv"
    dates.write_text("pdb,release_date\n6G10,2018-06-13\n")
    bench = tmp_path / "bench"
    for pdb in ("6G10", "9IP6"):
        (bench / pdb).mkdir(parents=True)

    monkeypatch.setattr(ptsm, "fetch_release_date", lambda p: "2025-01-01")
    sys.argv = ["plot_training_set_membership.py", "--benchmarks-dir", str(bench),
                "--dates-csv", str(dates), "--outdir", str(tmp_path / "p")]
    ptsm.main()

    cache = {r["pdb"]: r["release_date"] for r in csv.DictReader(open(dates))}
    assert cache["9IP6"] == "2025-01-01"


def test_a_failed_fetch_returns_none_rather_than_raising(monkeypatch, capsys):
    def boom(url, timeout=0):
        raise URLError("no network")
    monkeypatch.setattr(ptsm.urllib.request, "urlopen", boom)
    assert ptsm.fetch_release_date("9IP6") is None


def test_a_successful_fetch_parses_the_release_date(monkeypatch):
    class FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self):
            return json.dumps({"rcsb_accession_info":
                               {"initial_release_date": "2018-06-13T00:00:00Z"}}).encode()

    monkeypatch.setattr(ptsm.urllib.request, "urlopen",
                        lambda url, timeout=0: FakeResponse())
    monkeypatch.setattr(ptsm.json, "load",
                        lambda fh: json.loads(fh.read().decode()))
    assert ptsm.fetch_release_date("6G10") == "2018-06-13"


# ── extract_effector_template guards ──────────────────────────────────────

def test_template_extraction_handles_a_chain_already_named_a(tmp_path, monkeypatch):
    """No relabelling is needed when the effector is already chain A, and the
    code must not drop it while pruning."""
    monkeypatch.chdir(tmp_path)
    eet.extract_effector_template(str(REFS / "6G10.pdb"), "A")
    import gemmi
    st = gemmi.read_structure(str(tmp_path / "effector_template.cif"))
    assert [c.name for c in st[0]] == ["A"]


def test_template_extraction_rejects_a_missing_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(Exception):
        eet.extract_effector_template(str(tmp_path / "nope.pdb"), "B")

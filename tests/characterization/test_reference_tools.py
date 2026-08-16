"""Unit tests for the scripts that operate on the reference complexes.

These run against the committed references rather than synthetic input, since
their contract is with real crystal structures: multiple chains, insertion
codes, unmodelled termini and heteroatoms.
"""

import csv
import shutil
import sys

import compute_reference_similarity as crs
import extract_effector_template as eet
import plot_training_set_membership as ptsm
import pytest
from _analysis_common import MANIFEST, REPO_ROOT

pytestmark = pytest.mark.local_integration

REFS = REPO_ROOT / "data" / "complexes_for_benchmarking"

# compute_reference_similarity shells out to EMBOSS needle for the pairwise
# alignment. The container ships Python packages only, so these skip there
# rather than failing, in the same way the suite treats a missing esm or gemmi.
needs_needle = pytest.mark.skipif(
    shutil.which("needle") is None,
    reason="EMBOSS needle not on PATH")


@pytest.fixture
def two_refs(tmp_path):
    """A two-entry reference directory and matching manifest."""
    refs = tmp_path / "refs"
    refs.mkdir()
    for pdb in ("6G10", "6G11"):
        shutil.copy(REFS / f"{pdb}.pdb", refs / f"{pdb}.pdb")

    rows = [r for r in csv.DictReader(open(MANIFEST), delimiter="\t")
            if r["pdb"] in {"6G10", "6G11"}]
    man = tmp_path / "manifest.tsv"
    with open(man, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return refs, man


# ── extract_effector_template ─────────────────────────────────────────────

def test_template_keeps_only_the_effector_relabelled_to_a(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    eet.extract_effector_template(str(REFS / "6G10.pdb"), "B")

    import gemmi
    st = gemmi.read_structure(str(tmp_path / "effector_template.cif"))
    assert [c.name for c in st[0]] == ["A"]

    ref = gemmi.read_structure(str(REFS / "6G10.pdb"))
    n_eff = sum(1 for r in ref[0]["B"] if r.find_atom("CA", "*"))
    assert sum(1 for r in st[0]["A"] if r.find_atom("CA", "*")) == n_eff


def test_template_populates_the_fields_af3_requires(tmp_path, monkeypatch):
    """AF3's parser rejects a template without these, and gemmi does not
    write them when the source is a PDB."""
    monkeypatch.chdir(tmp_path)
    eet.extract_effector_template(str(REFS / "6G10.pdb"), "B")
    text = (tmp_path / "effector_template.cif").read_text()
    for field in ("_entity_poly_seq", "_entity_poly.pdbx_seq_one_letter_code",
                  "_atom_site.label_seq_id", "_pdbx_audit_revision_history"):
        assert field in text


def test_template_writes_both_pdb_and_cif(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    eet.extract_effector_template(str(REFS / "6G10.pdb"), "B")
    assert (tmp_path / "effector_template.pdb").is_file()
    assert (tmp_path / "effector_template.cif").is_file()


def test_template_exits_on_an_absent_chain(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        eet.extract_effector_template(str(REFS / "6G10.pdb"), "Z")


def test_template_main_runs_from_argv(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    sys.argv = ["extract_effector_template.py", str(REFS / "6G10.pdb"), "B"]
    eet.main()
    assert (tmp_path / "effector_template.cif").is_file()


def test_template_main_rejects_the_wrong_argument_count(monkeypatch):
    sys.argv = ["extract_effector_template.py"]
    with pytest.raises(SystemExit):
        eet.main()


# ── compute_reference_similarity ──────────────────────────────────────────

@needs_needle
def test_similarity_emits_one_row_per_pair(two_refs, tmp_path, capsys):
    refs, man = two_refs
    out = tmp_path / "sim.csv"
    sys.argv = ["compute_reference_similarity.py", "--refs-dir", str(refs),
                "--manifest", str(man), "--output", str(out)]
    crs.main()
    rows = list(csv.DictReader(open(out)))
    assert len(rows) == 1                      # one pair from two entries
    assert {rows[0]["pdb_a"], rows[0]["pdb_b"]} == {"6G10", "6G11"}


@needs_needle
def test_similarity_reports_both_identity_denominators(two_refs, tmp_path, capsys):
    """Percent identity has no single convention, so both are emitted rather
    than one being left implicit."""
    refs, man = two_refs
    out = tmp_path / "sim.csv"
    sys.argv = ["compute_reference_similarity.py", "--refs-dir", str(refs),
                "--manifest", str(man), "--output", str(out)]
    crs.main()
    row = next(iter(csv.DictReader(open(out))))
    for col in ("receptor_identity", "receptor_identity_over_alignment",
                "effector_identity", "effector_identity_over_alignment"):
        assert 0.0 <= float(row[col]) <= 100.0


@needs_needle
def test_similarity_separates_substitutions_from_gaps(two_refs, tmp_path, capsys):
    """These are crystal constructs, so a difference in modelled boundaries
    is a gap rather than a substitution."""
    refs, man = two_refs
    out = tmp_path / "sim.csv"
    sys.argv = ["compute_reference_similarity.py", "--refs-dir", str(refs),
                "--manifest", str(man), "--output", str(out)]
    crs.main()
    row = next(iter(csv.DictReader(open(out))))
    assert int(row["receptor_substitutions"]) >= 0
    assert int(row["receptor_gap_positions"]) >= 0


@needs_needle
def test_a_boundary_difference_is_a_gap_not_a_substitution(two_refs, tmp_path,
                                                           capsys):
    """6G10 and 6G11 carry the same Pikp-1 HMA sequence but are modelled to
    different boundaries. Identity therefore falls below 100 while the
    substitution count stays at zero, which is the distinction the separate
    gap column exists to preserve."""
    refs, man = two_refs
    out = tmp_path / "sim.csv"
    sys.argv = ["compute_reference_similarity.py", "--refs-dir", str(refs),
                "--manifest", str(man), "--output", str(out)]
    crs.main()
    row = next(iter(csv.DictReader(open(out))))
    assert int(row["receptor_substitutions"]) == 0
    assert int(row["receptor_gap_positions"]) > 0
    assert float(row["receptor_identity"]) < 100.0


# ── plot_training_set_membership ──────────────────────────────────────────

def test_exposure_counts_use_release_not_deposit_date(tmp_path, capsys):
    """A structure deposited before a cutoff but released after it was not
    available for training."""
    dates = tmp_path / "dates.csv"
    dates.write_text("pdb,release_date,deposit_date\n"
                     "6G10,2018-06-13,2018-03-20\n"
                     "9IP6,2025-01-01,2024-07-10\n")
    bench = tmp_path / "bench"
    for pdb in ("6G10", "9IP6"):
        (bench / pdb).mkdir(parents=True)

    sys.argv = ["plot_training_set_membership.py",
                "--benchmarks-dir", str(bench), "--dates-csv", str(dates),
                "--outdir", str(tmp_path / "plots")]
    ptsm.main()
    assert (tmp_path / "plots" / "training_set_membership.png").is_file()


def test_exposure_cutoffs_are_ordered_and_documented():
    """Every cutoff carries a source, so the table can be defended."""
    assert set(ptsm.MODEL_CUTOFFS) >= {
        "AlphaFold2-Multimer", "ColabFold", "Chai-1", "AlphaFold3",
        "Boltz-1", "Boltz-2", "ESMFold2"}
    for date, source in ptsm.MODEL_CUTOFFS.values():
        assert len(date) == 10 and source

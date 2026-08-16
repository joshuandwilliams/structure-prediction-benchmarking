"""Unit tests for the small helpers each pipeline stage depends on.

Covers three things the two benchmark runs cannot be trusted without:

* ``combine_metrics.map_model`` — turns a predictor directory name into the
  (model, msa) pair every figure groups by. A silent mislabel here moves
  points between boxes rather than erroring.
* ``combine_metrics.to_float`` — the numeric coercion in front of the
  selection key, where a mis-parsed blank would sort as a real value.
* ``af3_input`` — the JSON handed to AlphaFold 3, which must be
  template-free in both modes.

``extract_sequences`` is exercised through its gemmi path, and the constraint
extractors and YAML validator through their CLIs, matching how Nextflow calls
them.
"""

import json
import sys
from pathlib import Path

import af3_input
import combine_metrics as cmb
import extract_constraints_boltz1
import extract_constraints_boltz2
import pytest

pytestmark = pytest.mark.local_unit

REPO_ROOT = Path(__file__).resolve().parents[2]
BIN = REPO_ROOT / "bin"


# ── map_model ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("name,expected", [
    ("boltz1", ("boltz1", "no_msa")),
    ("boltz1_msa", ("boltz1", "msa")),
    ("boltz1_constrained", ("boltz1_constrained", "no_msa")),
    ("boltz2", ("boltz2", "no_msa")),
    ("boltz2_msa", ("boltz2", "msa")),
    ("boltz2_constrained", ("boltz2_constrained", "no_msa")),
    ("chai1", ("chai1", "no_msa")),
    ("af2m", ("af2m", "msa")),
    ("af3", ("af3", "msa")),
    ("af3_nomsa", ("af3", "no_msa")),
    ("colabfold", ("colabfold", "msa")),
    ("colabfold_nomsa", ("colabfold", "no_msa")),
    ("esmfold2", ("esmfold2", "no_msa")),
])
def test_map_model_covers_every_run1_arm(name, expected):
    assert cmb.map_model(name) == expected


@pytest.mark.parametrize("name,expected", [
    ("boltz2_template", ("boltz2_template", "no_msa")),
    ("boltz2_msa_template", ("boltz2_template", "msa")),
    ("boltz2_constrained_template", ("boltz2_constrained_template", "no_msa")),
])
def test_map_model_covers_every_run2_template_arm(name, expected):
    """The suffix heuristic keys off a TRAILING _msa/_nomsa, so
    'boltz2_msa_template' is only labelled msa because it is listed
    explicitly. Without that entry it silently becomes no_msa."""
    assert cmb.map_model(name) == expected


@pytest.mark.local_integration
def test_every_arm_in_main_nf_is_mapped():
    """MODEL_MAP must not drift behind main.nf's ALL_MODELS."""
    text = (REPO_ROOT / "main.nf").read_text()
    block = text.split("def ALL_MODELS = [", 1)[1].split("]", 1)[0]
    variants = [ln.split("'")[1] for ln in block.splitlines() if "'" in ln]
    assert len(variants) == 16
    unmapped = [a for a in variants if a not in cmb.MODEL_MAP]
    assert unmapped == [], f"variants missing from MODEL_MAP: {unmapped}"


def test_map_model_falls_back_for_an_unknown_name():
    assert cmb.map_model("newmodel_nomsa") == ("newmodel", "no_msa")
    assert cmb.map_model("newmodel_msa") == ("newmodel", "msa")
    assert cmb.map_model("newmodel") == ("newmodel", "no_msa")


def test_template_arms_stay_distinct_from_their_twins():
    """Run 2 must not collapse onto Run 1 in the aggregated table — the
    pairing between them is the whole point of the template comparison."""
    twins = [("boltz2", "boltz2_template"),
             ("boltz2_msa", "boltz2_msa_template"),
             ("boltz2_constrained", "boltz2_constrained_template")]
    for free, templated in twins:
        assert cmb.map_model(free)[0] != cmb.map_model(templated)[0]
        assert cmb.map_model(free)[1] == cmb.map_model(templated)[1]


# ── to_float ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("raw", ["", "  ", "NA", "na", "nan", "None", "null", None])
def test_to_float_treats_blanks_and_sentinels_as_missing(raw):
    assert cmb.to_float(raw) is None


@pytest.mark.parametrize("raw,expected", [
    ("0", 0.0), ("1.5", 1.5), ("-3.25", -3.25), (" 88.4 ", 88.4), ("1e2", 100.0),
])
def test_to_float_parses_real_numbers(raw, expected):
    assert cmb.to_float(raw) == pytest.approx(expected)


def test_to_float_rejects_non_numeric_text():
    assert cmb.to_float("high") is None


def test_zero_is_missing_distinct_from_none():
    """0.0 must not be confused with 'missing': the selection rule skips
    non-positive pLDDT, so conflating them would change which structure is
    published as a model's representative."""
    assert cmb.to_float("0") == 0.0
    assert cmb.to_float("") is None


# ── af3_input ─────────────────────────────────────────────────────────────

def _run_af3_input(tmp_path, mode, monkeypatch=None):
    """Call af3_input's main in-process so coverage records it."""
    out = tmp_path / "input.json"
    argv = ["af3_input.py", "MKFLVAAA", "GTALPPWW", mode, "--output", str(out)]
    old = sys.argv
    sys.argv = argv
    try:
        af3_input.main()
    finally:
        sys.argv = old
    return json.loads(out.read_text())


@pytest.mark.parametrize("mode", ["full", "nomsa"])
def test_af3_input_has_the_two_chains_and_five_seeds(tmp_path, mode):
    d = _run_af3_input(tmp_path, mode)
    assert d["dialect"] == "alphafold3"
    assert d["modelSeeds"] == [42, 123, 456, 789, 1024]
    chains = [s["protein"] for s in d["sequences"]]
    assert [c["id"] for c in chains] == ["A", "B"]
    assert chains[0]["sequence"] == "MKFLVAAA"
    assert chains[1]["sequence"] == "GTALPPWW"


@pytest.mark.parametrize("mode", ["full", "nomsa"])
def test_af3_input_never_injects_a_template(tmp_path, mode):
    """Both AF3 variants are template-free. In nomsa mode the field is present
    and empty; in full mode AF3 runs its own search, disabled in the module
    by --max_template_date."""
    d = _run_af3_input(tmp_path, mode)
    for s in d["sequences"]:
        assert s["protein"].get("templates", []) == []


def test_af3_nomsa_blanks_both_msa_fields(tmp_path):
    d = _run_af3_input(tmp_path, "nomsa")
    for s in d["sequences"]:
        assert s["protein"]["unpairedMsa"] == ""
        assert s["protein"]["pairedMsa"] == ""


def test_af3_full_mode_omits_msa_fields_so_the_pipeline_builds_them(tmp_path):
    d = _run_af3_input(tmp_path, "full")
    for s in d["sequences"]:
        assert "unpairedMsa" not in s["protein"]
        assert "pairedMsa" not in s["protein"]


def test_af3_input_rejects_an_unknown_mode(tmp_path):
    sys.argv = ["af3_input.py", "MKFL", "GTAL", "sideways",
                "--output", str(tmp_path / "x.json")]
    with pytest.raises(SystemExit):
        af3_input.main()


# ── extract_sequences ─────────────────────────────────────────────────────

def test_extract_sequences_reads_both_chains(make_pdb_file):
    gemmi = pytest.importorskip("gemmi", reason="extract_sequences needs gemmi")  # noqa: F841
    import extract_sequences as es

    atoms = ([{"chain": "A", "resseq": i + 1, "resname": r, "x": float(i),
               "y": 0.0, "z": 0.0}
              for i, r in enumerate(["MET", "LYS", "PHE", "LEU"])] +
             [{"chain": "B", "resseq": i + 1, "resname": r, "x": float(i),
               "y": 20.0, "z": 0.0}
              for i, r in enumerate(["GLY", "THR", "ALA"])])
    pdb = make_pdb_file(atoms, name="seqtest.pdb")

    # extract_sequences returns the chain list itself; main() wraps it in the
    # JSON envelope the workflow parses.
    by_id = {c["id"]: c for c in es.extract_sequences(str(pdb))}
    assert by_id["A"]["sequence"] == "MKFL"
    assert by_id["A"]["length"] == 4
    assert by_id["B"]["sequence"] == "GTA"
    assert by_id["B"]["length"] == 3


def test_extract_sequences_maps_selenomethionine_to_methionine(make_pdb_file):
    pytest.importorskip("gemmi", reason="extract_sequences needs gemmi")
    import extract_sequences as es

    atoms = [{"chain": "A", "resseq": i + 1, "resname": r, "x": float(i),
              "y": 0.0, "z": 0.0}
             for i, r in enumerate(["MSE", "LYS", "PHE"])]
    pdb = make_pdb_file(atoms, name="mse.pdb")
    chains = es.extract_sequences(str(pdb))
    assert chains[0]["sequence"] == "MKF"


# ── constraint extractors (CLI, as Nextflow calls them) ───────────────────

def _extract(module, pdb, *args, capsys=None):
    """Run an extractor's main in-process and return its stdout."""
    sys.argv = [module.__name__, str(pdb), "B", "C", *map(str, args)]
    module.main()
    return capsys.readouterr().out


def test_boltz1_extractor_emits_pocket_only_at_max_distance_6(
        make_pdb_file, two_chain_atoms, capsys):
    pdb = make_pdb_file(two_chain_atoms, name="c1.pdb")
    out = _extract(extract_constraints_boltz1, pdb, 8.0, 6.0, capsys=capsys)
    assert "- pocket:" in out
    assert "- contact:" not in out           # Boltz-1's schema rejects these
    assert "max_distance: 6.0" in out


def test_boltz1_extractor_emits_whatever_max_distance_it_is_given(
        make_pdb_file, two_chain_atoms, capsys):
    """The 6.0 requirement is enforced downstream by validate_boltz_yaml.py,
    not here, so this only pins that the value passes through."""
    pdb = make_pdb_file(two_chain_atoms, name="c1b.pdb")
    out = _extract(extract_constraints_boltz1, pdb, 8.0, 8.0, capsys=capsys)
    assert "max_distance: 8.0" in out


@pytest.mark.parametrize("module,args", [
    (extract_constraints_boltz1, (8.0, 6.0)),
    (extract_constraints_boltz2, (10.0, 50, 0.0, 8.0, 8.0)),
])
def test_extractors_exit_when_a_chain_has_no_ca(module, args, make_pdb_file, capsys):
    """Wrong chain IDs must fail loudly rather than emit an empty block."""
    single = [{"chain": "B", "resseq": 1, "x": 0.0, "y": 0.0, "z": 0.0}]
    pdb = make_pdb_file(single, name=f"{module.__name__}_nochain.pdb")
    with pytest.raises(SystemExit):
        _extract(module, pdb, *args, capsys=capsys)


@pytest.mark.parametrize("module,args", [
    (extract_constraints_boltz1, (1.0, 6.0)),
    (extract_constraints_boltz2, (1.0, 50, 0.0, 1.0, 8.0)),
])
def test_extractors_exit_when_no_residues_are_close_enough(
        module, args, make_pdb_file, two_chain_atoms, capsys):
    pdb = make_pdb_file(two_chain_atoms, name=f"{module.__name__}_far.pdb")
    with pytest.raises(SystemExit):
        _extract(module, pdb, *args, capsys=capsys)


def test_boltz2_extractor_emits_pocket_and_contacts(
        make_pdb_file, two_chain_atoms, capsys):
    pdb = make_pdb_file(two_chain_atoms, name="c2.pdb")
    out = _extract(extract_constraints_boltz2, pdb, 10.0, 50, 0.0, 8.0, 8.0,
                   capsys=capsys)
    assert "- pocket:" in out
    assert "- contact:" in out


def test_boltz2_extractor_honours_the_contact_cap(
        make_pdb_file, two_chain_atoms, capsys):
    pdb = make_pdb_file(two_chain_atoms, name="c3.pdb")
    out = _extract(extract_constraints_boltz2, pdb, 10.0, 1, 0.0, 8.0, 8.0,
                   capsys=capsys)
    assert out.count("- contact:") <= 1


@pytest.mark.parametrize("module,args", [
    (extract_constraints_boltz1, (8.0,)),
    (extract_constraints_boltz2, (10.0, 50, 0.0, 8.0)),
])
def test_extractors_reject_the_wrong_argument_count(module, args, capsys):
    with pytest.raises(SystemExit):
        _extract(module, "nope.pdb", *args, capsys=capsys)

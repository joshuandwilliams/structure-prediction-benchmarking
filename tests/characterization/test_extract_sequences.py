"""Unit tests for sequence extraction.

Every prediction starts from the sequences this produces, and the shell env
file it writes is interpolated straight into the Nextflow process scripts.
"""

import json
import sys

import extract_sequences as es
import pytest

pytestmark = pytest.mark.local_unit


@pytest.fixture
def two_chain_pdb(make_pdb_file):
    atoms = ([{"chain": "A", "resseq": i + 1, "resname": r,
               "x": float(i), "y": 0.0, "z": 0.0}
              for i, r in enumerate(["MET", "LYS", "PHE", "LEU"])] +
             [{"chain": "B", "resseq": i + 1, "resname": r,
               "x": float(i), "y": 20.0, "z": 0.0}
              for i, r in enumerate(["GLY", "THR", "ALA"])])
    return make_pdb_file(atoms, name="pair.pdb")


def _run(pdb, tmp_path, *extra):
    out, env = tmp_path / "seq.json", tmp_path / "seq.env"
    sys.argv = ["extract_sequences.py", str(pdb),
                "--output", str(out), "--env", str(env), *extra]
    es.main()
    return json.loads(out.read_text()), env.read_text()


def test_sequences_are_read_in_one_letter_code(two_chain_pdb, tmp_path, capsys):
    data, _ = _run(two_chain_pdb, tmp_path)
    by = {c["id"]: c for c in data["chains"]}
    assert by["A"]["sequence"] == "MKFL"
    assert by["B"]["sequence"] == "GTA"
    assert by["A"]["length"] == 4


def test_the_env_file_exposes_both_chains_to_the_shell(two_chain_pdb, tmp_path,
                                                       capsys):
    """Nextflow interpolates these directly into the predictor scripts."""
    _, env = _run(two_chain_pdb, tmp_path)
    assert 'CHAIN_A_SEQ="MKFL"' in env
    assert 'CHAIN_B_SEQ="GTA"' in env
    assert "CHAIN_A_LEN=4" in env
    assert "CHAIN_B_LEN=3" in env


def test_chains_flag_pins_which_chain_is_a(two_chain_pdb, tmp_path, capsys):
    """The reference may name its chains anything, so the caller chooses."""
    _, env = _run(two_chain_pdb, tmp_path, "--chains", "B", "A")
    assert 'CHAIN_A_SEQ="GTA"' in env
    assert 'CHAIN_B_SEQ="MKFL"' in env


def test_chains_can_come_from_the_environment(two_chain_pdb, tmp_path,
                                              monkeypatch, capsys):
    monkeypatch.setenv("BENCHMARK_CHAINS", "B A")
    _, env = _run(two_chain_pdb, tmp_path)
    assert 'CHAIN_A_SEQ="GTA"' in env


def test_the_json_carries_the_shape_main_nf_parses(two_chain_pdb, tmp_path,
                                                   capsys):
    """main.nf reads data.chains with id, sequence and length per entry."""
    data, _ = _run(two_chain_pdb, tmp_path)
    assert set(data) == {"chains", "chain_ids", "num_chains"}
    assert data["chain_ids"] == ["A", "B"]
    assert data["num_chains"] == 2
    assert set(data["chains"][0]) == {"id", "sequence", "length"}


def test_a_file_with_no_protein_exits(tmp_path, make_pdb_file, capsys):
    empty = tmp_path / "empty.pdb"
    empty.write_text("END\n")
    sys.argv = ["extract_sequences.py", str(empty)]
    with pytest.raises(SystemExit):
        es.main()


def test_selenomethionine_maps_to_methionine(make_pdb_file, tmp_path, capsys):
    """Crystal structures substitute MSE for MET, and the predictors expect M."""
    atoms = [{"chain": "A", "resseq": i + 1, "resname": r,
              "x": float(i), "y": 0.0, "z": 0.0}
             for i, r in enumerate(["MSE", "LYS", "PHE"])]
    data, _ = _run(make_pdb_file(atoms, name="mse.pdb"), tmp_path)
    assert data["chains"][0]["sequence"] == "MKF"


def test_output_defaults_to_stdout(two_chain_pdb, capsys):
    sys.argv = ["extract_sequences.py", str(two_chain_pdb)]
    es.main()
    assert "MKFL" in capsys.readouterr().out

"""Unit tests for bin/extract_sequences_from_fasta.py — the FASTA-mode input
adapter. Pure standard library, runs on a laptop.
"""

import json
import subprocess
import sys
from pathlib import Path

import extract_sequences_from_fasta as esf
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "bin" / "extract_sequences_from_fasta.py"

pytestmark = pytest.mark.local_unit


# ── parse_fasta ───────────────────────────────────────────────────────────

def test_parse_two_entries(tmp_path):
    p = tmp_path / "in.fasta"
    p.write_text(">receptor\nMKFL\n>effector\nGTAL\n")
    entries = esf.parse_fasta(str(p))
    assert entries == [("receptor", "MKFL"), ("effector", "GTAL")]


def test_parse_multiline_and_whitespace_and_case(tmp_path):
    p = tmp_path / "in.fasta"
    p.write_text(">a\nmk fl\nGT al\n\n>b\nqrst\n")
    entries = esf.parse_fasta(str(p))
    assert entries == [("a", "MKFLGTAL"), ("b", "QRST")]


def test_parse_strips_trailing_stop_codon(tmp_path):
    p = tmp_path / "in.fasta"
    p.write_text(">a\nMKFL*\n")
    assert esf.parse_fasta(str(p)) == [("a", "MKFL")]


def test_parse_sequence_before_header_raises(tmp_path):
    p = tmp_path / "in.fasta"
    p.write_text("MKFL\n>a\nGTAL\n")
    with pytest.raises(ValueError):
        esf.parse_fasta(str(p))


def test_parse_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        esf.parse_fasta("/no/such/file.fasta")


# ── validate_sequence ─────────────────────────────────────────────────────

def test_validate_accepts_protein():
    esf.validate_sequence("ok", "MKFLGTAL")  # no raise


def test_validate_rejects_empty():
    with pytest.raises(ValueError):
        esf.validate_sequence("empty", "")


def test_validate_rejects_non_aa():
    with pytest.raises(ValueError):
        esf.validate_sequence("dna-ish", "MKFL123")


# ── CLI end-to-end ────────────────────────────────────────────────────────

def _run(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        capture_output=True, text=True,
    )


def test_cli_emits_sequences_json(tmp_path):
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">rec\nMKFLMKFL\n>eff\nGTAL\n")
    out = tmp_path / "sequences.json"
    r = _run(str(fasta), "--output", str(out), "--chains", "A", "B")
    assert r.returncode == 0, r.stderr
    data = json.loads(out.read_text())
    assert data["num_chains"] == 2
    assert data["chain_ids"] == ["A", "B"]
    assert data["chains"][0] == {
        "id": "A", "sequence": "MKFLMKFL", "length": 8, "fasta_header": "rec"
    }
    assert data["chains"][1]["length"] == 4


def test_cli_rejects_wrong_entry_count(tmp_path):
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">only\nMKFL\n")
    r = _run(str(fasta), "--output", str(tmp_path / "x.json"))
    assert r.returncode == 1
    assert "exactly 2 entries" in r.stderr


def test_cli_rejects_identical_chain_ids(tmp_path):
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">rec\nMKFL\n>eff\nGTAL\n")
    r = _run(str(fasta), "--chains", "A", "A")
    assert r.returncode == 1
    assert "must differ" in r.stderr


def test_cli_writes_env_aliases(tmp_path):
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">rec\nMKFL\n>eff\nGTALGTAL\n")
    env = tmp_path / "sequences.env"
    r = _run(str(fasta), "--output", str(tmp_path / "s.json"),
             "--env", str(env), "--chains", "A", "B")
    assert r.returncode == 0, r.stderr
    text = env.read_text()
    assert 'CHAIN_A_SEQ="MKFL"' in text
    assert 'CHAIN_B_SEQ="GTALGTAL"' in text
    assert "CHAIN_B_LEN=8" in text

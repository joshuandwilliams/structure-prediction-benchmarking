"""Unit tests for the Nextflow trace parser.

This replaced ~150 lines of Groovy in main.nf's onComplete hook, which could
not be tested at all. standalone_elapsed_s is the column every cost figure
uses, so the MSA-folding rule is pinned here.
"""

import csv

import pytest
import trace_to_runtime_csv as t

pytestmark = pytest.mark.local_unit

HEADER = ("task_id\tprocess\tstatus\texit\trealtime\t%cpu\trss\tvmem\t"
          "peak_rss\tpeak_vmem\tqueue\n")


def _trace(*rows):
    return [HEADER] + [r if r.endswith("\n") else r + "\n" for r in rows]


def _row(proc, realtime="60000", status="COMPLETED"):
    return (f"1\t{proc}\t{status}\t0\t{realtime}\t100%\t"
            f"1 GB\t2 GB\t1 GB\t2 GB\tjic-gpu")


# ── Duration parsing ──────────────────────────────────────────────────────

@pytest.mark.parametrize("text,expected", [
    ("60000", 60),                 # plain milliseconds, the usual DSL2 form
    ("1000", 1),
    ("45.3s", 45),
    ("2m 10s", 130),
    ("1h 23m 45s", 5025),
    ("500ms", 0),
])
def test_parse_duration_handles_both_emitted_formats(text, expected):
    assert t.parse_duration_s(text) == expected


@pytest.mark.parametrize("text", ["", "-", None, "not a duration"])
def test_parse_duration_returns_none_on_missing(text):
    assert t.parse_duration_s(text) is None


# ── Memory parsing ────────────────────────────────────────────────────────

@pytest.mark.parametrize("text,expected", [
    ("1 GB", 1.0), ("512 MB", 0.5), ("2 TB", 2048.0), ("1048576 KB", 1.0),
])
def test_parse_mem_converts_to_gb(text, expected):
    assert t.parse_mem_gb(text) == pytest.approx(expected)


@pytest.mark.parametrize("text", ["", "-", "0", None, "lots"])
def test_parse_mem_returns_none_on_missing(text):
    assert t.parse_mem_gb(text) is None


# ── Row selection ─────────────────────────────────────────────────────────

def test_only_predictor_processes_are_kept():
    rows = t.rows_from_trace(_trace(
        _row("BOLTZ2"), _row("METRICS_BOLTZ2"), _row("COLABFOLD_SEARCH"),
        _row("EXTRACT_SEQUENCES"), _row("AGGREGATE_RESULTS")))
    assert [r["model"] for r in rows] == ["boltz2"]


def test_every_template_variant_is_recognised():
    rows = t.rows_from_trace(_trace(
        _row("BOLTZ2_TEMPLATE"), _row("BOLTZ2_MSA_TEMPLATE"),
        _row("BOLTZ2_CONSTRAINED_TEMPLATE")))
    assert len(rows) == 3


# ── The MSA-folding rule ──────────────────────────────────────────────────

def test_msa_search_time_is_added_only_to_variants_that_need_it():
    """The property every cost comparison rests on.

    Variants that cannot run without the shared search carry its cost;
    everything else has standalone equal to elapsed.
    """
    rows = t.rows_from_trace(_trace(
        _row("COLABFOLD_SEARCH", "120000"),   # 120 s
        _row("BOLTZ2", "60000"),              # no MSA
        _row("BOLTZ2_MSA", "60000"),
        _row("BOLTZ2_MSA_TEMPLATE", "60000"),
        _row("COLABFOLD", "60000"),
        _row("ESMFOLD2", "60000")))
    by = {r["model"]: r for r in rows}

    assert by["boltz2"]["standalone_elapsed_s"] == 60
    assert by["esmfold2"]["standalone_elapsed_s"] == 60
    for m in ("boltz2_msa", "boltz2_msa_template", "colabfold"):
        assert by[m]["standalone_elapsed_s"] == 180, m


def test_no_msa_search_row_leaves_standalone_equal_to_elapsed():
    rows = t.rows_from_trace(_trace(_row("BOLTZ2_MSA", "60000")))
    assert rows[0]["standalone_elapsed_s"] == rows[0]["elapsed_s"] == 60


def test_longest_msa_search_wins_when_several_ran():
    rows = t.rows_from_trace(_trace(
        _row("COLABFOLD_SEARCH", "60000"), _row("COLABFOLD_SEARCH", "120000"),
        _row("BOLTZ2_MSA", "60000")))
    assert rows[0]["standalone_elapsed_s"] == 180


# ── Output shape ──────────────────────────────────────────────────────────

def test_hms_matches_the_seconds_column():
    rows = t.rows_from_trace(_trace(_row("BOLTZ2", "5025000")))
    assert rows[0]["elapsed_s"] == 5025
    assert rows[0]["elapsed_hms"] == "01:23:45"


def test_missing_runtime_leaves_the_row_blank_rather_than_zero():
    """A zero would be read as a free prediction; blank drops out of medians."""
    rows = t.rows_from_trace(_trace(_row("BOLTZ2", "-")))
    assert rows[0]["elapsed_s"] == ""
    assert rows[0]["standalone_elapsed_s"] == ""


def test_written_csv_has_the_expected_columns(tmp_path, monkeypatch):
    trace = tmp_path / "trace.txt"
    trace.write_text("".join(_trace(_row("BOLTZ2"))))
    out = tmp_path / "rt.csv"
    monkeypatch.setattr("sys.argv", ["x", str(trace), str(out)])
    t.main()
    with open(out) as fh:
        got = list(csv.DictReader(fh))
    assert list(got[0]) == t.FIELDS
    assert got[0]["model"] == "boltz2"

"""Error and edge-case branches across bin/.

These are the paths that fire on malformed input rather than in a normal run.
They matter because a benchmark that silently swallows a bad reference, or
crashes halfway through 18 targets, wastes GPU hours before anyone notices.
"""

import sys

import _analysis_common as ac
import _constraint_geometry as geom
import _structure as S
import build_per_prediction_csv as bpp
import build_runtime_csv as brc
import numpy as np
import pytest
import trace_to_runtime_csv as ttr

pytestmark = pytest.mark.local_unit


# ── _structure ────────────────────────────────────────────────────────────

def test_read_chain_raises_for_an_absent_chain(make_pdb_file, two_chain_atoms):
    pdb = make_pdb_file(two_chain_atoms)
    with pytest.raises(KeyError, match="Z"):
        S.read_chain(pdb, "Z")


def test_read_chains_skips_non_amino_acid_residues(make_pdb_file):
    """Waters, ions and ligands are not part of the complex being scored."""
    atoms = [{"chain": "A", "resseq": 1, "resname": "ALA",
              "x": 0.0, "y": 0.0, "z": 0.0},
             {"chain": "A", "resseq": 2, "resname": "HOH",
              "x": 5.0, "y": 0.0, "z": 0.0},
             {"chain": "A", "resseq": 3, "resname": "ZN",
              "x": 9.0, "y": 0.0, "z": 0.0}]
    chains = S.read_chains(make_pdb_file(atoms, name="het.pdb"))
    assert chains["A"][1] == "A"          # the alanine only


def test_read_chains_returns_empty_for_a_file_with_no_protein(tmp_path):
    p = tmp_path / "empty.pdb"
    p.write_text("END\n")
    assert S.read_chains(p) == {}


def test_read_chains_returns_empty_when_gemmi_finds_no_models(tmp_path):
    """gemmi yields a structure with zero models for a file it cannot make
    sense of, rather than raising, so indexing st[0] threw IndexError and took
    the whole benchmark process down. This is not hypothetical: gemmi 0.6.5 in
    the metrics container does exactly that with the mmCIF the esm library
    writes, which crashed METRICS_ESMFOLD2 on a live run."""
    p = tmp_path / "nomodels.cif"
    p.write_text("data_empty\n#\n_entry.id empty\n#\n")
    assert S.read_chains(p) == {}


def test_kabsch_handles_empty_input():
    rmsd, R, t, n = S.kabsch(np.zeros((0, 3)), np.zeros((0, 3)))
    assert n == 0 and np.isnan(rmsd)


def test_kabsch_truncates_to_the_shorter_array():
    a = np.random.default_rng(0).normal(size=(10, 3))
    rmsd, _, _, n = S.kabsch(a, a[:6])
    assert n == 6


def test_matched_indices_returns_empty_on_empty_sequences():
    assert S.matched_indices("", "ACDE") == ([], [])
    assert S.matched_indices("ACDE", "") == ([], [])


def test_contact_pairs_returns_empty_when_a_side_is_empty():
    assert S.contact_pairs(np.zeros((0, 3)), np.ones((3, 3)), 8.0) == []
    assert S.contact_pairs(np.ones((3, 3)), np.zeros((0, 3)), 8.0) == []


def test_contact_pairs_are_sorted_nearest_first():
    a = np.array([[0.0, 0, 0]])
    b = np.array([[5.0, 0, 0], [1.0, 0, 0], [3.0, 0, 0]])
    dists = [d for _, _, d in S.contact_pairs(a, b, 8.0)]
    assert dists == sorted(dists)


def test_one_letter_maps_unknown_residues_to_x():
    assert S.one_letter("XYZ") == "X"
    assert S.one_letter("HOH") == "X"


# ── _constraint_geometry ──────────────────────────────────────────────────

def test_pocket_and_contacts_are_empty_when_a_chain_is_missing():
    assert geom.pocket_residues({}, {1: (0, 0, 0)}, 8.0) == []
    assert geom.contact_pairs({1: (0, 0, 0)}, {}, 8.0, 50) == []


def test_distance_matches_the_euclidean_norm():
    assert geom.distance((0, 0, 0), (3, 4, 0)) == pytest.approx(5.0)


# ── trace_to_runtime_csv ──────────────────────────────────────────────────

def test_trace_with_only_a_header_exits(tmp_path):
    t = tmp_path / "trace.txt"
    t.write_text("task_id\tprocess\n")
    sys.argv = ["trace_to_runtime_csv.py", str(t), str(tmp_path / "o.csv")]
    with pytest.raises(SystemExit):
        ttr.main()


def test_an_unreadable_trace_exits(tmp_path):
    sys.argv = ["trace_to_runtime_csv.py", str(tmp_path / "nope.txt"),
                str(tmp_path / "o.csv")]
    with pytest.raises(SystemExit):
        ttr.main()


def test_a_trace_with_no_predictor_rows_exits(tmp_path):
    t = tmp_path / "trace.txt"
    t.write_text("process\tstatus\texit\trealtime\t%cpu\trss\tvmem\t"
                 "peak_rss\tpeak_vmem\tqueue\n"
                 "METRICS_BOLTZ2\tCOMPLETED\t0\t1000\t10%\t1 GB\t1 GB\t"
                 "1 GB\t1 GB\tjic-medium\n")
    sys.argv = ["trace_to_runtime_csv.py", str(t), str(tmp_path / "o.csv")]
    with pytest.raises(SystemExit):
        ttr.main()


def test_wrong_argument_count_exits():
    sys.argv = ["trace_to_runtime_csv.py"]
    with pytest.raises(SystemExit):
        ttr.main()


# ── The aggregation scripts ───────────────────────────────────────────────

def test_per_prediction_exits_on_a_missing_benchmarks_dir(tmp_path):
    sys.argv = ["build_per_prediction_csv.py", "--benchmarks-dir",
                str(tmp_path / "nope"), "--output", str(tmp_path / "o.csv")]
    with pytest.raises(SystemExit):
        bpp.main()


def test_runtime_csv_exits_on_a_missing_benchmarks_dir(tmp_path):
    sys.argv = ["build_runtime_csv.py", "--benchmarks-dir",
                str(tmp_path / "nope"), "--output", str(tmp_path / "o.csv")]
    with pytest.raises(SystemExit):
        brc.main()


def test_runtime_csv_reports_targets_with_no_stats_file(tmp_path, capsys):
    (tmp_path / "6G10").mkdir()
    sys.argv = ["build_runtime_csv.py", "--benchmarks-dir", str(tmp_path),
                "--output", str(tmp_path / "o.csv")]
    with pytest.raises(SystemExit):
        brc.main()


# ── _analysis_common ──────────────────────────────────────────────────────

def test_readable_log_replaces_the_decade_ticks():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ac.readable_log(ax.yaxis, [1, 2, 5, 10])
    labels = [t.get_text() for t in ax.get_yticklabels()]
    plt.close(fig)
    assert labels == ["1", "2", "5", "10"]


def test_msa_legend_labels_both_blocks():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    leg = ac.msa_legend(ax)
    texts = [t.get_text() for t in leg.get_texts()]
    plt.close(fig)
    assert texts == ["No MSA", "MSA"]


def test_read_chain_returns_coordinates_and_sequence(make_pdb_file, two_chain_atoms):
    coords, seq = S.read_chain(make_pdb_file(two_chain_atoms), "B")
    assert coords.shape == (3, 3)
    assert len(seq) == 3


def test_a_residue_without_a_ca_is_skipped(make_pdb_file):
    atoms = [{"chain": "A", "resseq": 1, "resname": "ALA", "name": "CB",
              "x": 0.0, "y": 0.0, "z": 0.0},
             {"chain": "A", "resseq": 2, "resname": "LYS",
              "x": 3.0, "y": 0.0, "z": 0.0}]
    chains = S.read_chains(make_pdb_file(atoms, name="noca.pdb"))
    assert chains["A"][1] == "K"


def test_only_the_first_altloc_of_a_residue_is_kept(make_pdb_file):
    """Counting a residue twice would corrupt both the sequence and the
    coordinate array."""
    atoms = [{"chain": "A", "resseq": 1, "resname": "ALA",
              "x": 0.0, "y": 0.0, "z": 0.0},
             {"chain": "A", "resseq": 1, "resname": "SER",
              "x": 0.5, "y": 0.0, "z": 0.0},
             {"chain": "A", "resseq": 2, "resname": "LYS",
              "x": 3.0, "y": 0.0, "z": 0.0}]
    chains = S.read_chains(make_pdb_file(atoms, name="alt.pdb"))
    assert len(chains["A"][1]) == 2


def test_matched_indices_returns_empty_when_alignment_fails(monkeypatch):
    """A failure inside Biopython must degrade to no pairing rather than
    propagate and abort a whole benchmark."""
    class Boom:
        def align(self, a, b):
            raise RuntimeError("aligner exploded")

    monkeypatch.setattr(S, "aligner", lambda: Boom())
    assert S.matched_indices("ACDE", "ACDE") == ([], [])

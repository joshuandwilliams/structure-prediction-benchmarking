"""Build synthetic predictions with known interface labels.

The classifier decides which of two receptor surfaces an effector sits on. To
test that decision we need poses whose correct answer is known in advance,
which the published best_models/ trees cannot give us. Every pose here is
derived from a committed reference structure by a rigid transform, so the
expected label follows from the construction rather than from a benchmark run.

Three pose classes are produced per target.

    correct   the target's own reference complex, unmodified
    pia       the frame effector superposed onto this receptor, which puts it
              on the AVR-Pia surface
    other     the true effector pushed far along the receptor surface normal,
              contacting neither patch

Both scripts enumerate MODEL_MAP and look for <PDB>__<predictor tag>.pdb, so
each pose class borrows a real tag rather than a descriptive name. The tag
carries no meaning here beyond giving the classifier something to iterate, and
POSE_TAGS is the mapping back to what each file actually contains.

Output is the flat layout both scripts accept via --pred-dir. Nothing here
touches the network or a GPU.
"""

import os
import sys

import gemmi
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO, "bin"))

import _structure  # noqa: E402

FRAME = "6Q76"
REFS = os.path.join(REPO, "data", "complexes_for_benchmarking")

# Kept small on purpose. The classifier is per-target, so three targets
# exercise every branch, and a 14-target fixture set would be committed bulk
# for no extra coverage.
TARGETS = ["6G10", "6FU9", FRAME]

# Pose class -> predictor tag it is written under. boltz2 carries the correct
# pose because common_wrong_interface defaults to --model boltz2 --msa no_msa
# and resolves that to the "boltz2" tag.
POSE_TAGS = {"correct": "boltz2", "pia": "boltz1", "other": "chai1"}

OFFSET_A = 40.0     # far enough that no residue pair falls inside any cutoff


def _read(path):
    st = gemmi.read_structure(path)
    st.setup_entities()
    return st


def _chain_ca(st, chain_id):
    out = []
    for res in st[0][chain_id]:
        at = res.find_atom("CA", "*")
        if at is not None:
            out.append(np.array([at.pos.x, at.pos.y, at.pos.z]))
    return np.array(out)


def _shift_chain(st, chain_id, vec):
    for res in st[0][chain_id]:
        for at in res:
            at.pos = gemmi.Position(
                at.pos.x + vec[0], at.pos.y + vec[1], at.pos.z + vec[2])


def _transform_chain(st, chain_id, R, t):
    """Apply a kabsch transform using the same convention as _structure.

    apply_transform is `coords @ R + t`, the row-vector form. Writing it as
    `R @ v + t` here silently applies the transpose, which superposes onto a
    mirrored frame and throws the effector tens of angstroms off.
    """
    for res in st[0][chain_id]:
        for at in res:
            v = np.array([at.pos.x, at.pos.y, at.pos.z])
            w = v @ R + t
            at.pos = gemmi.Position(*w)


def _write(st, out_path):
    st.setup_entities()
    with open(out_path, "w") as fh:
        fh.write(st.make_pdb_string())


def build(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    frame_ref = os.path.join(REFS, f"{FRAME}.pdb")
    frame_rec, frame_rec_seq = _structure.read_chain(frame_ref, "A")

    written = []
    for pdb in TARGETS:
        ref = os.path.join(REFS, f"{pdb}.pdb")
        if not os.path.isfile(ref):
            sys.exit(f"ERROR: missing reference {ref}")

        rec, rec_seq = _structure.read_chain(ref, "A")

        # correct: the reference itself is by definition on the true surface.
        _write(_read(ref), os.path.join(out_dir, f"{pdb}__{POSE_TAGS['correct']}.pdb"))
        written.append(f"{pdb}__{POSE_TAGS['correct']}.pdb")

        # pia: superpose the frame receptor onto this receptor, apply that same
        # transform to the frame effector, and graft it in as chain B. The
        # effector then occupies the AVR-Pia surface of THIS receptor.
        ia, ib = _structure.matched_indices(frame_rec_seq, rec_seq)
        if len(ia) < 20:
            sys.exit(f"ERROR: {pdb} aligns to the frame at only {len(ia)} residues")
        _, R, t, _ = _structure.kabsch(frame_rec[list(ia)], rec[list(ib)])

        st = _read(ref)
        frame_st = _read(frame_ref)
        _transform_chain(frame_st, "B", R, t)
        del st[0]["B"]
        st[0].add_chain(frame_st[0]["B"].clone())
        _write(st, os.path.join(out_dir, f"{pdb}__{POSE_TAGS['pia']}.pdb"))
        written.append(f"{pdb}__{POSE_TAGS['pia']}.pdb")

        # other: displace the true effector clear of both surfaces.
        st = _read(ref)
        away = rec.mean(0) - _chain_ca(st, "B").mean(0)
        norm = np.linalg.norm(away)
        away = (away / norm) if norm > 1e-6 else np.array([1.0, 0.0, 0.0])
        _shift_chain(st, "B", -away * OFFSET_A)
        _write(st, os.path.join(out_dir, f"{pdb}__{POSE_TAGS['other']}.pdb"))
        written.append(f"{pdb}__{POSE_TAGS['other']}.pdb")

    return written


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "interface_predictions"
    for name in build(target):
        print(f"  wrote {name}")

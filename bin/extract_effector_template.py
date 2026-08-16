#!/usr/bin/env python3
"""Extract the effector chain as a structural template for AF3 and ColabFold.

AF3's template parser requires _entity_poly_seq, the one-letter sequence code,
integer label_seq_id and a revision-history date, none of which
gemmi.make_mmcif_document populates from a PDB source. This prunes to the
effector chain, relabels it to A, then fills those fields explicitly and
validates them before writing.

Usage:
    extract_effector_template.py <ref_pdb> <effector_chain>
"""

import sys

import gemmi


def extract_effector_template(ref_pdb, effector_chain_id):
    st = gemmi.read_structure(ref_pdb)

    chain_names = {c.name for model in st for c in model}
    if effector_chain_id not in chain_names:
        print(f"ERROR: chain '{effector_chain_id}' not found in {ref_pdb}. "
              f"Available: {sorted(chain_names)}", file=sys.stderr)
        sys.exit(1)

    # Prune every chain except the effector across all models
    for model in st:
        for cname in [c.name for c in model if c.name != effector_chain_id]:
            model.remove_chain(cname)

    # Remove waters — mixing polymer + non-polymer rows in a single entity
    # confuses AF3's template featuriser
    st.remove_waters()

    # Relabel surviving chain to 'A' (single-chain template convention
    # expected by ColabFold --custom-template-path)
    if effector_chain_id != "A":
        for model in st:
            for chain in model:
                if chain.name == effector_chain_id:
                    chain.name = "A"

    st.name = "effector_template"

    # Build polymer Entity (required before populating sequence metadata)
    st.setup_entities()

    # Explicitly populate Entity.full_sequence and Residue.label_seq.
    # gemmi reads PDB files without these fields; without them the CIF
    # writer emits '.' placeholders and AF3 fails with:
    #   ValueError: invalid literal for int() with base 10: '.'
    for model in st:
        for chain in model:
            polymer = chain.get_polymer()
            entity = st.get_entity_of(polymer)
            if entity is None:   # pragma: no cover
                print(f"ERROR: No Entity for chain '{chain.name}' after "
                      "setup_entities() — input PDB may be malformed",
                      file=sys.stderr)
                sys.exit(1)
            entity.full_sequence = [r.name for r in polymer]
            for i, r in enumerate(polymer, start=1):
                r.label_seq = i

    # AF3 requires _pdbx_audit_revision_history.revision_date; structures
    # derived from PDB coordinates do not carry it.  A fixed historic date
    # is fine — AF3 is not date-filtering when templates are injected directly.
    st.info["_pdbx_database_status.recvd_initial_deposition_date"] = "2024-01-01"

    # Write mmCIF with MmcifOutputGroups(True) to get all fields
    groups = gemmi.MmcifOutputGroups(True)
    doc = st.make_mmcif_document(groups)
    block = doc[0]
    loop = block.init_mmcif_loop("_pdbx_audit_revision_history.", [
        "ordinal", "data_content_type", "major_revision",
        "minor_revision", "revision_date",
    ])
    loop.add_row(["1", "'Structure model'", "1", "0", "2024-01-01"])
    doc.write_file("effector_template.cif")

    # Validate — fail fast here rather than 30+ minutes into a GPU job
    text = open("effector_template.cif").read()

    if "_entity_poly_seq.entity_id" not in text:   # pragma: no cover - gemmi writes this
        print("ERROR: _entity_poly_seq missing from generated CIF — "
              "Entity.full_sequence was not populated", file=sys.stderr)
        sys.exit(1)

    if "_pdbx_audit_revision_history.revision_date" not in text:   # pragma: no cover - gemmi writes this
        print("ERROR: _pdbx_audit_revision_history.revision_date missing — "
              "AF3 will fail with 'The structure must have a release date'",
              file=sys.stderr)
        sys.exit(1)

    check_doc = gemmi.cif.read("effector_template.cif")
    check_vals = list(check_doc[0].find_values(
        "_pdbx_audit_revision_history.revision_date"
    ))
    if not check_vals:   # pragma: no cover
        print("ERROR: revision_date in text but not discoverable by mmCIF "
              "parser — AF3 will fail", file=sys.stderr)
        sys.exit(1)

    for line in text.splitlines():
        if line.startswith(("ATOM ", "HETATM ")):
            cols = line.split()
            if len(cols) >= 9 and cols[8] == ".":   # pragma: no cover
                print("ERROR: _atom_site.label_seq_id is '.' — "
                      "Residue.label_seq was not populated", file=sys.stderr)
                sys.exit(1)

    # Count polymer residues for logging
    n_polymer = sum(1 for chain in st[0] for res in chain
                    if res.entity_type == gemmi.EntityType.Polymer)
    print(f"Extracted chain {effector_chain_id} -> chain A: "
          f"{n_polymer} polymer residues — mmCIF validated OK")
    print("Written: effector_template.cif")

    # Write PDB for ColabFold --custom-template-path
    st.write_pdb("effector_template.pdb")
    print("Written: effector_template.pdb")


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <ref_pdb> <effector_chain>",
              file=sys.stderr)
        sys.exit(1)
    extract_effector_template(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()

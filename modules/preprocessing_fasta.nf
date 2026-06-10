/*
 * =============================================================================
 * Preprocessing module — FASTA input variant
 * =============================================================================
 * Alternative to preprocessing.nf::EXTRACT_SEQUENCES for runs where the user
 * supplies a two-entry FASTA file instead of a reference PDB.  Emits an
 * identically-shaped sequences.json so every downstream module is unchanged.
 *
 * Because there is no experimental reference structure in FASTA mode, this
 * process also emits a placeholder reference.pdb file containing only a
 * comment line.  compute_metrics.py treats a missing/empty reference as a
 * warning and skips the RMSD/DockQ block, while still producing all
 * confidence-based metrics (pLDDT, ipTM, pTM, etc.).  This means the
 * existing metrics.nf and aggregate.nf modules need NO modifications.
 *
 * The two FASTA entries are interpreted in file order:
 *   entry 1 -> receptor (assigned chain ID receptor_chain, default 'A')
 *   entry 2 -> effector (assigned chain ID effector_chain, default 'B')
 */


/*
 * EXTRACT_SEQUENCES_FROM_FASTA
 * ----------------------------
 * Parse a two-entry FASTA file via extract_sequences_from_fasta.py and
 * emit sequences.json + a placeholder reference.pdb.  Output channels
 * mirror EXTRACT_SEQUENCES so main.nf can choose between the two
 * preprocessing paths transparently.
 */
process EXTRACT_SEQUENCES_FROM_FASTA {
    tag "${params.project_name}"
    label 'cpu'

    publishDir "${params.outdir}", mode: 'copy',
        pattern: '{sequences.json,reference.pdb,input.fasta}'

    input:
    path input_fasta
    val  receptor_chain
    val  effector_chain

    output:
    path "sequences.json", emit: sequences_json
    path "reference.pdb",  emit: reference_pdb
    path "input.fasta",    emit: input_fasta_copy

    script:
    """
    set -euo pipefail

    # Keep an unambiguous copy of the input FASTA alongside sequences.json.
    cp ${input_fasta} input.fasta

    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/extract_sequences_from_fasta.py \\
            input.fasta \\
            --output sequences.json \\
            --env sequences.env \\
            --chains ${receptor_chain} ${effector_chain}

    # Placeholder reference.pdb so the metrics process can stage a file with
    # this name.  compute_metrics.py checks os.path.exists() on the path it
    # is given and, when no real reference structure is supplied, falls back
    # to confidence-only metrics.  We give it a single REMARK line so the
    # file is non-empty but contains no ATOM records — the parser will read
    # zero Cα atoms and skip the structural-RMSD block with a warning.
    cat > reference.pdb << 'PDBEOF'
REMARK   1 PLACEHOLDER REFERENCE — FASTA INPUT MODE
REMARK   1 NO EXPERIMENTAL STRUCTURE WAS PROVIDED.  STRUCTURAL RMSD AND
REMARK   1 DOCKQ METRICS WILL BE SKIPPED BY compute_metrics.py.  CONFIDENCE
REMARK   1 METRICS (pLDDT / ipTM / pTM / etc.) ARE STILL COMPUTED.
END
PDBEOF
    """
}

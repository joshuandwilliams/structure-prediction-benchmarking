/*
 * =============================================================================
 * ColabFold module
 * =============================================================================
 * Two variants sharing colabfold.img:
 *
 *   COLABFOLD        — uses the shared COLABFOLD_SEARCH a3m MSAs as input
 *                      to colabfold_batch.  Because colabfold_batch accepts
 *                      a directory of a3m files as input, we stage the
 *                      shared msa/ dir from COLABFOLD_SEARCH directly.
 *   COLABFOLD_NOMSA  — runs colabfold_batch with --msa-mode single_sequence,
 *                      skipping all MSA search entirely.  Takes a FASTA
 *                      directly (no upstream COLABFOLD_SEARCH).
 *
 * Both use alphafold2_multimer_v3 weights, 5 models × 5 seeds, 20 recycles.
 */


/*
 * COLABFOLD
 * ---------
 * ColabFold prediction using the shared MMseqs2 MSAs produced by
 * COLABFOLD_SEARCH.  colabfold_batch reads a directory of a3m files.
 */
process COLABFOLD {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/colabfold", mode: 'copy',
        pattern: 'output/**'

    input:
    path msa_dir    // staged chain_A.a3m / chain_B.a3m directory

    output:
    path "output",  emit: prediction_dir

    script:
    """
    set -euo pipefail

    mkdir -p output

    echo "=== Structure Prediction (colabfold_batch) ==="
    singularity exec --nv ${params.colabfold_container} colabfold_batch \\
        ${msa_dir} \\
        output \\
        --model-type alphafold2_multimer_v3 \\
        --num-models 5 \\
        --num-recycle 20 \\
        --num-seeds 5 \\
        --rank iptm

    echo "ColabFold output files:"
    find output -type f | head -20 || true
    """
}


/*
 * COLABFOLD_NOMSA
 * ---------------
 * ColabFold in single-sequence mode via --msa-mode single_sequence.  No
 * MSA search is performed; input is a joined-style FASTA (chain_A:chain_B
 * in one entry) that colabfold_batch accepts directly.
 */
process COLABFOLD_NOMSA {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/colabfold_nomsa", mode: 'copy',
        pattern: '{input.fasta,output/**}'

    input:
    val receptor_seq
    val effector_seq

    output:
    path "output",      emit: prediction_dir
    path "input.fasta"

    script:
    """
    set -euo pipefail

    # ─── Joined-format FASTA (chain_A:chain_B) ────────────────────────────
    # colabfold_batch accepts this directly when --msa-mode single_sequence.
    cat > input.fasta << EOF
>complex
${receptor_seq}:${effector_seq}
EOF

    mkdir -p output

    echo "=== Structure Prediction (colabfold_batch --msa-mode single_sequence) ==="
    singularity exec --nv ${params.colabfold_container} colabfold_batch \\
        input.fasta \\
        output \\
        --model-type alphafold2_multimer_v3 \\
        --msa-mode single_sequence \\
        --num-models 5 \\
        --num-recycle 20 \\
        --num-seeds 5 \\
        --rank iptm

    echo "ColabFold (no MSA) output files:"
    find output -type f | head -20 || true
    """
}

/*
 * =============================================================================
 * AlphaFold2-Multimer module
 * =============================================================================
 * Uses the NBI HPC AlphaFold 2.3.2 installation via `source package`.
 *
 * Unlike the original bench_af2m.sh, this process does NOT need to snapshot
 * and restore the pre-source-package environment before invoking the metrics
 * container — COMPUTE_METRICS runs in its own separate Nextflow process with
 * its own clean environment, so PATH/LD_LIBRARY_PATH/PYTHONPATH pollution
 * stays contained inside this work directory.
 *
 * Resource preset matches the original SBATCH header: 20 CPU, 64 GB RAM,
 * 1 GPU, 24h walltime.
 */


/*
 * AF2M
 * ----
 * Run AlphaFold2-Multimer 2.3.2 with full MSA search (reduced_dbs preset)
 * and real template search capped at 2020-05-14 for consistent benchmarking.
 * Produces 5 models × 5 multimer predictions each = 25 structures by default.
 */
process AF2M {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/af2m", mode: 'copy',
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

    export TF_FORCE_GPU_ALLOW_GROWTH=true
    export TF_FORCE_UNIFIED_MEMORY=1
    export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

    # ─── Write multimer FASTA (separate-style) ────────────────────────────
    cat > input.fasta << EOF
>chain_A
${receptor_seq}
>chain_B
${effector_seq}
EOF

    DATA_DIR="${params.af2_data_dir}"
    mkdir -p output

    # ─── Load AlphaFold 2.3.2 from HPC source package ────────────────────
    source package ${params.af2_package_id}

    echo "Running AF2-Multimer prediction (with MSA)..."

    run_alphafold.py \\
        --fasta_paths=\${PWD}/input.fasta \\
        --data_dir="\${DATA_DIR}" \\
        --output_dir=\${PWD}/output \\
        --max_template_date=2020-05-14 \\
        --model_preset=multimer \\
        --db_preset=reduced_dbs \\
        --small_bfd_database_path="\${DATA_DIR}/small_bfd/bfd-first_non_consensus_sequences.fasta" \\
        --uniref90_database_path="\${DATA_DIR}/uniref90/uniref90.fasta" \\
        --mgnify_database_path="\${DATA_DIR}/mgnify/mgy_clusters_2022_05.fa" \\
        --pdb_seqres_database_path="\${DATA_DIR}/pdb_seqres/pdb_seqres.txt" \\
        --uniprot_database_path="\${DATA_DIR}/uniprot/uniprot.fasta" \\
        --template_mmcif_dir="\${DATA_DIR}/pdb_mmcif/mmcif_files" \\
        --obsolete_pdbs_path="\${DATA_DIR}/pdb_mmcif/obsolete.dat" \\
        --num_multimer_predictions_per_model=5 \\
        --use_gpu_relax

    echo "AF2-Multimer output files:"
    find output -type f | head -20 || true
    """
}

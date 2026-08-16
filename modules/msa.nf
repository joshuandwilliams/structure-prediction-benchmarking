/*
 * Shared MMseqs2 search, run once per target and reused by every variant that
 * needs an MSA (boltz1_msa, boltz2_msa, boltz2_msa_template, colabfold).
 *
 * The input FASTA carries three records so one search yields both forms: the
 * per-chain unpaired a3m that Boltz consumes, and the ':'-joined complex that
 * colabfold_search returns as a single paired a3m for ColabFold.
 */

/*
 * Run MMseqs2 against colabfold_databases to produce paired a3m files for
 * chain_A and chain_B.  The input FASTA is written in "separate" style so
 * colabfold_search names outputs chain_A.a3m / chain_B.a3m directly.
 */
process COLABFOLD_SEARCH {
    tag "${params.project_name}"
    label 'cpu'

    publishDir "${params.outdir}/msa", mode: 'copy'

    input:
    val  receptor_seq
    val  effector_seq

    output:
    path "msa",                emit: msa_dir
    path "msa/chain_A.a3m",    emit: chain_a_a3m
    path "msa/chain_B.a3m",    emit: chain_b_a3m
    path "msa/complex.a3m",    emit: complex_a3m

    script:
    """
    # colabfold_search produces one a3m per FASTA record, named after the
    # header. We request THREE:
    #   chain_A / chain_B  — per-chain (unpaired) MSAs; Boltz consumes these.
    #   complex            — a ':'-joined query, which colabfold_search treats
    #                        as a COMPLEX and returns as a single paired a3m.
    # ColabFold needs that paired complex a3m: handing it the per-chain dir
    # instead makes colabfold_batch fold each chain as a separate monomer.
    cat > input.fasta << 'FASTA_EOF'
>chain_A
${receptor_seq}
>chain_B
${effector_seq}
>complex
${receptor_seq}:${effector_seq}
FASTA_EOF

    if [ ! -d "${params.colabfold_db}" ]; then
        echo "ERROR: ColabFold databases not found at ${params.colabfold_db}" >&2
        exit 1
    fi

    mkdir -p msa

    singularity exec ${params.colabfold_container} \\
        env MMSEQS_IGNORE_INDEX=1 colabfold_search \\
            --mmseqs mmseqs \\
            --use-env 1 \\
            --use-templates 0 \\
            --threads ${task.cpus} \\
            --prefilter-mode 1 \\
            input.fasta \\
            ${params.colabfold_db} \\
            msa

    if [ ! -f msa/chain_A.a3m ] || [ ! -f msa/chain_B.a3m ] || [ ! -f msa/complex.a3m ]; then
        echo "ERROR: Expected MSA files (chain_A.a3m, chain_B.a3m, complex.a3m) not found in msa/" >&2
        ls -la msa/ >&2
        exit 1
    fi

    echo "MSA files:"
    wc -l msa/chain_A.a3m msa/chain_B.a3m msa/complex.a3m
    """
}

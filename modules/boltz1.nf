/*
 * Boltz-1: single-sequence baseline, with a ColabFold MSA, and with a
 * pocket-only constraint. Boltz-1's schema forces max_distance to 6.0 and
 * rejects contact constraints, so the constrained variant is pocket only.
 *
 * --no_kernels is required because cuequivariance_torch is absent from the
 * image. Seeds run serially inside the process to respect the GPU queue.
 */

/*
 * Unconstrained Boltz-1 baseline: single-sequence MSA (just the query), no
 * constraints, 5 seeds × 5 diffusion samples = 25 structures.
 */
process BOLTZ1 {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/boltz1", mode: 'copy',
        pattern: '{input.fasta,all_outputs/**}'

    input:
    val receptor_seq
    val effector_seq
    val rec_chain
    val eff_chain

    output:
    path "all_outputs",  emit: prediction_dir
    path "input.fasta"

    script:
    """
    set -euo pipefail

    mkdir -p msa
    cat > msa/A.a3m << EOF
>query
${receptor_seq}
EOF
    cat > msa/B.a3m << EOF
>query
${effector_seq}
EOF

    cat > input.fasta << EOF
>A|protein|\${PWD}/msa/A.a3m
${receptor_seq}
>B|protein|\${PWD}/msa/B.a3m
${effector_seq}
EOF

    echo "Input FASTA:"
    cat input.fasta
    echo ""

    bash ${projectDir}/bin/run_boltz_seeds.sh \\
        input.fasta ${params.benchmark_container} ${projectDir} "Boltz-1" \\
        --model boltz1
    """
}

/*
 * Boltz-1 with a real MSA.  Receives chain_A.a3m and chain_B.a3m from the
 * shared COLABFOLD_SEARCH process and uses them directly in the Boltz YAML
 * input.
 */
process BOLTZ1_MSA {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/boltz1_msa", mode: 'copy',
        pattern: '{input.yaml,all_outputs/**}'

    input:
    val  receptor_seq
    val  effector_seq
    val  rec_chain
    val  eff_chain
    path chain_a_a3m
    path chain_b_a3m

    output:
    path "all_outputs",  emit: prediction_dir
    path "input.yaml"

    script:
    """
    set -euo pipefail

    # chain_a_a3m/chain_b_a3m are already staged by Nextflow in the work dir.
    # We reference them by their staged names directly in the YAML.

    cat > input.yaml << YAMLEOF
version: 1
sequences:
  - protein:
      id: ${rec_chain}
      sequence: ${receptor_seq}
      msa: \${PWD}/${chain_a_a3m.name}
  - protein:
      id: ${eff_chain}
      sequence: ${effector_seq}
      msa: \${PWD}/${chain_b_a3m.name}
YAMLEOF

    echo "Input YAML:"
    cat input.yaml
    echo ""

    bash ${projectDir}/bin/run_boltz_seeds.sh \\
        input.yaml ${params.benchmark_container} ${projectDir} "Boltz-1 MSA" \\
        --model boltz1
    """
}

/*
 * Boltz-1 with a pocket constraint derived from the reference complex. NOTE:
 * Boltz-1 schema.py enforces max_distance == 6.0 for ALL constraint blocks
 * and silently rejects contact constraints, so only a pocket block is
 * generated. bin/extract_constraints_boltz1.py enforces the exactly-6.0
 * invariant, and bin/validate_boltz_yaml.py --model boltz1 double-checks it
 * before boltz predict runs.
 */
process BOLTZ1_CONSTRAINED {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/boltz1_constrained", mode: 'copy',
        pattern: '{input.yaml,constraints.yaml,all_outputs/**}'

    input:
    val  receptor_seq
    val  effector_seq
    val  rec_chain
    val  eff_chain
    path reference_pdb

    output:
    path "all_outputs",      emit: prediction_dir
    path "input.yaml"
    path "constraints.yaml"

    script:
    // Fixed values (Boltz-1 schema requirement)
    def pocket_cutoff       = 8.0
    def pocket_max_distance = 6.0  // MUST be exactly 6.0
    """
    set -euo pipefail

    echo "=== Extracting pocket constraint from reference PDB ==="
    echo "  Reference PDB:   ${reference_pdb}"
    echo "  Pocket cutoff:   ${pocket_cutoff} Å"
    echo "  Pocket max_dist: ${pocket_max_distance} Å (must be exactly 6.0 for Boltz-1)"
    echo "  Contact constraints: OMITTED (not supported by Boltz-1)"
    echo ""

    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/extract_constraints_boltz1.py \\
            ${reference_pdb} ${rec_chain} ${eff_chain} \\
            ${pocket_cutoff} ${pocket_max_distance} \\
        > constraints.yaml

    N_POCKET=\$(grep -c "^  - pocket:" constraints.yaml || true)
    echo "Generated constraints: \${N_POCKET} pocket, 0 contact (not supported by Boltz-1)"

    if [ "\${N_POCKET}" -eq 0 ]; then
        echo "ERROR: No pocket constraint generated. Check chain IDs in ${reference_pdb}" >&2
        exit 1
    fi

    mkdir -p msa
    cat > msa/A.a3m << EOF
>query
${receptor_seq}
EOF
    cat > msa/B.a3m << EOF
>query
${effector_seq}
EOF

    cat > input.yaml << YAMLEOF
version: 1
sequences:
  - protein:
      id: ${rec_chain}
      sequence: ${receptor_seq}
      msa: \${PWD}/msa/A.a3m
  - protein:
      id: ${eff_chain}
      sequence: ${effector_seq}
      msa: \${PWD}/msa/B.a3m
YAMLEOF

    cat constraints.yaml >> input.yaml

    echo "Full input.yaml:"
    cat input.yaml
    echo ""

    # Validate YAML (max_distance==6.0, no contact blocks)
    echo "=== Validating YAML ==="
    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/validate_boltz_yaml.py input.yaml --model boltz1

    bash ${projectDir}/bin/run_boltz_seeds.sh \\
        input.yaml ${params.benchmark_container} ${projectDir} "Boltz-1 constrained" \\
        --model boltz1
    """
}

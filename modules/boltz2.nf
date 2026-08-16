/*
 * Boltz-2: single-sequence baseline, with a ColabFold MSA, and with pocket
 * plus dense contact constraints, each in a template-free and an
 * effector-template form.
 *
 * The template variants supply the reference effector chain as a structural
 * template (force: true, threshold 1.0, enforced by --use_potentials), matching
 * what structure-negative-steering gives Boltz-2. They are not part of the
 * model comparison. Each pairs with its template-free twin to measure what the
 * known effector fold is worth.
 *
 * --no_kernels is required because cuequivariance_torch is absent from the
 * image. Seeds run serially inside the process to respect the GPU queue.
 */

/*
 * Unconstrained Boltz-2 baseline: single-sequence MSA, no constraints, no
 * structural template.  Nothing derived from the reference complex reaches
 * the model, so this variant is directly comparable with every other
 * predictor.
 */
process BOLTZ2 {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/boltz2", mode: 'copy',
        pattern: '{input.yaml,all_outputs/**}'

    input:
    val receptor_seq
    val effector_seq
    val rec_chain
    val eff_chain

    output:
    path "all_outputs",  emit: prediction_dir
    path "input.yaml"

    script:
    """
    set -euo pipefail

    mkdir -p msa
    cat > msa/chain_A.a3m << EOF
>query
${receptor_seq}
EOF
    cat > msa/chain_B.a3m << EOF
>query
${effector_seq}
EOF

    cat > input.yaml << YAMLEOF
version: 1
sequences:
  - protein:
      id: ${rec_chain}
      sequence: ${receptor_seq}
      msa: \${PWD}/msa/chain_A.a3m
  - protein:
      id: ${eff_chain}
      sequence: ${effector_seq}
      msa: \${PWD}/msa/chain_B.a3m
YAMLEOF

    echo "Input YAML:"
    cat input.yaml
    echo ""

    bash ${projectDir}/bin/run_boltz_seeds.sh \\
        input.yaml ${params.benchmark_container} ${projectDir} "Boltz-2"
    """
}

/*
 * Boltz-2 with the colabfold_search a3m MSAs shared from the
 * COLABFOLD_SEARCH process.  Identical to BOLTZ2 except for the MSA source.
 */
process BOLTZ2_MSA {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/boltz2_msa", mode: 'copy',
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
        input.yaml ${params.benchmark_container} ${projectDir} "Boltz-2 MSA"
    """
}

/*
 * Boltz-2 with pocket + dense contact constraints derived from the reference
 * complex via bin/extract_constraints_boltz2.py.  No template — constraints
 * only — so the reference structure is never directly given to the model
 * (the reference IS the answer we are scoring against).
 */
process BOLTZ2_CONSTRAINED {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/boltz2_constrained", mode: 'copy',
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
    def contact_cutoff      = 10.0
    def contact_max         = 50
    def contact_tolerance   = 0.0
    def pocket_cutoff       = 8.0
    def pocket_max_distance = 8.0
    """
    set -euo pipefail

    echo "=== Extracting constraints from reference PDB ==="
    echo "  Reference PDB:     ${reference_pdb}"
    echo "  Contact cutoff:    ${contact_cutoff} Å"
    echo "  Contact max:       ${contact_max}"
    echo "  Contact tolerance: ${contact_tolerance} Å"
    echo "  Pocket cutoff:     ${pocket_cutoff} Å"
    echo "  Pocket max_dist:   ${pocket_max_distance} Å"
    echo ""

    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/extract_constraints_boltz2.py \\
            ${reference_pdb} ${rec_chain} ${eff_chain} \\
            ${contact_cutoff} ${contact_max} ${contact_tolerance} \\
            ${pocket_cutoff} ${pocket_max_distance} \\
        > constraints.yaml

    N_CONTACTS=\$(grep -c "^  - contact:" constraints.yaml || true)
    N_POCKET=\$(grep -c "^  - pocket:" constraints.yaml || true)
    echo "Generated constraints: \${N_POCKET} pocket, \${N_CONTACTS} contact pairs"

    if [ "\${N_CONTACTS}" -eq 0 ] && [ "\${N_POCKET}" -eq 0 ]; then
        echo "ERROR: No constraints generated." >&2
        exit 1
    fi

    mkdir -p msa
    cat > msa/chain_A.a3m << EOF
>query
${receptor_seq}
EOF
    cat > msa/chain_B.a3m << EOF
>query
${effector_seq}
EOF

    cat > input.yaml << YAMLEOF
version: 1
sequences:
  - protein:
      id: ${rec_chain}
      sequence: ${receptor_seq}
      msa: \${PWD}/msa/chain_A.a3m
  - protein:
      id: ${eff_chain}
      sequence: ${effector_seq}
      msa: \${PWD}/msa/chain_B.a3m
YAMLEOF

    cat constraints.yaml >> input.yaml

    echo "Full input.yaml (first 40 lines):"
    head -40 input.yaml
    TOTAL_LINES=\$(wc -l < input.yaml)
    echo "  ... (\${TOTAL_LINES} lines total, \${N_CONTACTS} contact constraints)"
    echo ""

    echo "=== Validating YAML ==="
    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/validate_boltz_yaml.py input.yaml --model boltz2

    bash ${projectDir}/bin/run_boltz_seeds.sh \\
        input.yaml ${params.benchmark_container} ${projectDir} "Boltz-2 constrained"
    """
}

/*
 * =============================================================================
 * Template variants (Run 2)
 * =============================================================================
 * The three Boltz-2 variants above, re-run with the effector chain of the
 * reference complex supplied as a structural template (force: true,
 * threshold 1.0, enforced by --use_potentials).  The template block matches
 * the one structure-negative-steering's engine writes, so these variants and the
 * steering runs give the model the same structural information.
 *
 * These are NOT part of the model-vs-model comparison — every variant in that
 * comparison is template-free.  Their purpose is the within-model measurement
 * of what the known effector fold is worth: pair each against its
 * template-free twin (boltz2, boltz2_msa, boltz2_constrained), which differs
 * in nothing else.
 */

/* BOLTZ2 plus the effector template. Single-sequence MSA, no constraints. */
process BOLTZ2_TEMPLATE {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/boltz2_template", mode: 'copy',
        pattern: '{input.yaml,all_outputs/**}'

    input:
    val  receptor_seq
    val  effector_seq
    val  rec_chain
    val  eff_chain
    path effector_template_cif

    output:
    path "all_outputs",  emit: prediction_dir
    path "input.yaml"

    script:
    """
    set -euo pipefail

    mkdir -p msa
    cat > msa/chain_A.a3m << EOF
>query
${receptor_seq}
EOF
    cat > msa/chain_B.a3m << EOF
>query
${effector_seq}
EOF

    cat > input.yaml << YAMLEOF
version: 1
sequences:
  - protein:
      id: ${rec_chain}
      sequence: ${receptor_seq}
      msa: \${PWD}/msa/chain_A.a3m
  - protein:
      id: ${eff_chain}
      sequence: ${effector_seq}
      msa: \${PWD}/msa/chain_B.a3m
      templates:
        - cif: \${PWD}/${effector_template_cif}
          chain_id: ${eff_chain}
          template_id: A
          force: true
          threshold: 1.0
YAMLEOF

    echo "Input YAML:"
    cat input.yaml
    echo ""

    bash ${projectDir}/bin/run_boltz_seeds.sh \\
        input.yaml ${params.benchmark_container} ${projectDir} "Boltz-2 (template)"
    """
}

/* BOLTZ2_MSA plus the effector template. */
process BOLTZ2_MSA_TEMPLATE {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/boltz2_msa_template", mode: 'copy',
        pattern: '{input.yaml,all_outputs/**}'

    input:
    val  receptor_seq
    val  effector_seq
    val  rec_chain
    val  eff_chain
    path chain_a_a3m
    path chain_b_a3m
    path effector_template_cif

    output:
    path "all_outputs",  emit: prediction_dir
    path "input.yaml"

    script:
    """
    set -euo pipefail

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
      templates:
        - cif: \${PWD}/${effector_template_cif}
          chain_id: ${eff_chain}
          template_id: A
          force: true
          threshold: 1.0
YAMLEOF

    echo "Input YAML:"
    cat input.yaml
    echo ""

    bash ${projectDir}/bin/run_boltz_seeds.sh \\
        input.yaml ${params.benchmark_container} ${projectDir} "Boltz-2 MSA (template)"
    """
}

/*
 * BOLTZ2_CONSTRAINED plus the effector template — pocket + dense contact
 * constraints AND the known effector fold.  This is the combination the
 * steering engine's constrained runs use.
 */
process BOLTZ2_CONSTRAINED_TEMPLATE {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/boltz2_constrained_template", mode: 'copy',
        pattern: '{input.yaml,constraints.yaml,all_outputs/**}'

    input:
    val  receptor_seq
    val  effector_seq
    val  rec_chain
    val  eff_chain
    path reference_pdb
    path effector_template_cif

    output:
    path "all_outputs",      emit: prediction_dir
    path "input.yaml"
    path "constraints.yaml"

    script:
    def contact_cutoff      = 10.0
    def contact_max         = 50
    def contact_tolerance   = 0.0
    def pocket_cutoff       = 8.0
    def pocket_max_distance = 8.0
    """
    set -euo pipefail

    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/extract_constraints_boltz2.py \\
            ${reference_pdb} ${rec_chain} ${eff_chain} \\
            ${contact_cutoff} ${contact_max} ${contact_tolerance} \\
            ${pocket_cutoff} ${pocket_max_distance} \\
        > constraints.yaml

    N_CONTACTS=\$(grep -c "^  - contact:" constraints.yaml || true)
    N_POCKET=\$(grep -c "^  - pocket:" constraints.yaml || true)
    echo "Generated constraints: \${N_POCKET} pocket, \${N_CONTACTS} contact pairs"

    if [ "\${N_CONTACTS}" -eq 0 ] && [ "\${N_POCKET}" -eq 0 ]; then
        echo "ERROR: No constraints generated." >&2
        exit 1
    fi

    mkdir -p msa
    cat > msa/chain_A.a3m << EOF
>query
${receptor_seq}
EOF
    cat > msa/chain_B.a3m << EOF
>query
${effector_seq}
EOF

    cat > input.yaml << YAMLEOF
version: 1
sequences:
  - protein:
      id: ${rec_chain}
      sequence: ${receptor_seq}
      msa: \${PWD}/msa/chain_A.a3m
  - protein:
      id: ${eff_chain}
      sequence: ${effector_seq}
      msa: \${PWD}/msa/chain_B.a3m
      templates:
        - cif: \${PWD}/${effector_template_cif}
          chain_id: ${eff_chain}
          template_id: A
          force: true
          threshold: 1.0
YAMLEOF

    cat constraints.yaml >> input.yaml

    echo "=== Validating YAML ==="
    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/validate_boltz_yaml.py input.yaml --model boltz2

    bash ${projectDir}/bin/run_boltz_seeds.sh \\
        input.yaml ${params.benchmark_container} ${projectDir} "Boltz-2 constrained (template)"
    """
}

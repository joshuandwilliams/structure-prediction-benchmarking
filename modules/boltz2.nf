/*
 * =============================================================================
 * Boltz-2 module
 * =============================================================================
 * Three variants sharing the same unified `boltz predict` CLI (no --model
 * flag needed — Boltz-2 is the default in the unified package):
 *
 *   BOLTZ2              — no MSA, no constraints (single-sequence baseline)
 *   BOLTZ2_MSA          — with colabfold_search a3m MSA
 *   BOLTZ2_CONSTRAINED  — pocket + dense contact constraints from reference
 *
 * Gotchas preserved from the original bench scripts:
 *   1. --no_kernels is required because cuequivariance_torch is absent from
 *      Boltz1_Boltz2_Chai1.img.
 *   2. `boltz predict` exits 0 even on silent parse failures — but only the
 *      Boltz-1 scripts added explicit PDB-count guards.  For Boltz-2 we
 *      follow the same defensive pattern for consistency.
 *   3. Seeds are run SERIALLY inside the process (5 seeds: 42, 123, 456,
 *      789, 1024).
 */


/*
 * BOLTZ2
 * ------
 * Unconstrained Boltz-2 baseline: single-sequence MSA, no constraints.
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

    # ─── Single-sequence MSAs (query only) ────────────────────────────────
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

    run_seed() {
        local seed="\$1"
        local out_dir="\$2"

        singularity exec --nv \\
            --bind \${PWD}:\${PWD} \\
            ${params.benchmark_container} boltz predict \\
                input.yaml \\
                --out_dir "\${out_dir}" \\
                --recycling_steps 20 \\
                --diffusion_samples 5 \\
                --sampling_steps 20 \\
                --seed "\${seed}" \\
                --num_workers 0 \\
                --output_format pdb \\
                --write_full_pae \\
                --no_kernels \\
                --override

        local n_pdb
        n_pdb=\$(find "\${out_dir}" -name "*.pdb" 2>/dev/null | wc -l)
        if [ "\${n_pdb}" -eq 0 ]; then
            echo "ERROR: Seed \${seed} produced no PDB files." >&2
            exit 1
        fi
        echo "  Seed \${seed}: \${n_pdb} PDB(s) OK"
    }

    run_seed 42 output
    for SEED in 123 456 789 1024; do
        echo "Running Boltz-2 with seed \${SEED}..."
        run_seed "\${SEED}" "output_seed\${SEED}"
    done

    # ─── Aggregate outputs ────────────────────────────────────────────────
    bash ${projectDir}/bin/aggregate_seed_outputs.sh pdb npz json

    echo "Aggregated PDB files:"
    find all_outputs -name "*.pdb" | sort
    """
}


/*
 * BOLTZ2_MSA
 * ----------
 * Boltz-2 with the colabfold_search a3m MSAs shared from the COLABFOLD_SEARCH
 * process.  Identical to BOLTZ2 except for the MSA source.
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

    run_seed() {
        local seed="\$1"
        local out_dir="\$2"

        singularity exec --nv \\
            --bind \${PWD}:\${PWD} \\
            ${params.benchmark_container} boltz predict \\
                input.yaml \\
                --out_dir "\${out_dir}" \\
                --recycling_steps 20 \\
                --diffusion_samples 5 \\
                --sampling_steps 20 \\
                --seed "\${seed}" \\
                --num_workers 0 \\
                --output_format pdb \\
                --write_full_pae \\
                --use_potentials \\
                --no_kernels \\
                --override

        local n_pdb
        n_pdb=\$(find "\${out_dir}" -name "*.pdb" 2>/dev/null | wc -l)
        if [ "\${n_pdb}" -eq 0 ]; then
            echo "ERROR: Seed \${seed} produced no PDB files." >&2
            exit 1
        fi
        echo "  Seed \${seed}: \${n_pdb} PDB(s) OK"
    }

    run_seed 42 output
    for SEED in 123 456 789 1024; do
        echo "Running Boltz-2 MSA with seed \${SEED}..."
        run_seed "\${SEED}" "output_seed\${SEED}"
    done

    bash ${projectDir}/bin/aggregate_seed_outputs.sh pdb npz json

    echo "Aggregated PDB files:"
    find all_outputs -name "*.pdb" | sort
    """
}


/*
 * BOLTZ2_CONSTRAINED
 * ------------------
 * Boltz-2 with pocket + dense contact constraints derived from the
 * reference complex via bin/extract_constraints_boltz2.py.  No template
 * — constraints only — so the reference structure is never directly
 * given to the model (the reference IS the answer we are scoring against).
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

    # ─── Extract constraints from reference PDB ───────────────────────────
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

    # ─── Prepare Boltz YAML ───────────────────────────────────────────────
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

    # ─── Validate YAML ────────────────────────────────────────────────────
    echo "=== Validating YAML ==="
    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/validate_boltz_yaml.py input.yaml --model boltz2

    # ─── Run predictions ──────────────────────────────────────────────────
    run_seed() {
        local seed="\$1"
        local out_dir="\$2"

        singularity exec --nv \\
            --bind \${PWD}:\${PWD} \\
            ${params.benchmark_container} boltz predict \\
                input.yaml \\
                --out_dir "\${out_dir}" \\
                --recycling_steps 20 \\
                --diffusion_samples 5 \\
                --sampling_steps 20 \\
                --seed "\${seed}" \\
                --num_workers 0 \\
                --output_format pdb \\
                --write_full_pae \\
                --use_potentials \\
                --no_kernels \\
                --override

        local n_pdb
        n_pdb=\$(find "\${out_dir}" -name "*.pdb" 2>/dev/null | wc -l)
        if [ "\${n_pdb}" -eq 0 ]; then
            echo "ERROR: Seed \${seed} produced no PDB files." >&2
            exit 1
        fi
        echo "  Seed \${seed}: \${n_pdb} PDB(s) OK"
    }

    run_seed 42 output
    for SEED in 123 456 789 1024; do
        echo "Running Boltz-2 constrained with seed \${SEED}..."
        run_seed "\${SEED}" "output_seed\${SEED}"
    done

    # ─── Aggregate outputs ────────────────────────────────────────────────
    bash ${projectDir}/bin/aggregate_seed_outputs.sh pdb npz json

    echo "Aggregated PDB files:"
    find all_outputs -name "*.pdb" | sort
    """
}

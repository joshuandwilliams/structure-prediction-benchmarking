/*
 * =============================================================================
 * Boltz-1 module
 * =============================================================================
 * Three variants sharing the same unified `boltz predict --model boltz1` CLI:
 *
 *   BOLTZ1              — no MSA, no constraints (single-sequence baseline)
 *   BOLTZ1_MSA          — with colabfold_search a3m MSA
 *   BOLTZ1_CONSTRAINED  — pocket constraint only (Boltz-1 schema enforces
 *                         max_distance==6.0 and rejects contact constraints)
 *
 * Shared gotchas preserved from the original bench scripts:
 *   1. `boltz predict` exits 0 even on silent parse failures — the process
 *      loop verifies that each seed produced at least one PDB.
 *   2. --no_kernels is required because cuequivariance_torch is absent
 *      from benchmark_models.img.
 *   3. Seeds are run SERIALLY inside the process (five seeds: 42, 123, 456,
 *      789, 1024) to match the original bench scripts and respect the small
 *      jic-gpu queue.
 *   4. After all seeds complete, the outputs are aggregated into
 *      all_outputs/<seed_tag>/... so compute_metrics.py can recursively glob
 *      across every seed in one pass.
 *
 * The prediction_dir output is the aggregated staging dir and is what the
 * downstream metrics process consumes.
 */


/*
 * BOLTZ1
 * ------
 * Unconstrained Boltz-1 baseline: single-sequence MSA (just the query),
 * no constraints, 5 seeds × 5 diffusion samples = 25 structures.
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

    # ─── Input: bare-chain FASTA pointing at single-sequence "MSAs" ─────
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

    # ─── Run predictions serially across 5 seeds ──────────────────────────
    run_seed() {
        local seed="\$1"
        local out_dir="\$2"

        singularity exec --nv \\
            --bind \${PWD}:\${PWD} \\
            ${params.benchmark_container} boltz predict \\
                input.fasta \\
                --model boltz1 \\
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
            echo "ERROR: Seed \${seed} produced no PDB files — boltz silently skipped input." >&2
            exit 1
        fi
        echo "  Seed \${seed}: \${n_pdb} PDB(s) OK"
    }

    run_seed 42 output
    for SEED in 123 456 789 1024; do
        echo "Running Boltz-1 with seed \${SEED}..."
        run_seed "\${SEED}" "output_seed\${SEED}"
    done

    # ─── Aggregate all seed outputs into one staging dir ──────────────────
    bash ${projectDir}/bin/aggregate_seed_outputs.sh pdb npz json

    echo "Aggregated PDB files:"
    find all_outputs -name "*.pdb" | sort
    """
}


/*
 * BOLTZ1_MSA
 * ----------
 * Boltz-1 with a real MSA.  Receives chain_A.a3m and chain_B.a3m from
 * the shared COLABFOLD_SEARCH process and uses them directly in the
 * Boltz YAML input.
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

    # ─── Stage the shared MSAs locally so Boltz can read them ─────────────
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

    run_seed() {
        local seed="\$1"
        local out_dir="\$2"

        singularity exec --nv \\
            --bind \${PWD}:\${PWD} \\
            ${params.benchmark_container} boltz predict \\
                input.yaml \\
                --model boltz1 \\
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
            echo "ERROR: Seed \${seed} produced no PDB files — boltz silently skipped input." >&2
            exit 1
        fi
        echo "  Seed \${seed}: \${n_pdb} PDB(s) OK"
    }

    run_seed 42 output
    for SEED in 123 456 789 1024; do
        echo "Running Boltz-1 MSA with seed \${SEED}..."
        run_seed "\${SEED}" "output_seed\${SEED}"
    done

    # ─── Aggregate outputs ────────────────────────────────────────────────
    bash ${projectDir}/bin/aggregate_seed_outputs.sh pdb npz json

    echo "Aggregated PDB files:"
    find all_outputs -name "*.pdb" | sort
    """
}


/*
 * BOLTZ1_CONSTRAINED
 * ------------------
 * Boltz-1 with a pocket constraint derived from the reference complex.
 *
 * NOTE: Boltz-1 schema.py enforces max_distance == 6.0 for ALL constraint
 * blocks and silently rejects contact constraints, so only a pocket block
 * is generated. bin/extract_constraints_boltz1.py enforces the exactly-6.0
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

    # ─── Extract pocket constraint from reference PDB ─────────────────────
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

    # ─── Prepare Boltz YAML (single-sequence MSAs + pocket constraint) ───
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

    # ─── Validate YAML (max_distance==6.0, no contact blocks) ────────────
    echo "=== Validating YAML ==="
    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/validate_boltz_yaml.py input.yaml --model boltz1

    # ─── Run predictions ──────────────────────────────────────────────────
    run_seed() {
        local seed="\$1"
        local out_dir="\$2"

        singularity exec --nv \\
            --bind \${PWD}:\${PWD} \\
            ${params.benchmark_container} boltz predict \\
                input.yaml \\
                --model boltz1 \\
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
            echo "  boltz predict exited 0 but silently skipped the input." >&2
            echo "  Check stderr for 'Failed to process' or 'Max distance != 6.0'." >&2
            exit 1
        fi
        echo "  Seed \${seed}: \${n_pdb} PDB(s) OK"
    }

    run_seed 42 output
    for SEED in 123 456 789 1024; do
        echo "Running Boltz-1 constrained with seed \${SEED}..."
        run_seed "\${SEED}" "output_seed\${SEED}"
    done

    # ─── Aggregate outputs ────────────────────────────────────────────────
    bash ${projectDir}/bin/aggregate_seed_outputs.sh pdb npz json

    echo "Aggregated PDB files:"
    find all_outputs -name "*.pdb" | sort
    """
}

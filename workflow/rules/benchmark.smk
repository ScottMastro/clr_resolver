# benchmark: simulate CLR reads from known region assemblies, run the resolve
# pipeline on them, and score read phasing, copy count, and bp accuracy.
#
# Depends on: resolve (the sim_<sample> dataset flows through every resolve
# rule unchanged via the DAG; no resolve rule is redefined here).
#
# Pipeline (5 rules; the resolve sub-pipeline runs in between):
#   simulate_clr -> [resolve on sim_<sample>] ->
#   benchmark_truth ; score_phasing ; measure_bp_accuracy ->
#   summarize_benchmark
#
# Run:
#   snakemake -s workflow/Snakefile benchmark

BM_SCRIPTS = "workflow/scripts/benchmark"
BENCH = f"{WORK_DIR}/benchmark"


rule benchmark:
    input:
        expand(f"{RESULTS_DIR}/benchmark/{{sample}}/summary.tsv",
               sample=SAMPLES),
        expand(f"{RESULTS_DIR}/benchmark/{{sample}}/summary.txt",
               sample=SAMPLES),


rule simulate_clr:
    """Simulate diploid PacBio CLR reads from a sample's region assemblies."""
    input:
        script=f"{BM_SCRIPTS}/simulate_clr.sh",
        hap1=lambda w: sample_hap_fa(w.sample, 1),
        hap2=lambda w: sample_hap_fa(w.sample, 2),
    output:
        reads=f"{WORK_DIR}/sim/{{sample}}/clr_reads.fq",
        maf1=f"{WORK_DIR}/sim/{{sample}}/hap1_0001.maf.gz",
        maf2=f"{WORK_DIR}/sim/{{sample}}/hap2_0001.maf.gz",
        ref1=f"{WORK_DIR}/sim/{{sample}}/hap1_0001.ref",
        ref2=f"{WORK_DIR}/sim/{{sample}}/hap2_0001.ref",
    params:
        outdir=lambda w: f"{WORK_DIR}/sim/{w.sample}",
    threads: 4
    resources:
        runtime=30,
    log:
        "logs/simulate_clr/{sample}.log",
    shell:
        r"""
        source config/tools.sh
        source config/params.sh
        bash {input.script} {input.hap1} {input.hap2} {params.outdir} \
            "$SIM_DEPTH" "$SIM_LEN_MEAN" "$SIM_LEN_SD" \
            "$PBSIM" "$PBSIM_MODEL" &> {log}
        """


rule benchmark_truth:
    """Per-read truth: which orange copy each simulated read came from."""
    input:
        script=f"{BM_SCRIPTS}/clr_truth.py",
        maf1=rules.simulate_clr.output.maf1,
        maf2=rules.simulate_clr.output.maf2,
        gfa=GRAPH_GFA,
        node_color=NODE_COLOR,
    output:
        truth=f"{BENCH}/{{sample}}/truth.tsv",
    params:
        color="O",
    log:
        "logs/benchmark_truth/{sample}.log",
    shell:
        r"""
        python3 {input.script} --maf-h1 {input.maf1} --maf-h2 {input.maf2} \
            --gfa {input.gfa} --node-color {input.node_color} \
            --sample {wildcards.sample} --color {params.color} \
            --out {output.truth} &> {log}
        """


rule score_phasing:
    """Read phasing accuracy of the resolve clustering versus truth."""
    input:
        script=f"{BM_SCRIPTS}/score.py",
        pred=f"{WORK_DIR}/untangle/sim_{{sample}}/pred.tsv",
        truth=rules.benchmark_truth.output.truth,
    output:
        phasing=f"{BENCH}/{{sample}}/phasing.tsv",
    log:
        "logs/score_phasing/{sample}.log",
    shell:
        r"""
        python3 {input.script} --pred {input.pred} --truth {input.truth} \
            --out-tsv {output.phasing} &> {log}
        """


rule measure_bp_accuracy:
    """Align resolved copies back to the truth haplotypes for per-copy identity."""
    input:
        script=f"{BM_SCRIPTS}/compare_to_truth.py",
        copies=f"{RESULTS_DIR}/resolve/sim_{{sample}}/copies.fasta",
        ref1=rules.simulate_clr.output.ref1,
        ref2=rules.simulate_clr.output.ref2,
    output:
        bp=f"{BENCH}/{{sample}}/bp_accuracy.tsv",
    log:
        "logs/measure_bp_accuracy/{sample}.log",
    shell:
        r"""
        source config/tools.sh
        python3 {input.script} --copies {input.copies} \
            --truth-h1 {input.ref1} --truth-h2 {input.ref2} \
            --minimap2 "$MINIMAP2" --out {output.bp} &> {log}
        """


rule summarize_benchmark:
    """Combine phasing, copy count, and bp accuracy into the sample summary."""
    input:
        script=f"{BM_SCRIPTS}/summarize_benchmark.py",
        phasing=rules.score_phasing.output.phasing,
        bp=rules.measure_bp_accuracy.output.bp,
        manifest=f"{RESULTS_DIR}/resolve/sim_{{sample}}/copies.tsv",
        truth=rules.benchmark_truth.output.truth,
    output:
        summary_tsv=f"{RESULTS_DIR}/benchmark/{{sample}}/summary.tsv",
        summary_txt=f"{RESULTS_DIR}/benchmark/{{sample}}/summary.txt",
    log:
        "logs/summarize_benchmark/{sample}.log",
    shell:
        r"""
        python3 {input.script} --phasing {input.phasing} --bp {input.bp} \
            --copies-manifest {input.manifest} --truth {input.truth} \
            --sample {wildcards.sample} --out-tsv {output.summary_tsv} \
            --out-txt {output.summary_txt} &> {log}
        """

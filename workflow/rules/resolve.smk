# resolve: separate the orange (MUC20) paralog copies in a CLR read set and
# assemble each copy with Flye.
#
# Depends on: nothing in workflow/ (common.smk constants only).
#
# Pipeline (8 rules):
#   align_reads_to_graph -> clr_extract_reads -> clr_call_markers ->
#   build_marker_matrix -> cluster_copies -> split_pools (checkpoint) ->
#   flye_assemble_pool -> collect_copies
#
# Run:
#   snakemake -s workflow/Snakefile resolve

SCRIPTS = "workflow/scripts/resolve"
WORK = WORK_DIR
RESULTS = RESULTS_DIR
COLOR = "O"


rule resolve:
    input:
        expand(f"{RESULTS}/resolve/{{dataset}}/copies.fasta", dataset=DATASETS),
        expand(f"{RESULTS}/resolve/{{dataset}}/copies.tsv", dataset=DATASETS),


rule align_reads_to_graph:
    """Align the CLR reads to the pangenome graph, producing a gzipped GAF."""
    input:
        reads=resolve_reads_fastq,
        gfa=GRAPH_GFA,
    output:
        gaf=f"{WORK}/align/{{dataset}}.gaf.gz",
    params:
        preset="vg",
    threads: 8
    resources:
        mem_mb=8000,
        runtime=120,
    log:
        "logs/align_reads_to_graph/{dataset}.log",
    shell:
        r"""
        source config/tools.sh
        # GraphAligner has no gzip output: write the plain GAF, then compress.
        gaf_gz="{output.gaf}"
        plain="${{gaf_gz%.gz}}"
        "$GRAPHALIGNER" -g {input.gfa} -f {input.reads} -a "$plain" \
            -x {params.preset} -t {threads} &> {log}
        gzip -f "$plain"
        """


rule clr_extract_reads:
    """Per-read table of the CLR reads scoped to the orange colour block."""
    input:
        script=f"{SCRIPTS}/clr_reads.py",
        gaf=rules.align_reads_to_graph.output.gaf,
        node_color=NODE_COLOR,
    output:
        reads=f"{WORK}/untangle/{{dataset}}/reads.tsv",
    params:
        color=COLOR,
    log:
        "logs/clr_extract_reads/{dataset}.log",
    shell:
        r"""
        source config/params.sh
        python3 {input.script} --gaf {input.gaf} \
            --node-color {input.node_color} --color {params.color} \
            --min-c-nodes "$MIN_C_NODES" --out {output.reads} &> {log}
        """


rule clr_call_markers:
    """Bubble allele calls for the orange block from the CLR alignments."""
    input:
        script=f"{SCRIPTS}/clr_markers.py",
        gaf=rules.align_reads_to_graph.output.gaf,
        bubbles=BUBBLES,
        node_color=NODE_COLOR,
        gfa=GRAPH_GFA,
    output:
        markers=f"{WORK}/untangle/{{dataset}}/markers.tsv",
        calls=f"{WORK}/untangle/{{dataset}}/read_markers.tsv",
    params:
        color=COLOR,
    resources:
        mem_mb=8000,
    log:
        "logs/clr_call_markers/{dataset}.log",
    shell:
        r"""
        python3 {input.script} --gaf {input.gaf} --bubbles {input.bubbles} \
            --node-color {input.node_color} --gfa {input.gfa} \
            --color {params.color} --out-markers {output.markers} \
            --out-calls {output.calls} &> {log}
        """


rule build_marker_matrix:
    """Per-barcode marker matrix and structural-coverage table."""
    input:
        script=f"{SCRIPTS}/build_matrix.py",
        reads=rules.clr_extract_reads.output.reads,
        read_markers=rules.clr_call_markers.output.calls,
        markers=rules.clr_call_markers.output.markers,
        node_color=NODE_COLOR,
    output:
        matrix=f"{WORK}/untangle/{{dataset}}/matrix.tsv",
        struct=f"{WORK}/untangle/{{dataset}}/struct.tsv",
    params:
        color=COLOR,
    log:
        "logs/build_marker_matrix/{dataset}.log",
    shell:
        r"""
        python3 {input.script} --reads {input.reads} \
            --read-markers {input.read_markers} --markers {input.markers} \
            --node-color {input.node_color} --color {params.color} \
            --out-matrix {output.matrix} --out-struct {output.struct} &> {log}
        """


rule cluster_copies:
    """Cluster barcodes into orange paralog copies (whatshap polyphase).

    Run with the whatshap-env interpreter, not `python3`: the script imports
    whatshap.polyphase, which is installed only in that environment."""
    input:
        script=f"{SCRIPTS}/cluster_whatshap.py",
        matrix=rules.build_marker_matrix.output.matrix,
    output:
        pred=f"{WORK}/untangle/{{dataset}}/pred.tsv",
    params:
        ploidy=lambda w: dataset_ploidy(w.dataset),
    log:
        "logs/cluster_copies/{dataset}.log",
    shell:
        r"""
        source config/tools.sh
        source config/params.sh
        "$WHATSHAP_PYTHON" {input.script} --matrix {input.matrix} \
            --out {output.pred} --collapse --ploidy {params.ploidy} \
            --err "$WHATSHAP_ERR" &> {log}
        """


checkpoint split_pools:
    """Split the read FASTQ into one FASTQ per resolved copy.

    A checkpoint: the pool count is data-dependent (at most the ploidy), so
    the downstream Flye fan-out is only known once this has run."""
    input:
        script=f"{SCRIPTS}/split_pools.py",
        pred=rules.cluster_copies.output.pred,
        reads=resolve_reads_fastq,
    output:
        pools_dir=directory(f"{WORK}/pools/{{dataset}}"),
    log:
        "logs/split_pools/{dataset}.log",
    shell:
        r"""
        python3 {input.script} --pred {input.pred} --reads {input.reads} \
            --out-dir {output.pools_dir} &> {log}
        """


rule flye_assemble_pool:
    """Assemble one per-copy read pool with Flye into a single contig."""
    input:
        pool=f"{WORK}/pools/{{dataset}}/{{pool}}.fq",
    output:
        asm=f"{WORK}/flye/{{dataset}}/{{pool}}/assembly.fasta",
    params:
        outdir=lambda w: f"{WORK}/flye/{w.dataset}/{w.pool}",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        "logs/flye_assemble_pool/{dataset}.{pool}.log",
    shell:
        r"""
        source config/tools.sh
        # Flye locates its bundled flye-minimap2 / flye-samtools on PATH;
        # add its env bin dir so calling the binary by absolute path works.
        export PATH="$(dirname "$FLYE"):$PATH"
        # A pool too thin to assemble (Flye: "No disjointigs were assembled")
        # is not a fatal error -- it means the cluster held too few reads to
        # be a real copy. Emit an empty assembly so collect_copies records it
        # as a 0-contig pool instead of aborting the whole run.
        if "$FLYE" --pacbio-raw {input.pool} --out-dir {params.outdir} \
                --threads {threads} &> {log}; then
            :
        else
            echo "[WARN] Flye failed on {wildcards.dataset}/{wildcards.pool} --"\
                 "emitting empty assembly" >> {log}
            mkdir -p {params.outdir}
            : > {output.asm}
        fi
        """


def aggregate_pool_assemblies(wildcards):
    """Resolve the split_pools checkpoint, list one assembly.fasta per pool."""
    ckpt = checkpoints.split_pools.get(dataset=wildcards.dataset)
    pools = glob_wildcards(
        os.path.join(ckpt.output.pools_dir, "{pool}.fq")).pool
    return expand(f"{WORK}/flye/{{dataset}}/{{pool}}/assembly.fasta",
                  dataset=wildcards.dataset, pool=sorted(set(pools)))


rule collect_copies:
    """Gather every per-pool Flye assembly into copies.fasta + a manifest."""
    input:
        script=f"{SCRIPTS}/collect_copies.py",
        asms=aggregate_pool_assemblies,
    output:
        fasta=f"{RESULTS}/resolve/{{dataset}}/copies.fasta",
        manifest=f"{RESULTS}/resolve/{{dataset}}/copies.tsv",
    log:
        "logs/collect_copies/{dataset}.log",
    shell:
        r"""
        python3 {input.script} --assemblies {input.asms} \
            --out-fasta {output.fasta} --out-manifest {output.manifest} &> {log}
        """

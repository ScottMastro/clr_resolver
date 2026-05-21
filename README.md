# clr_resolve

Isolated, end-to-end pipeline that resolves the **orange / MUC20** paralog
copies of the MUC pangenome locus from PacBio **CLR** long reads, and a
benchmark that measures how well it does.

MUC20 sits in a tandem segmental duplication: each haplotype carries several
near-identical orange copies. `clr_resolve` separates the CLR reads by copy,
assembles each copy independently, and so reconstructs the duplication that a
whole-region assembler would otherwise collapse.

## What it does

### `resolve` — production

Takes a set of pre-filtered CLR reads and:

1. **aligns** them to the MUC pangenome graph (GraphAligner → GAF);
2. **separates the orange paralog copies** — the orange untangling chain
   `clr_reads` → `clr_markers` → `build_matrix` → `cluster_whatshap`, which
   calls graph-bubble alleles per read and clusters reads by copy with
   WhatsHap polyphase;
3. **splits** the reads into one pool per resolved copy;
4. **assembles** each pool with Flye — one contig per copy.

Output: `results/resolve/<dataset>/copies.fasta` (the per-copy contigs) and
`copies.tsv` (a manifest: copy id, source pool, contig length).

### `benchmark` — testing

For each sample, simulates CLR reads (pbsim3) from a known pair of region
assemblies, runs `resolve` on them, and scores the result against truth:

- **read phasing accuracy** — fraction of reads clustered to their true copy;
- **copy count correct** — resolved copies vs distinct truth copies;
- **bp accuracy** — each resolved copy aligned back to the truth assembly
  (minimap2, identity and length ratio).

Output: `results/benchmark/<sample>/summary.tsv` (one machine-loadable row)
and `summary.txt` (a human-readable report).

The benchmark reuses every `resolve` rule unchanged: a sample's simulated read
set enters the pipeline as the dataset `sim_<sample>`.

#### Validation result — HG00133

| Metric | Result |
|---|---|
| Read phasing accuracy | 96.4% (275 reads) |
| Copy count | 6 resolved / 6 expected — correct |
| bp accuracy | 98.0% mean identity, 6/6 copies one contig each |

## Layout

```
clr_resolve/
├── install.sh                  # create the conda envs (cold-machine setup)
├── workflow/
│   ├── Snakefile               # entry point; includes the rule modules
│   ├── rules/
│   │   ├── common.smk          # constants, wildcard constraints, input fns
│   │   ├── resolve.smk         # production pipeline (target rule `resolve`)
│   │   └── benchmark.smk       # testing pipeline (target rule `benchmark`)
│   ├── scripts/
│   │   ├── resolve/            # align → untangle → split → assemble scripts
│   │   └── benchmark/          # simulate, truth, score, summarise scripts
│   └── envs/                   # conda recreation recipes (one per tool)
├── config/
│   ├── workflow.yaml           # graph inputs, datasets, benchmark samples
│   ├── params.sh               # tunables (read cutoff, ploidy, sim depth …)
│   └── tools.sh                # site-specific tool paths
├── data/                       # node_color.tsv + graph.bubbles.json (the
│                               # orange painting / bubbles of the graph);
│                               # drop the graph GFA here too (see Setup)
├── results/                    # per-dataset copies, per-sample summaries
├── work/                       # cache of intermediates (gitignored)
└── logs/                       # per-rule logs (gitignored)
```

## Setup

### 1. Tools

`clr_resolve` calls external tools (Snakemake, GraphAligner, Flye, WhatsHap,
pbsim3, minimap2) through small conda environments. On a fresh machine:

```bash
./install.sh
```

This creates one conda env per `workflow/envs/*.yaml` — including a
`snakemake` env for the workflow engine — skipping any that already exist.
`config/tools.sh` then resolves each binary from those envs (`conda info
--base` is honoured, so any conda install location works); GraphAligner and
minimap2 fall back to a `PATH` lookup if their env is absent. Any tool path
can be overridden by pre-setting its variable (`FLYE`, `PBSIM`, …) in the
environment.

`conda` must already be installed.

### 2. Graph

The pangenome graph GFA (`muc-hprc-pggb-v.1.0.gfa`, ~268 MB) is too large to
ship in the repo. Obtain it and place it at `data/muc-hprc-pggb-v.1.0.gfa`
(or point `graph_gfa` in `config/workflow.yaml` at its location). The matching
orange painting (`node_color.tsv`) and bubble decomposition
(`graph.bubbles.json`) are vendored under `data/`.

## Run

From the `clr_resolve/` directory, with the `snakemake` env active:

```bash
conda activate snakemake

# resolve copies for a production read set
snakemake -s workflow/Snakefile --cores 8 resolve

# benchmark the configured samples
snakemake -s workflow/Snakefile --cores 8 benchmark
```

To resolve a real CLR read set, add it to `datasets:` in
`config/workflow.yaml` (`<name>: <path-to-fastq>`); the dataset name must not
start with `sim_`, which is reserved for benchmark-simulated reads.

To benchmark a sample, add it to `samples:` and give its two region-assembly
FASTAs under `sample_haps:` (see the example in `config/workflow.yaml`).

## Configuration

- **`config/workflow.yaml`** — the graph GFA, painting and bubble inputs; the
  `datasets:` map of production read sets; the benchmark `samples:` and their
  `sample_haps:` region assemblies.
- **`config/params.sh`** — `MIN_C_NODES` (orange-read admission cutoff),
  `PLOIDY` (expected copy count), `WHATSHAP_ERR` (clustering error rate), and
  the pbsim depth / read-length knobs. Rules `source` this file, so after
  editing it force the affected rule, e.g.
  `snakemake --forcerun cluster_copies benchmark`.
- **`config/tools.sh`** — per-host tool paths (see Setup).

## Scope

The pipeline ends at raw Flye contigs (≈96–99% identity per copy). Scaffolding
the copies into whole haplotypes and polishing to ~99.9% are deliberately out
of scope.

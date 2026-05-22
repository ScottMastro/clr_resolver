# common: shared constants, wildcard constraints, and dataset/sample input
# functions for the clr_resolve pipeline.
#
# Depends on: nothing in workflow/.
#
# Included first by the Snakefile, so its globals are visible in every
# analysis .smk. Holds no rules.

import os
import sys

from snakemake.exceptions import WorkflowError

WORK_DIR = "work"
RESULTS_DIR = "results"

GRAPH_GFA = config["graph_gfa"]
NODE_COLOR = config["node_color"]
BUBBLES = config["bubbles"]

# benchmark-simulated read sets enter pipeline A as `sim_<sample>` datasets
SIM_PREFIX = "sim_"

wildcard_constraints:
    dataset=r"[A-Za-z0-9_]+",
    sample=r"[A-Za-z0-9]+",
    pool=r"wh_[0-9]+",


def _items(value):
    """Comma/whitespace string or list from config -> sorted unique list."""
    if isinstance(value, (list, tuple)):
        items = list(value)
    else:
        items = str(value or "").replace(",", " ").split()
    return sorted(set(i for i in items if i))


SAMPLES = _items(config.get("samples", ""))
SAMPLE_HAPS = dict(config.get("sample_haps", {}))
USER_DATASETS = dict(config.get("datasets", {}))
DATASETS = sorted(set(USER_DATASETS) | {SIM_PREFIX + s for s in SAMPLES})

SAMPLE_PLOIDY = {s: int(k) for s, k in config.get("sample_ploidy", {}).items()}
DEFAULT_PLOIDY = int(config.get("ploidy", 6))

sys.stderr.write(
    f"[INFO] common: {len(DATASETS)} dataset(s) {DATASETS} · "
    f"{len(SAMPLES)} benchmark sample(s) {SAMPLES}\n")


def sim_sample(dataset):
    """`sim_<sample>` -> `<sample>`; None for a non-simulated dataset."""
    if dataset.startswith(SIM_PREFIX):
        return dataset[len(SIM_PREFIX):]
    return None


def resolve_reads_fastq(wildcards):
    """The CLR FASTQ feeding align_reads_to_graph for a {dataset}.

    sim_* datasets resolve to the simulate_clr output (which wires the DAG
    edge into the benchmark sub-pipeline); every other dataset resolves to a
    path declared in config `datasets:`."""
    s = sim_sample(wildcards.dataset)
    if s is not None:
        return f"{WORK_DIR}/sim/{s}/clr_reads.fq"
    if wildcards.dataset not in USER_DATASETS:
        raise WorkflowError(
            f"dataset '{wildcards.dataset}' is not a sim_* dataset and is not "
            f"declared in config datasets:; known = {sorted(USER_DATASETS)}")
    return USER_DATASETS[wildcards.dataset]


def dataset_ploidy(dataset):
    """Expected orange copy number (cluster_whatshap ploidy) for a {dataset}:
    the benchmark sample's painting-truth k where known, else the default."""
    s = sim_sample(dataset)
    if s is not None and s in SAMPLE_PLOIDY:
        return SAMPLE_PLOIDY[s]
    return DEFAULT_PLOIDY


def sample_hap_fa(sample, hap):
    """Region-assembly FASTA for haplotype `hap` (1 or 2) of a benchmark sample."""
    haps = SAMPLE_HAPS.get(sample)
    if not haps or len(haps) < 2:
        raise WorkflowError(
            f"sample '{sample}' needs a 2-entry sample_haps: config block")
    return haps[hap - 1]

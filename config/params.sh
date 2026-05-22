#!/usr/bin/env bash
# Shell-level tunables for clr_resolve. Sourced by rules; values are passed to
# scripts as explicit CLI flags, never read from the environment by a script.

# --- resolve: orange (MUC20) paralog separation ---
export MIN_C_NODES=6      # clr_reads.py: keep reads touching >= N orange nodes
# cluster_whatshap.py per-call error rate. With orange-only markers accuracy
# plateaus across [0.01, 0.07] (~96-97%) and collapses below 0.01; 0.07 is the
# whatshap default and the established orange setting.
export WHATSHAP_ERR=0.07
# Expected copy number (ploidy) is per-dataset -- see `ploidy` / `sample_ploidy`
# in config/workflow.yaml, resolved by dataset_ploidy() in common.smk.

# --- benchmark: CLR read simulation (pbsim3) ---
export SIM_DEPTH=20         # per-haplotype fold coverage
export SIM_LEN_MEAN=15000   # mean read length (bp)
export SIM_LEN_SD=8000      # read-length standard deviation (bp)

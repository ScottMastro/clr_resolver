#!/usr/bin/env bash
# Simulate diploid PacBio CLR reads from a pair of region assemblies (pbsim3).
#
# Usage:
#   simulate_clr.sh <HAP1_FA> <HAP2_FA> <OUT_DIR> <DEPTH> <LEN_MEAN> <LEN_SD> \
#                   <PBSIM> <PBSIM_MODEL>
#
# Writes into OUT_DIR: clr_reads.fq (both haplotypes' reads merged, renamed
# @h1_*/@h2_* so the source haplotype is recoverable), and pbsim's per-
# haplotype hap{1,2}_0001.maf.gz (true read spans) and hap{1,2}_0001.ref.
set -euo pipefail

HAP1=${1:?args}
HAP2=${2:?args}
OUT_DIR=${3:?args}
DEPTH=${4:?args}
LEN_MEAN=${5:?args}
LEN_SD=${6:?args}
PBSIM=${7:?args}
PBSIM_MODEL=${8:?args}

mkdir -p "$OUT_DIR"
ABS1=$(readlink -f "$HAP1")
ABS2=$(readlink -f "$HAP2")
cd "$OUT_DIR"

# pbsim3 emits <prefix>_0001.{fq.gz|fq|fastq}, .maf(.gz) and .ref; which of
# those are gzipped varies by build, so each is normalised below.
emit_fastq() {  # haplotype index -> stdout: the plain FASTQ
    local i=$1
    if [[ -f "hap${i}_0001.fq.gz" ]]; then
        zcat "hap${i}_0001.fq.gz"
    elif [[ -f "hap${i}_0001.fq" ]]; then
        cat "hap${i}_0001.fq"
    elif [[ -f "hap${i}_0001.fastq" ]]; then
        cat "hap${i}_0001.fastq"
    else
        echo "[ERROR] simulate_clr: no pbsim FASTQ for haplotype $i" >&2
        exit 1
    fi
}

i=0
for ABS in "$ABS1" "$ABS2"; do
    i=$((i + 1))
    "$PBSIM" --strategy wgs --method errhmm --errhmm "$PBSIM_MODEL" \
        --depth "$DEPTH" --length-mean "$LEN_MEAN" --length-sd "$LEN_SD" \
        --genome "$ABS" --prefix "hap${i}"
    if [[ -f "hap${i}_0001.maf" ]]; then
        gzip -f "hap${i}_0001.maf"
    fi
done

: > clr_reads.fq
for i in 1 2; do
    emit_fastq "$i" \
        | awk -v h="$i" 'NR%4==1{print "@h"h"_"substr($1,2)} NR%4!=1{print}' \
        >> clr_reads.fq
done

echo "[INFO] simulate_clr: $(($(wc -l < clr_reads.fq) / 4)) CLR reads" >&2
echo "[OK] simulate_clr: $OUT_DIR/clr_reads.fq" >&2

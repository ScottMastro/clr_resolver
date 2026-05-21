#!/usr/bin/env bash
# Create the conda environments clr_resolve needs (cold-machine setup).
#
# Usage:
#   ./install.sh
#
# Creates one conda env per workflow/envs/*.yaml (graphaligner, minimap2,
# flye, whatshap-env, pbsim), skipping any that already exist. After it
# finishes, confirm config/tools.sh resolves to the right binaries for this
# host -- by default it reads them from these envs under `conda info --base`.
set -euo pipefail
cd "$(dirname "$0")"

if ! command -v conda >/dev/null 2>&1; then
    echo "[ERROR] install: conda not found on PATH" >&2
    exit 1
fi

existing=$(conda env list | awk 'NF && $1 !~ /^#/ {print $1}')

for yaml in workflow/envs/*.yaml; do
    name=$(sed -n 's/^name:[[:space:]]*//p' "$yaml" | head -1)
    if [[ -z "$name" ]]; then
        echo "[WARN] install: $yaml has no name: field, skipped" >&2
        continue
    fi
    if grep -qx "$name" <<<"$existing"; then
        echo "[SKIP] install: env '$name' already exists" >&2
        continue
    fi
    echo "[INSTALL] install: creating env '$name' from $yaml" >&2
    conda env create -f "$yaml"
done

echo "[OK] install: conda envs ready -- review config/tools.sh" >&2

#!/usr/bin/env bash
# Site-specific tool resolution for clr_resolve. Replicators may pre-set any
# of these variables in the environment to override the defaults below.
#
# Defaults read each tool from the conda env that ./install.sh creates from
# workflow/envs/*.yaml. GraphAligner and minimap2 fall back to a PATH lookup
# when their env is absent (any install works); flye, whatshap, and pbsim are
# taken strictly from their env.

_conda_base="$(conda info --base 2>/dev/null || echo "$HOME/anaconda3")"

_env_bin() {  # <env> <binary> -> path inside that conda env
    echo "$_conda_base/envs/$1/bin/$2"
}

_resolve() {  # <env> <binary> -> env path if present, else a PATH lookup
    local p
    p="$(_env_bin "$1" "$2")"
    if [[ -x "$p" ]]; then
        echo "$p"
    else
        command -v "$2" 2>/dev/null || echo "$p"
    fi
}

export GRAPHALIGNER="${GRAPHALIGNER:-$(_resolve graphaligner GraphAligner)}"
export MINIMAP2="${MINIMAP2:-$(_resolve minimap2 minimap2)}"
export FLYE="${FLYE:-$(_env_bin flye flye)}"

# cluster_whatshap.py imports whatshap.polyphase -- run it with this python.
export WHATSHAP_PYTHON="${WHATSHAP_PYTHON:-$(_env_bin whatshap-env python)}"

# pbsim3 long-read simulator and its ERRHMM error model (benchmark only).
export PBSIM="${PBSIM:-$(_env_bin pbsim pbsim)}"
export PBSIM_MODEL="${PBSIM_MODEL:-$_conda_base/envs/pbsim/data/ERRHMM-RSII.model}"

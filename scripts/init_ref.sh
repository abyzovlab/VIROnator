#!/bin/bash

USAGE_1="Usage: `basename $0` -i input.fa -v viral.fa -o output.fa [-b /path/to/bwa]"
combined_refs="/research/labs/res-users/m277455/genomes/refs/combined_refs"

ARG_REF=""
ARG_OUT=""
ARG_HVR=""
BWA_BIN="bwa"  # default: auto-discover from PATH

# read input args
while [[ "$#" -gt 0 ]]; do case $1 in
  -i|--inref) ARG_REF="$2"; shift;;
  -v|--viral) ARG_HVR="$2"; shift;;
  -o|--outref) ARG_OUT="$2"; shift;;
  -b|--bwa) BWA_BIN="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

# Multi-stage auto-discovery for bwa executable
resolve_bwa() {
    local target="$1"
    if [ -n "$target" ] && command -v "$target" &>/dev/null; then
        echo "$target"
        return 0
    fi
    if command -v bwa &>/dev/null; then
        which bwa
        return 0
    fi
    for path in \
        "/usr/local/bin/bwa" \
        "/usr/bin/bwa" \
        "/opt/conda/bin/bwa" \
        "/opt/conda/envs/vironator/bin/bwa" \
        "$HOME/miniconda3/bin/bwa" \
        "$HOME/miniconda3/envs/vironator/bin/bwa" \
        "$HOME/anaconda3/bin/bwa" \
        "$HOME/anaconda3/envs/vironator/bin/bwa"
    do
        if [ -x "$path" ]; then
            echo "$path"
            return 0
        fi
    done
    local found
    found=$(find /opt /usr /home /root -name "bwa" -type f 2>/dev/null | head -n 1)
    if [ -n "$found" ] && [ -x "$found" ]; then
        echo "$found"
        return 0
    fi
    return 1
}

BWA_RESOLVED=$(resolve_bwa "${BWA_BIN}")
if [ $? -eq 0 ] && [ -n "${BWA_RESOLVED}" ]; then
    BWA_BIN="${BWA_RESOLVED}"
    echo "Using BWA binary: ${BWA_BIN}"
else
    echo "ERROR: Could not locate bwa executable." >&2
    exit 127
fi

if [ "$ARG_REF" == "" ]; then
  echo
  echo "-i input missing"
  echo $USAGE_1
  echo
  exit 1
fi
if [ "$ARG_OUT" == "" ]; then
  echo
  echo "-o input missing"
  echo $USAGE_1
  echo
  exit 1
fi
if [ "$ARG_HVR" == "" ]; then
  HVR="/research/labs/res-users/m277455/genomes/refs/viral_refs/HumanViral_Reference_02-07-2022.fa"
  echo
  echo "Using default viral reference sequences:"
  echo $HVR
  echo
else
  HVR=$ARG_HVR
  echo
  echo "Using user-specified viral reference sequences:"
  echo $HVR
  echo
fi

# scripts
make_json="/research/labs/res-users/m277455/exogene_python/make_exogene_json.py"

# samtools
samtools faidx $ARG_REF
ARG_REF_NO_EXTENSION=$(basename $ARG_REF | cut -f 1 -d '.')
n_contigs=$(wc -l < /research/labs/res-users/m277455/genomes/refs/human_refs/$ARG_REF_NO_EXTENSION.fa.fai)

ARG_OUT_NO_EXTENSION=$(basename $ARG_OUT | cut -f 1 -d '.')
cat $ARG_REF $HVR > $ARG_OUT
samtools faidx $ARG_OUT
python $make_json $HVR $n_contigs $combined_refs/$ARG_OUT_NO_EXTENSION.exogene.json
"${BWA_BIN}" index $ARG_OUT

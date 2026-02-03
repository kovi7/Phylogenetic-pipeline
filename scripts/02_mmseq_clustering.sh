#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 input_fasta output_dir"
  exit 1
fi

IN_FASTA=$1
OUTDIR=$2

# Hyperparameters (override via env vars)
THREADS=${THREADS:-8}
COV=${COV:-0.8}

if ! command -v mmseqs >/dev/null 2>&1; then
  echo "Error: mmseqs not found in PATH"
  exit 2
fi

if [[ ! -s "$IN_FASTA" ]]; then
  echo "Error: input FASTA not found or empty: $IN_FASTA"
  exit 2
fi

if [[ -d "$OUTDIR" && -n "$(ls -A "$OUTDIR" 2>/dev/null)" ]]; then
  echo "Error: output directory is not empty: $OUTDIR"
  exit 3
fi

mkdir -p "$OUTDIR"

echo "Clustering"
mmseqs easy-cluster \
  "$IN_FASTA" \
  "$OUTDIR/clu" \
  "$OUTDIR/tmp" \
  --threads "$THREADS" \
  -c "$COV" \
  -v 1

echo "MMseqs2 done → $OUTDIR"
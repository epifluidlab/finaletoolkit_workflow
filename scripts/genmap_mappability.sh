#!/usr/bin/env bash
#
# genmap_mappability.sh — compute a CONTINUOUS k-mer mappability bigWig for any
# genome using GenMap, suitable for the workflow's `mappability_file`.
#
#   Usage:
#     scripts/genmap_mappability.sh <genome.fa> <chrom.sizes> <out.bw> [k=45] [e=0]
#
# The track value at each position is 1 / (number of occurrences of the k-mer),
# i.e. continuous in [0,1] (1 = uniquely mappable). The workflow filters bins by
# the per-bin MEAN of this track.
#
# Requirements: genmap, bedGraphToBigWig  (conda install -c bioconda genmap ucsc-bedgraphtobigwig)
# This is memory/time-heavy for large genomes.
#

set -euo pipefail

FA="${1:?genome fasta}"
CHROMSIZES="${2:?chrom.sizes}"
OUT="${3:?output .bw}"
K="${4:-45}"
E="${5:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-8}}"

command -v genmap          >/dev/null || { echo "ERROR: genmap not in PATH" >&2; exit 3; }
command -v bedGraphToBigWig >/dev/null || { echo "ERROR: bedGraphToBigWig not in PATH" >&2; exit 3; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/genmap.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
IDX="$WORK/index"
PREFIX="$WORK/map"

echo "[$(date)] genmap index ($FA) ..."
genmap index -F "$FA" -I "$IDX"

echo "[$(date)] genmap map -K $K -E $E (T=$THREADS) ..."
genmap map -K "$K" -E "$E" -I "$IDX" -O "$PREFIX" -bg -T "$THREADS"

echo "[$(date)] bedGraph -> bigWig ..."
LC_ALL=C sort -k1,1 -k2,2n "${PREFIX}.bedgraph" > "$WORK/map.sorted.bedgraph"
bedGraphToBigWig "$WORK/map.sorted.bedgraph" "$CHROMSIZES" "$OUT"

echo "[$(date)] done -> $OUT"

#!/usr/bin/env bash
#
# ftk_exc.sh — submit the workflow to SLURM via slurm_profile, in the background.
#
#   ./ftk_exc.sh [configfile] [tmpdir]
#
#   configfile : params YAML to run (default: params.yaml)
#   tmpdir     : scratch dir for per-job temporaries (default: ./tmp)
#
# Output/progress is written to snakemake.log.
#
set -euo pipefail

CONFIG="${1:-params.yaml}"
TMPDIR_PATH="${2:-./tmp}"
mkdir -p "$TMPDIR_PATH"

snakemake --profile slurm_profile \
  --configfile "$CONFIG" \
  --rerun-incomplete \
  --default-resources "tmpdir='$TMPDIR_PATH'" \
  > snakemake.log 2>&1 &

echo "submitted in background (PID $!) using $CONFIG"
echo "follow progress with:  tail -f snakemake.log"

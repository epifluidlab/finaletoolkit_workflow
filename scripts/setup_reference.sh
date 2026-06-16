#!/usr/bin/env bash
#
# setup_reference.sh — build the per-genome supplement/reference files that the
# FinaleToolkit workflow needs, into one directory, with the exact filenames the
# example configs (params.hg38.yaml / params.t2t-chm13.yaml) expect.
#
# Everything except the (multi-GB) mappability bigWig is built live from public
# sources (UCSC, Boyle-Lab Blacklist, BEDbase). The mappability track is fetched
# from Zenodo, supplied with --mappability, or built with scripts/genmap_mappability.sh.
#
#   Usage:
#     scripts/setup_reference.sh <hg38|t2t-chm13> [supplement_dir] [bin_size_kb]
#
#   Examples:
#     scripts/setup_reference.sh t2t-chm13 supplement 500
#     scripts/setup_reference.sh hg38      supplement 500
#
#   Options (env or flags):
#     --mappability <file.bw>   Use a local mappability bigWig instead of downloading.
#     --blacklist <file.bed>    Use a local blacklist BED (cleaned to BED3) instead of downloading.
#     --no-mappability          Skip the mappability file entirely.
#
# Produces (in <supplement_dir>):
#     <g>.chrom.sizes  <g>.2bit  <g>.<N>kb.bins  <g>.blacklist.bed
#     <g>.gap.bed      <g>.45mer.mappability.bw     (g = hg38 | chm13)
#
set -euo pipefail

# --------------------------------------------------------------------------- #
# Hosting / public sources
# --------------------------------------------------------------------------- #
# Continuous 45-mer mappability bigWigs (multi-GB) are hosted on Zenodo
# (FinaleToolkit Dataset, concept DOI 10.5281/zenodo.14284132). Override with the
# FTK_*_URL env vars if needed.
ZENODO_CHM13_MAPPABILITY_URL="${FTK_CHM13_MAPPABILITY_URL:-https://zenodo.org/api/records/20724659/files/chm13.45mer.mappability.bw/content}"
ZENODO_HG38_MAPPABILITY_URL="${FTK_HG38_MAPPABILITY_URL:-https://zenodo.org/api/records/20724659/files/hg38.45mer.mappability.bw/content}"

# Blacklists.
HG38_BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
# T2T-CHM13: excluderanges "T2T.excluderanges" (Boyle-Lab Blacklist: High-Signal +
# Low-Mappability) via the BEDbase mirror — the canonical Bioconductor AnnotationHub
# host (Azure blob) is unreachable from many HPC networks. BEDbase record id pinned:
#   5e6171d695d1deeaad4894658a8dc240  (T2T.excluderanges)
CHM13_BLACKLIST_URL="https://api.bedbase.org/v1/files/files/5/e/5e6171d695d1deeaad4894658a8dc240.bed.gz"

# --------------------------------------------------------------------------- #
# Args
# --------------------------------------------------------------------------- #
LOCAL_MAPPABILITY=""
LOCAL_BLACKLIST=""
SKIP_MAPPABILITY=0
POS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mappability)    LOCAL_MAPPABILITY="${2:?--mappability needs a path}"; shift 2 ;;
    --blacklist)      LOCAL_BLACKLIST="${2:?--blacklist needs a path}"; shift 2 ;;
    --no-mappability) SKIP_MAPPABILITY=1; shift ;;
    -h|--help)        grep '^#' "$0" | sed 's/^#//'; exit 0 ;;
    *)                POS+=("$1"); shift ;;
  esac
done
set -- "${POS[@]:-}"

GENOME="${1:-}"
SUPP_DIR="${2:-supplement}"
BIN_KB="${3:-500}"

case "$GENOME" in
  hg38)      G=hg38;  UCSC_DB=hg38; CHROM_RE='chr([1-9]|1[0-9]|2[0-2]|X)';   BLACKLIST_URL="$HG38_BLACKLIST_URL" ;;
  t2t-chm13) G=chm13; UCSC_DB=hs1;  CHROM_RE='chr([1-9]|1[0-9]|2[0-2]|X|M)'; BLACKLIST_URL="$CHM13_BLACKLIST_URL" ;;
  *) echo "ERROR: genome must be 'hg38' or 't2t-chm13' (got '${GENOME:-<none>}')" >&2; exit 2 ;;
esac

mkdir -p "$SUPP_DIR"
BIN_BP=$(( BIN_KB * 1000 ))
BINS="$SUPP_DIR/${G}.${BIN_KB}kb.bins"

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
require() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH (conda activate finaletoolkit_workflow)" >&2; exit 3; }; }
have()    { [[ -s "$1" ]]; }
fetch()   {
  local url="$1" out="$2"
  if have "$out"; then echo "  [skip] $(basename "$out") exists"; return 0; fi
  echo "  [get ] $(basename "$out")  <-  $url"
  curl -fSL --retry 3 -o "$out.part" "$url"
  mv "$out.part" "$out"
}
# stdin -> sorted BED3 (drops header/track lines and non-numeric coords)
clean_bed3() { awk 'BEGIN{OFS="\t"} $1!~/^(#|track|browser)/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2<$3 {print $1,$2,$3}' | LC_ALL=C sort -k1,1 -k2,2n; }

require curl
require bedtools

echo "== FinaleToolkit reference setup: $GENOME -> $SUPP_DIR (bins=${BIN_KB}kb) =="

# --------------------------------------------------------------------------- #
# 2bit + chrom.sizes (UCSC, both genomes)
# --------------------------------------------------------------------------- #
fetch "https://hgdownload.soe.ucsc.edu/goldenPath/${UCSC_DB}/bigZips/${UCSC_DB}.2bit"        "$SUPP_DIR/${G}.2bit"
fetch "https://hgdownload.soe.ucsc.edu/goldenPath/${UCSC_DB}/bigZips/${UCSC_DB}.chrom.sizes" "$SUPP_DIR/${G}.chrom.sizes"

# --------------------------------------------------------------------------- #
# Interval bins (autosomes + chrX, plus chrM for T2T-CHM13)
# --------------------------------------------------------------------------- #
if ! have "$BINS"; then
  echo "  [gen ] ${G}.${BIN_KB}kb.bins"
  grep -wE "$CHROM_RE" "$SUPP_DIR/${G}.chrom.sizes" > "$SUPP_DIR/.${G}.genome"
  bedtools makewindows -g "$SUPP_DIR/.${G}.genome" -w "$BIN_BP" > "$BINS"
  rm -f "$SUPP_DIR/.${G}.genome"
fi

# --------------------------------------------------------------------------- #
# Blacklist -> clean sorted BED3 (local override, else public source per genome)
# --------------------------------------------------------------------------- #
if ! have "$SUPP_DIR/${G}.blacklist.bed"; then
  if [[ -n "$LOCAL_BLACKLIST" ]]; then
    echo "  [copy] ${G}.blacklist.bed <- $LOCAL_BLACKLIST"
    clean_bed3 < "$LOCAL_BLACKLIST" > "$SUPP_DIR/${G}.blacklist.bed"
  else
    fetch "$BLACKLIST_URL" "$SUPP_DIR/${G}.blacklist.bed.gz"
    zcat "$SUPP_DIR/${G}.blacklist.bed.gz" | clean_bed3 > "$SUPP_DIR/${G}.blacklist.bed"
    rm -f "$SUPP_DIR/${G}.blacklist.bed.gz"
  fi
fi

# --------------------------------------------------------------------------- #
# Assembly gaps
# --------------------------------------------------------------------------- #
if [[ "$G" == "chm13" ]]; then
  echo "  [gen ] chm13.gap.bed (empty; T2T-CHM13 is gap-free)"
  : > "$SUPP_DIR/chm13.gap.bed"
elif ! have "$SUPP_DIR/hg38.gap.bed"; then
  require finaletoolkit
  echo "  [gen ] hg38.gap.bed (finaletoolkit gap-bed)"
  finaletoolkit gap-bed hg38 "$SUPP_DIR/hg38.gap.bed"
fi

# --------------------------------------------------------------------------- #
# Mappability bigWig
# --------------------------------------------------------------------------- #
MAPP="$SUPP_DIR/${G}.45mer.mappability.bw"
if [[ "$G" == "chm13" ]]; then MAPP_URL="$ZENODO_CHM13_MAPPABILITY_URL"; else MAPP_URL="$ZENODO_HG38_MAPPABILITY_URL"; fi

if [[ "$SKIP_MAPPABILITY" == "1" ]]; then
  echo "  [skip] mappability (--no-mappability)"
elif [[ -n "$LOCAL_MAPPABILITY" ]]; then
  echo "  [copy] mappability <- $LOCAL_MAPPABILITY"
  cp -f "$LOCAL_MAPPABILITY" "$MAPP"
elif [[ -n "$MAPP_URL" ]]; then
  fetch "$MAPP_URL" "$MAPP"   # continuous 45-mer track from Zenodo (multi-GB)
else
  cat >&2 <<EOF
  [note] Mappability track not fetched. Set the Zenodo URL for $G at the top of
         this script (ZENODO_${G^^}_MAPPABILITY_URL), or pass
         --mappability <${G}.45mer.mappability.bw>, or build it with GenMap:
           scripts/genmap_mappability.sh <${G}.fa> $SUPP_DIR/${G}.chrom.sizes \\
               $MAPP 45 0
EOF
fi

echo "== done. files in $SUPP_DIR: =="
ls -1 "$SUPP_DIR" | sed 's/^/   /'

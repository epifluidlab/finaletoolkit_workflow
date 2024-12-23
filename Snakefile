# MUST INSTALL:
# pip install finaletoolkit
# conda install bedtools htslib

# Run snakemake --configfile [path-to-yaml]
import os
import glob

OUTPUT_DIR = config.get("output_dir", "output")
INPUT_DIR = config.get("input_dir", "input")
MAPQ = int(config.get("mapq", 0))
BLACKLIST = config.get("blacklist", None)
COMPRESS = config.get("compress", False)

# Finale Toolkit parameters
FRAG_LENGTH_BINS = config.get("frag_length_bins", False)
FRAG_LENGTH_BINS_BIN_SIZE = int(config.get("frag_length_bins_bin_size", 10)) if FRAG_LENGTH_BINS else None
FRAG_LENGTH_BINS_MIN = int(config.get("frag_length_bins_min", 50)) if FRAG_LENGTH_BINS else None
FRAG_LENGTH_BINS_MAX = int(config.get("frag_length_bins_max", 350)) if FRAG_LENGTH_BINS else None
FRAG_LENGTH_BINS_CHROM = config.get("frag_length_bins_chrom", None) if FRAG_LENGTH_BINS else None
FRAG_LENGTH_BINS_START = int(config.get("frag_length_bins_start", 0)) if FRAG_LENGTH_BINS else None
FRAG_LENGTH_BINS_END = int(config.get("frag_length_bins_end", 0)) if FRAG_LENGTH_BINS else None
FRAG_LENGTH_BINS_OUTPUT_TSV = config.get("frag_length_bins_output_tsv", None) if FRAG_LENGTH_BINS else None
FRAG_LENGTH_BINS_HISTOGRAM = config.get("frag_length_bins_histogram", None) if FRAG_LENGTH_BINS else None

FRAG_LENGTH_INTERVALS = config.get("frag_length_intervals", False)
FRAG_LENGTH_INTERVALS_MIN = int(config.get("frag_length_intervals_min", 50)) if FRAG_LENGTH_INTERVALS else None
FRAG_LENGTH_INTERVALS_MAX = int(config.get("frag_length_intervals_max", 350)) if FRAG_LENGTH_INTERVALS else None
FRAG_LENGTH_INTERVALS_OUTPUT_BED = config.get("frag_length_intervals_output_bed", None) if FRAG_LENGTH_INTERVALS else None
FRAG_LENGTH_INTERVALS_WORKERS = int(config.get("frag_length_intervals_workers", 1)) if FRAG_LENGTH_INTERVALS else None

if not OUTPUT_DIR:
    raise ValueError("output_dir must be specified in the configfile")
if not INPUT_DIR:
    raise ValueError("input_dir must be specified in the configfile")

bed_files = glob.glob(os.path.join(INPUT_DIR, "*.bed.gz"))
bam_files = glob.glob(os.path.join(INPUT_DIR, "*.bam"))
cram_files = glob.glob(os.path.join(INPUT_DIR, "*.cram"))
input_files = bed_files + bam_files + cram_files
SAMPLES = [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in input_files]

rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "{sample}.filtered.bed.gz" if COMPRESS else "{sample}.filtered.bed"), sample=SAMPLES),
        expand(FRAG_LENGTH_BINS_OUTPUT_TSV if FRAG_LENGTH_BINS else [], sample=SAMPLES),
        expand(FRAG_LENGTH_BINS_HISTOGRAM if FRAG_LENGTH_BINS else [], sample=SAMPLES),
        expand(FRAG_LENGTH_INTERVALS_OUTPUT_BED if FRAG_LENGTH_INTERVALS else [], sample=SAMPLES),

rule filter_bed:
    input:
        os.path.join(INPUT_DIR, "{sample}.bed.gz")
    output:
        os.path.join(OUTPUT_DIR, "{sample}.tmp.bed")
    run:
        if MAPQ > 0:
            shell(
                """
                zcat {input} | awk -F '\\t' '$4 >= {MAPQ} {{print $1 "\\t" $2 "\\t" $3 "\\t" $4}}' > {output}
                """
            )
        else:
            shell(f"zcat {input} > {output}")

rule blacklist_filter:
    input:
        bed=os.path.join(OUTPUT_DIR, "{sample}.tmp.bed"),
    output:
        os.path.join(OUTPUT_DIR, "{sample}.filtered.bed")
    run:
        if BLACKLIST:
            shell(f"bedtools intersect -v -a {input.bed} -b {BLACKLIST} > {output}")
        else:
            shell(f"cp {input.bed} {output}")

if COMPRESS:
    rule compress_bed:
        input:
            os.path.join(OUTPUT_DIR, "{sample}.filtered.bed")
        output:
            os.path.join(OUTPUT_DIR, "{sample}.filtered.bed.gz")
        run:
            shell(f"bgzip -c {input} > {output}")

# Finale Toolkit Rules
if FRAG_LENGTH_BINS: #only make rules if the main switch is true
    rule frag_length_bins:
        input:
            bam_or_cram = lambda wildcards: [f for f in bam_files + cram_files if wildcards.sample in os.path.basename(f)][0]
        output:
            tsv=FRAG_LENGTH_BINS_OUTPUT_TSV,
            histogram=FRAG_LENGTH_BINS_HISTOGRAM
        params:
            bin_size=FRAG_LENGTH_BINS_BIN_SIZE,
            min_len=FRAG_LENGTH_BINS_MIN,
            max_len=FRAG_LENGTH_BINS_MAX,
            chrom=FRAG_LENGTH_BINS_CHROM,
            start=FRAG_LENGTH_BINS_START,
            end=FRAG_LENGTH_BINS_END
        shell:
            """
            finaletoolkit frag-length-bins {input.bam_or_cram} \
              -q {MAPQ} \
              --bin-size {params.bin_size} \
              -p midpoint \
              -min {params.min_len} \
              -max {params.max_len} \
              -c {params.chrom} \
              -S {params.start} \
              -E {params.end} \
              -o {output.tsv} \
              --histogram-path {output.histogram} \
              -v
            """

if FRAG_LENGTH_INTERVALS: #only make rules if the main switch is true
    rule frag_length_intervals:
        input:
            bam_or_cram = lambda wildcards: [f for f in bam_files + cram_files if wildcards.sample in os.path.basename(f)][0]
        output:
            bed=FRAG_LENGTH_INTERVALS_OUTPUT_BED
        params:
            min_len=FRAG_LENGTH_INTERVALS_MIN,
            max_len=FRAG_LENGTH_INTERVALS_MAX
        shell:
            """
            finaletoolkit frag-length-intervals {input.bam_or_cram} {output.bed}.bins \
              -q {MAPQ} \
              -p midpoint \
              -min {params.min_len} \
              -max {params.max_len} \
              -o {output.bed} \
              -w {FRAG_LENGTH_INTERVALS_WORKERS} \
              -v
            """

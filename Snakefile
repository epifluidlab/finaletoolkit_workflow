import os
import glob
import subprocess

def check_tools():
    tools = {
    "finaletoolkit": "finaletoolkit --help",
    "bedtools": "bedtools --help",
    "htslib": "htsfile --help",
    "samtools": "samtools --help"
    }
    missing_tools = []
    for tool, command in tools.items():
        try:
            subprocess.run(command, shell=True, check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing_tools.append(tool)
    if missing_tools:
        raise SystemExit(f"Error: The following tools are not installed: {', '.join(missing_tools)}. Please install them. (finaletoolkit via pip, bedtools, samtools, htslib via conda or system package manager)")

check_tools() #run the check

OUTPUT_DIR = config.get("output_dir", "output")

INPUT_DIR = config.get("input_dir", "input")

# Should contain interval, blacklist, .chrom.sizes, .2bit, and other files that supplement the input files
SUPPLEMENT_DIR = config.get("supplement_dir", "supplement")

MAPQ = int(config.get("mapq", 0))

BLACKLIST = config.get("blacklist", None)

# Finale Toolkit parameters
FINALETOOLKIT = config.get("use_finaletoolkit", False) # Specify whether you are using one of Finaletoolkit's features
frag_length_bins_mapq = config.get("frag_length_bins_mapq", -1)
frag_length_bins_bin_size = config.get("frag_length_bins_bin_size", -1)
frag_length_bins_min_len = config.get("frag_length_bins_min_len", -1)
frag_length_bins_max_len = config.get("frag_length_bins_max_len", -1)
frag_length_bins_chrom = config.get("frag_length_bins_chrom", "-1")
frag_length_bins_start = config.get("frag_length_bins_start", -1)
frag_length_bins_end = config.get("frag_length_bins_end", -1)

bed_files = glob.glob(os.path.join(INPUT_DIR, "*.bed.gz"))
bam_files = glob.glob(os.path.join(INPUT_DIR, "*.bam"))
cram_files = glob.glob(os.path.join(INPUT_DIR, "*.cram"))

SAMPLES_BED = [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in bed_files]
SAMPLES_BAM = [os.path.splitext(os.path.basename(f))[0] for f in bam_files]
SAMPLES_CRAM = [os.path.splitext(os.path.basename(f))[0] for f in cram_files]
SAMPLES_ALL = SAMPLES_BED + SAMPLES_BAM + SAMPLES_CRAM

rule all:
    input:
        # BED output
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.bed.gz") if not FINALETOOLKIT else [], sample=SAMPLES_BED),
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.bed.gz.tbi") if not FINALETOOLKIT else [], sample=SAMPLES_BED),
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.bed.gz.tsv"), sample=SAMPLES_BED, ext=["bam", "bed.gz", "cram"]) if FINALETOOLKIT else [],
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.bed.gz.png"), sample=SAMPLES_BED, ext=["bam", "bed.gz", "cram"]) if FINALETOOLKIT else [],

        # BAM output
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.bam") if not FINALETOOLKIT else [], sample=SAMPLES_BAM),
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.bam.bai") if not FINALETOOLKIT else [], sample=SAMPLES_BAM),
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.bam.tsv"), sample=SAMPLES_BAM, ext=["bam", "bed.gz", "cram"]) if FINALETOOLKIT else [],
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.bam.png"), sample=SAMPLES_BAM, ext=["bam", "bed.gz", "cram"]) if FINALETOOLKIT else [],

        # CRAM output
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.cram") if not FINALETOOLKIT else [], sample=SAMPLES_CRAM),
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.cram.crai") if not FINALETOOLKIT else [], sample=SAMPLES_CRAM),
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.cram.tsv"), sample=SAMPLES_CRAM, ext=["bam", "bed.gz", "cram"]) if FINALETOOLKIT else [],
        expand(os.path.join(OUTPUT_DIR, "{sample}.final.cram.png"), sample=SAMPLES_CRAM, ext=["bam", "bed.gz", "cram"]) if FINALETOOLKIT else [],

# STEP 0: Dramatic performance improvements by compressing and indexing the blacklist file
if BLACKLIST:
    blacklist_input_path = os.path.join(SUPPLEMENT_DIR, BLACKLIST)

    if not os.path.exists(blacklist_input_path):
        raise FileNotFoundError(f"Blacklist file not found: {blacklist_input_path}. The blacklist file must be in the supplement directory.")

    if BLACKLIST.endswith(".gz"):
        blacklist_output_gz_path = blacklist_input_path
        blacklist_output_tbi_path = blacklist_output_gz_path + ".tbi"
        if not os.path.exists(blacklist_output_tbi_path):
            shell(f"tabix -p bed '{blacklist_output_gz_path}'")
    else:
        blacklist_output_gz_path = os.path.join(SUPPLEMENT_DIR, BLACKLIST + ".gz")
        blacklist_output_tbi_path = blacklist_output_gz_path + ".tbi"
        if not os.path.exists(blacklist_output_gz_path):
            print(blacklist_input_path,blacklist_output_gz_path)
            shell(f"bgzip -@ 4 -c {blacklist_input_path} > {blacklist_output_gz_path}")
        if not os.path.exists(blacklist_output_tbi_path):
            shell(f"tabix -p bed '{blacklist_output_gz_path}'")
    BLACKLIST = blacklist_output_gz_path

# STEP 1: Run filtering based on MAPQ scores
rule filter_bed_mapq:
    input:
        raw_bed=os.path.join(INPUT_DIR, "{sample}.bed.gz"),
        raw_tbi=os.path.join(INPUT_DIR, "{sample}.bed.gz.tbi")
    output:
        bed=os.path.join(OUTPUT_DIR, "{sample}.filtered.bed.gz"),
        tbi=os.path.join(OUTPUT_DIR, "{sample}.filtered.bed.gz.tbi")
    threads: 4 
    run:
        if MAPQ > 0:
            shell(f"zcat {input.raw_bed} | awk -F '\\t' '$4 >= {MAPQ}' | bgzip -c -@ {threads} > {output.bed}")
            shell(f"tabix -p bed {output.bed}")
        else:
            shell(f"cp {input.raw_bed} {output.bed}")
            shell(f"cp {input.raw_tbi} {output.tbi}")

rule filter_bam_mapq:
    input:
        raw_bam=os.path.join(INPUT_DIR, "{sample}.bam"),
        raw_bai=os.path.join(INPUT_DIR, "{sample}.bam.bai")
    output:
        bam=os.path.join(OUTPUT_DIR, "{sample}.filtered.bam"),
        bai=os.path.join(OUTPUT_DIR, "{sample}.filtered.bam.bai")
    threads: 4
    run:
        if MAPQ > 0:
            shell(f"samtools view -b -q {MAPQ} -@ {threads} {input.raw_bam} > {output.bam}")
            shell(f"samtools index {output.bam}")
        else:
            shell(f"cp {input.raw_bam} {output.bam}")
            shell(f"cp {input.raw_bai} {output.bai}")

rule filter_cram_mapq:
    input:
        raw_cram=os.path.join(INPUT_DIR, "{sample}.cram"),
        raw_crai=os.path.join(INPUT_DIR, "{sample}.cram.crai")
    output:
        cram=os.path.join(OUTPUT_DIR, "{sample}.filtered.cram"),
        crai=os.path.join(OUTPUT_DIR, "{sample}.filtered.cram.crai")
    threads: 4
    run:
        if MAPQ > 0:
            shell(f"samtools view -b -q {MAPQ} -@ {threads} {input.raw_cram} > {output.cram}")
            shell(f"samtools index {output.cram}")
        else:
            shell(f"cp {input.raw_cram} {output.cram}")
            shell(f"cp {input.raw_crai} {output.crai}")

# STEP 2: Remove regions based off of blacklist file
rule blacklist_bed:
    input:
        bed=os.path.join(OUTPUT_DIR, "{sample}.filtered.bed.gz"),
        tbi=os.path.join(OUTPUT_DIR, "{sample}.filtered.bed.gz.tbi"),
    output:
        bed=os.path.join(OUTPUT_DIR, "{sample}.final.bed.gz"),
        tbi=os.path.join(OUTPUT_DIR, "{sample}.final.bed.gz.tbi"),
    params:
        blacklist=BLACKLIST
    threads: 4
    shell:
        """
        if [[ "{params.blacklist}" == "None" ]]; then
            cp {input.bed} {output.bed}
            cp {input.tbi} {output.tbi}
        else
            bedtools subtract -a {input.bed} -b {params.blacklist} | bgzip -@ {threads} > {output.bed}
            tabix -p bed {output.bed}
        fi
        """

rule blacklist_bam:
    input:
        bam=os.path.join(OUTPUT_DIR, "{sample}.filtered.bam"),
        bai=os.path.join(OUTPUT_DIR, "{sample}.filtered.bam.bai"),
    output:
        bam=os.path.join(OUTPUT_DIR, "{sample}.final.bam"),
        bai=os.path.join(OUTPUT_DIR, "{sample}.final.bam.bai"),
    params:
        blacklist=BLACKLIST
    threads: 4
    shell:
        """
        if [[ "{params.blacklist}" == "None" ]]; then
            cp {input.bam} {output.bam}
            cp {input.bai} {output.bai}
        else
            samtools view -b -L {params.blacklist} -@ {threads} {input.bam} > {output.bam}
            samtools index {output.bam}
        fi
        """

rule blacklist_cram:
    input:
        cram=os.path.join(OUTPUT_DIR, "{sample}.filtered.cram"),
        crai=os.path.join(OUTPUT_DIR, "{sample}.filtered.cram.crai"),
    output:
        cram=os.path.join(OUTPUT_DIR, "{sample}.final.cram"),
        crai=os.path.join(OUTPUT_DIR, "{sample}.final.cram.crai"),
    params:
        blacklist=BLACKLIST
    threads: 4 
    shell:
        """
        if [[ "{params.blacklist}" == "None" ]]; then
            cp {input.cram} {output.cram}
            cp {input.crai} {output.crai}
        else
            samtools view -b -L {params.blacklist} -@ {threads} {input.cram} > {output.cram}
            samtools index {output.cram}
        fi
        """

# STEP 3: Finaletoolkit commands
rule finaletoolkit_frag_length_bins:
    input:
        lambda wildcards: os.path.join(OUTPUT_DIR, f"{wildcards.sample}.final.{wildcards.ext}")
    output:
        tsv=os.path.join(OUTPUT_DIR, "{sample}.final.{ext}.tsv"),
        png=os.path.join(OUTPUT_DIR, "{sample}.final.{ext}.png")
    run:
        command_parts = [f"finaletoolkit frag-length-bins {input}"]
        if frag_length_bins_mapq != -1:
            command_parts.append(f"-q {frag_length_bins_mapq}")
        if frag_length_bins_bin_size != -1:
            command_parts.append(f"--bin-size {frag_length_bins_bin_size}")
        command_parts.append("-p midpoint")
        if frag_length_bins_min_len != -1:
            command_parts.append(f"-min {frag_length_bins_min_len}")
        if frag_length_bins_max_len != -1:
            command_parts.append(f"-max {frag_length_bins_max_len}")
        if frag_length_bins_chrom != "-1":
            command_parts.append(f"-c {frag_length_bins_chrom}")
        if frag_length_bins_start != -1:
            command_parts.append(f"-S {frag_length_bins_start}")
        if frag_length_bins_end != -1:
            command_parts.append(f"-E {frag_length_bins_end}")
        command_parts.append(f"-o {output.tsv}")
        command_parts.append(f"--histogram-path {output.png}")
        command_parts.append("-v")
        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")

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

output_dir = config.get("output_dir", "output")

input_dir = config.get("input_dir", "input")

# Should contain interval, blacklist, .chrom.sizes, .2bit, and other files that supplement the input files
SUPPLEMENT_DIR = config.get("supplement_dir", "supplement")

MAPQ = int(config.get("mapq", 0))

BLACKLIST = config.get("blacklist", None)

# Finale Toolkit parameters

# frag-length-bins
frag_length_bins = config.get("frag_length_bins", False)
frag_length_bins_mapq = config.get("frag_length_bins_mapq", None)
frag_length_bins_policy = config.get("frag_length_bins_policy", None)
frag_length_bins_bin_size = config.get("frag_length_bins_bin_size", None)
frag_length_bins_min_len = config.get("frag_length_bins_min_len", None)
frag_length_bins_max_len = config.get("frag_length_bins_max_len", None)
frag_length_bins_chrom = config.get("frag_length_bins_chrom", "")
frag_length_bins_start = config.get("frag_length_bins_start", None)
frag_length_bins_end = config.get("frag_length_bins_end", None)

# frag-length-intervals
frag_length_intervals = config.get("frag_length_intervals", False)
frag_length_interval_file = config.get("frag_length_interval_file", False)
frag_length_intervals_mapq = config.get("frag_length_intervals_mapq", None)
frag_length_intervals_policy = config.get("frag_length_intervals_policy", None)
frag_length_intervals_min_len = config.get("frag_length_intervals_min_len", None)
frag_length_intervals_max_len = config.get("frag_length_intervals_max_len", None)
frag_length_intervals_workers = config.get("frag_length_intervals_workers", 1)

# coverage
coverage = config.get("coverage", False)
coverage_interval_file = config.get("coverage_interval_file", False)
coverage_mapq = config.get("coverage_mapq", None)
coverage_intersect_policy = config.get("coverage_intersect_policy", None)
coverage_min_len = config.get("coverage_min_len", None)
coverage_max_len = config.get("coverage_max_len", None)
coverage_workers = config.get("coverage_workers", 1)
coverage_normalize = config.get("coverage_normalize", False)
coverage_scale_factor = config.get("coverage_scale_factor", None)

# end-motifs
end_motifs = config.get("end_motifs", False)
end_motifs_refseq_file = config.get("end_motifs_refseq_file", None)
end_motifs_kmer_length = config.get("end_motifs_kmer_length", 4)
end_motifs_min_len = config.get("end_motifs_min_len", None)
end_motifs_max_len = config.get("end_motifs_max_len", None)
end_motifs_no_both_strands = config.get("end_motifs_no_both_strands", True)
end_motifs_negative_strand = config.get("end_motifs_negative_strand", False)
end_motifs_mapq = config.get("end_motifs_mapq", 20)
end_motifs_workers = config.get("end_motifs_workers", 1)

# interval-end-motifs
interval_end_motifs = config.get("interval_end_motifs", False)
interval_end_motifs_refseq_file = config.get("interval_end_motifs_refseq_file", None)
interval_end_motifs_interval_file = config.get("interval_end_motifs_interval_file", None)
interval_end_motifs_kmer_length = config.get("interval_end_motifs_kmer_length", 4)
interval_end_motifs_min_len = config.get("interval_end_motifs_min_len", None)
interval_end_motifs_max_len = config.get("interval_end_motifs_max_len", None)
interval_end_motifs_no_both_strands = config.get("interval_end_motifs_no_both_strands", True)
interval_end_motifs_negative_strand = config.get("interval_end_motifs_negative_strand", False)
interval_end_motifs_mapq = config.get("interval_end_motifs_mapq", 20)
interval_end_motifs_workers = config.get("interval_end_motifs_workers", None)

using_finaletoolkit = frag_length_bins or frag_length_intervals or coverage or end_motifs

bed_files = glob.glob(os.path.join(input_dir, "*.bed.gz"))
bam_files = glob.glob(os.path.join(input_dir, "*.bam"))
cram_files = glob.glob(os.path.join(input_dir, "*.cram"))

sample_bed = [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in bed_files]
sample_bam = [os.path.splitext(os.path.basename(f))[0] for f in bam_files]
sample_cram = [os.path.splitext(os.path.basename(f))[0] for f in cram_files]

def shortcut(endings, samples, condition):
    return expand(os.path.join(output_dir, "{sample}.{ending}"), sample=samples, ending=endings) if condition else []

rule all:
    input:
        # BED output
        shortcut(["final.bed.gz", "final.bed.gz.tbi"], sample_bed, not using_finaletoolkit),
        # BAM output
        shortcut(["final.bam", "final.bam.bai"], sample_bam, not using_finaletoolkit),
        # CRAM output
        shortcut(["final.cram", "final.cram.crai"], sample_cram, not using_finaletoolkit),

        # Finale Toolkit outputs (separate bins and intervals)

        # (BED)
        shortcut(["bed.gz.frag_length_bins.tsv", "bed.gz.frag_length_bins.png"], sample_bed, frag_length_bins),
        shortcut(["bed.gz.frag_length_intervals.bed"], sample_bed, frag_length_intervals),
        shortcut(["bed.gz.coverage.bed"], sample_bed, coverage),
        shortcut(["bed.gz.end_motifs.tsv"], sample_bed, end_motifs),
        shortcut(["bed.gz.interval_end_motifs.tsv"], sample_bed, interval_end_motifs),

        # (BAM)
        shortcut(["bam.frag_length_bins.tsv","bam.frag_length_bins.png"], sample_bam, frag_length_bins),
        shortcut(["bam.frag_length_intervals.bed"], sample_bam, frag_length_intervals),
        shortcut(["bam.coverage.bed"], sample_bam, coverage),
        shortcut(["bam.end_motifs.tsv"], sample_bam, end_motifs),
        shortcut(["bam.interval_end_motifs.tsv"], sample_bam, interval_end_motifs),

        # (CRAM)
        shortcut(["cram.frag_length_bins.tsv","cram.frag_length_bins.png"], sample_cram, frag_length_bins),
        shortcut(["cram.frag_length_intervals.bed"], sample_cram, frag_length_intervals),
        shortcut(["cram.coverage.bed"], sample_cram, coverage),
        shortcut(["cram.end_motifs.tsv"], sample_cram, end_motifs),
        shortcut(["cram.interval_end_motifs.tsv"], sample_cram, interval_end_motifs),

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
        raw_bed=os.path.join(input_dir, "{sample}.bed.gz"),
        raw_tbi=os.path.join(input_dir, "{sample}.bed.gz.tbi")
    output:
        bed=os.path.join(output_dir, "{sample}.filtered.bed.gz"),
        tbi=os.path.join(output_dir, "{sample}.filtered.bed.gz.tbi")
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
        raw_bam=os.path.join(input_dir, "{sample}.bam"),
        raw_bai=os.path.join(input_dir, "{sample}.bam.bai")
    output:
        bam=os.path.join(output_dir, "{sample}.filtered.bam"),
        bai=os.path.join(output_dir, "{sample}.filtered.bam.bai")
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
        raw_cram=os.path.join(input_dir, "{sample}.cram"),
        raw_crai=os.path.join(input_dir, "{sample}.cram.crai")
    output:
        cram=os.path.join(output_dir, "{sample}.filtered.cram"),
        crai=os.path.join(output_dir, "{sample}.filtered.cram.crai")
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
        bed=os.path.join(output_dir, "{sample}.filtered.bed.gz"),
        tbi=os.path.join(output_dir, "{sample}.filtered.bed.gz.tbi"),
    output:
        bed=os.path.join(output_dir, "{sample}.final.bed.gz"),
        tbi=os.path.join(output_dir, "{sample}.final.bed.gz.tbi"),
    params:
        blacklist=BLACKLIST
    threads: 4
    shell:
        """
        if [[ "{params.blacklist}" == "None" ]]; then
            mv {input.bed} {output.bed}
            mv {input.tbi} {output.tbi}
        else
            bedtools subtract -a {input.bed} -b {params.blacklist} | bgzip -@ {threads} > {output.bed}
            tabix -p bed {output.bed}
            rm {input.bed}
            rm {input.tbi}
        fi
        """

rule blacklist_bam:
    input:
        bam=os.path.join(output_dir, "{sample}.filtered.bam"),
        bai=os.path.join(output_dir, "{sample}.filtered.bam.bai"),
    output:
        bam=os.path.join(output_dir, "{sample}.final.bam"),
        bai=os.path.join(output_dir, "{sample}.final.bam.bai"),
    params:
        blacklist=BLACKLIST
    threads: 4
    shell:
        """
        if [[ "{params.blacklist}" == "None" ]]; then
            mv {input.bam} {output.bam}
            mv {input.bai} {output.bai}
        else
            samtools view -b -L {params.blacklist} -@ {threads} {input.bam} > {output.bam}
            samtools index {output.bam}
            rm {input.bam}
            rm {input.bai}
        fi
        """

rule blacklist_cram:
    input:
        cram=os.path.join(output_dir, "{sample}.filtered.cram"),
        crai=os.path.join(output_dir, "{sample}.filtered.cram.crai"),
    output:
        cram=os.path.join(output_dir, "{sample}.final.cram"),
        crai=os.path.join(output_dir, "{sample}.final.cram.crai"),
    params:
        blacklist=BLACKLIST
    threads: 4 
    shell:
        """
        if [[ "{params.blacklist}" == "None" ]]; then
            mv {input.cram} {output.cram}
            mv {input.crai} {output.crai}
        else
            samtools view -b -L {params.blacklist} -@ {threads} {input.cram} > {output.cram}
            samtools index {output.cram}
            rm {input.cram}
            rm {input.crai}
        fi
        """

# STEP 3: Finaletoolkit commands. Using "is not None" to avoid falsy values accidentaly bypassing flags

rule frag_length_bins:
    input:
        lambda wildcards: os.path.join(output_dir, f"{wildcards.sample}.final.{wildcards.ext}")
    output:
        tsv=os.path.join(output_dir, "{sample}.{ext}.frag_length_bins.tsv"),
        png=os.path.join(output_dir, "{sample}.{ext}.frag_length_bins.png")
    run:
        command_parts = [f"finaletoolkit frag-length-bins {input}"]
        if frag_length_bins_mapq is not None:
            command_parts.append(f"-q {frag_length_bins_mapq}")
        if frag_length_bins_bin_size is not None:
            command_parts.append(f"--bin-size {frag_length_bins_bin_size}")
        if frag_length_bins_policy is not None:
            command_parts.append(f"-p {frag_length_bins_policy}")
        if frag_length_bins_min_len is not None:
            command_parts.append(f"-min {frag_length_bins_min_len}")
        if frag_length_bins_max_len is not None:
            command_parts.append(f"-max {frag_length_bins_max_len}")
        if frag_length_bins_chrom is not None and frag_length_bins_chrom != "":
            command_parts.append(f"-c {frag_length_bins_chrom}")
        if frag_length_bins_start is not None:
            command_parts.append(f"-S {frag_length_bins_start}")
        if frag_length_bins_end is not None:
            command_parts.append(f"-E {frag_length_bins_end}")
        command_parts.append(f"-o {output.tsv}")
        command_parts.append(f"--histogram-path {output.png} -v")
        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")

rule frag_length_intervals:
    input:
        data=lambda wildcards:  os.path.join(output_dir, f"{wildcards.sample}.final.{wildcards.ext}"),
        intervals=os.path.join(SUPPLEMENT_DIR, f"{frag_length_interval_file}")
    output:
        os.path.join(output_dir, "{sample}.{ext}.frag_length_intervals.bed")
    run:
        command_parts = [
            "finaletoolkit",
            "frag-length-intervals",
            input.data,
            input.intervals,
        ]

        if frag_length_intervals_min_len is not None:
            command_parts.append(f"-min {frag_length_intervals_min_len}")
        if frag_length_intervals_max_len is not None:
            command_parts.append(f"-max {frag_length_intervals_max_len}")

        if frag_length_intervals_policy is not None:
            command_parts.append(f"-p {frag_length_intervals_policy}")

        if frag_length_intervals_mapq is not None:
            command_parts.append(f"-q {frag_length_intervals_mapq}")

        if frag_length_intervals_workers is not None:
            command_parts.append(f"-w {frag_length_intervals_workers}")

        command_parts.append(f"-o {output} -v")

        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")

rule coverage:
    input:
        data=lambda wildcards:  os.path.join(output_dir, f"{wildcards.sample}.final.{wildcards.ext}"),
        intervals=os.path.join(SUPPLEMENT_DIR, f"{coverage_interval_file}")
    output:
        os.path.join(output_dir, "{sample}.{ext}.coverage.bed")
    run:
        command_parts = [
            "finaletoolkit",
            "coverage",
            input.data,
            input.intervals,
        ]

        if coverage_min_len is not None:
            command_parts.append(f"-min {coverage_min_len}")
        if coverage_max_len is not None:
            command_parts.append(f"-max {coverage_max_len}")

        if coverage_normalize:
            command_parts.append("-n")
            if coverage_scale_factor is not None:
                command_parts.append(f"-s {coverage_scale_factor}")

        if coverage_intersect_policy is not None:
            command_parts.append(f"-p {coverage_intersect_policy}")

        if coverage_mapq is not None:
            command_parts.append(f"-q {coverage_mapq}")

        if coverage_workers is not None:
            command_parts.append(f"-w {coverage_workers}")

        command_parts.append(f"-o {output} -v")

        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")

rule finaletoolkit_end_motifs:
    input:
        data=lambda wildcards: os.path.join(output_dir, f"{wildcards.sample}.final.{wildcards.ext}"),
        refseq=os.path.join(SUPPLEMENT_DIR, f"{end_motifs_refseq_file}")
    output:
        os.path.join(output_dir, "{sample}.{ext}.end_motifs.tsv")
    run:
        command_parts = [
            "finaletoolkit",
            "end-motifs",
            input.data,
            input.refseq,
        ]

        if end_motifs_kmer_length is not None:
            command_parts.append(f"-k {end_motifs_kmer_length}")
        if end_motifs_min_len is not None:
            command_parts.append(f"-min {end_motifs_min_len}")
        if end_motifs_max_len is not None:
            command_parts.append(f"-max {end_motifs_max_len}")

        if end_motifs_no_both_strands:
            command_parts.append("-B")
            if end_motifs_negative_strand:
                command_parts.append("-n")

        if end_motifs_mapq is not None:
            command_parts.append(f"-q {end_motifs_mapq}")
        if end_motifs_workers is not None:
            command_parts.append(f"-w {end_motifs_workers}")

        command_parts.append("-o")
        command_parts.append(f"{output}")
        command_parts.append("-v")

        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")

rule finaletoolkit_interval_end_motifs:
    input:
        data=lambda wildcards:  os.path.join(output_dir, f"{wildcards.sample}.final.{wildcards.ext}"),
        refseq=os.path.join(SUPPLEMENT_DIR, interval_end_motifs_refseq_file),
        intervals=os.path.join(SUPPLEMENT_DIR, interval_end_motifs_interval_file)
    output:
        os.path.join(output_dir, "{sample}.{ext}.interval_end_motifs.tsv")
    run:
        command_parts = [
            "finaletoolkit",
            "interval-end-motifs",
            input.data,
            input.refseq,
            input.intervals,
        ]

        if interval_end_motifs_kmer_length is not None:
            command_parts.append(f"-k {interval_end_motifs_kmer_length}")
        if interval_end_motifs_min_len is not None:
            command_parts.append(f"-min {interval_end_motifs_min_len}")
        if interval_end_motifs_max_len is not None:
            command_parts.append(f"-max {interval_end_motifs_max_len}")

        if interval_end_motifs_no_both_strands:
            command_parts.append("-B")
            if interval_end_motifs_negative_strand:
                command_parts.append("-n")

        if interval_end_motifs_mapq is not None:
            command_parts.append(f"-q {interval_end_motifs_mapq}")
        if interval_end_motifs_workers is not None:
            command_parts.append(f"-w {interval_end_motifs_workers}")

        command_parts.append("-o")
        command_parts.append(f"{output}")
        command_parts.append("-v")

        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")

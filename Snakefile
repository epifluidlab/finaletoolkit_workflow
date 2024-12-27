import os
import glob
import subprocess

cmb = os.path.join
cnfg = config.get

def io(endings, samples, condition):
    return expand(cmb(out_dir, "{sample}.{ending}"), sample=samples, ending=endings) if condition else []

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

check_tools()

out_dir = cnfg("output_dir", "output")
in_dir = cnfg("input_dir", "input")
sup_dir = cnfg("supplement_dir", "supplement") # Should contain interval, blacklist, .chrom.sizes, .2bit, and other files that supplement the input files

mapq = int(cnfg("mapq", 0))

blacklist = cnfg("blacklist", None)

# Finale Toolkit parameters

# frag-length-bins
frag_length_bins = cnfg("frag_length_bins", False)
frag_length_bins_mapq = cnfg("frag_length_bins_mapq", None)
frag_length_bins_policy = cnfg("frag_length_bins_policy", None)
frag_length_bins_bin_size = cnfg("frag_length_bins_bin_size", None)
frag_length_bins_min_len = cnfg("frag_length_bins_min_len", None)
frag_length_bins_max_len = cnfg("frag_length_bins_max_len", None)
frag_length_bins_chrom = cnfg("frag_length_bins_chrom", "")
frag_length_bins_start = cnfg("frag_length_bins_start", None)
frag_length_bins_end = cnfg("frag_length_bins_end", None)

# frag-length-intervals
frag_length_intervals = cnfg("frag_length_intervals", False)
frag_length_interval_file = cnfg("frag_length_interval_file", False)
frag_length_intervals_mapq = cnfg("frag_length_intervals_mapq", None)
frag_length_intervals_policy = cnfg("frag_length_intervals_policy", None)
frag_length_intervals_min_len = cnfg("frag_length_intervals_min_len", None)
frag_length_intervals_max_len = cnfg("frag_length_intervals_max_len", None)
frag_length_intervals_workers = cnfg("frag_length_intervals_workers", 1)

# coverage
coverage = cnfg("coverage", False)
coverage_interval_file = cnfg("coverage_interval_file", False)
coverage_mapq = cnfg("coverage_mapq", None)
coverage_intersect_policy = cnfg("coverage_intersect_policy", None)
coverage_min_len = cnfg("coverage_min_len", None)
coverage_max_len = cnfg("coverage_max_len", None)
coverage_workers = cnfg("coverage_workers", 1)
coverage_normalize = cnfg("coverage_normalize", False)
coverage_scale_factor = cnfg("coverage_scale_factor", None)

# end-motifs
end_motifs = cnfg("end_motifs", False)
end_motifs_refseq_file = cnfg("end_motifs_refseq_file", None)
end_motifs_kmer_length = cnfg("end_motifs_kmer_length", 4)
end_motifs_min_len = cnfg("end_motifs_min_len", None)
end_motifs_max_len = cnfg("end_motifs_max_len", None)
end_motifs_no_both_strands = cnfg("end_motifs_no_both_strands", True)
end_motifs_negative_strand = cnfg("end_motifs_negative_strand", False)
end_motifs_mapq = cnfg("end_motifs_mapq", 20)
end_motifs_workers = cnfg("end_motifs_workers", 1)

# interval-end-motifs
interval_end_motifs = cnfg("interval_end_motifs", False)
interval_end_motifs_refseq_file = cnfg("interval_end_motifs_refseq_file", None)
interval_end_motifs_interval_file = cnfg("interval_end_motifs_interval_file", None)
interval_end_motifs_kmer_length = cnfg("interval_end_motifs_kmer_length", 4)
interval_end_motifs_min_len = cnfg("interval_end_motifs_min_len", None)
interval_end_motifs_max_len = cnfg("interval_end_motifs_max_len", None)
interval_end_motifs_no_both_strands = cnfg("interval_end_motifs_no_both_strands", True)
interval_end_motifs_negative_strand = cnfg("interval_end_motifs_negative_strand", False)
interval_end_motifs_mapq = cnfg("interval_end_motifs_mapq", 20)
interval_end_motifs_workers = cnfg("interval_end_motifs_workers", None)

# wps
wps = cnfg("wps", False)
wps_chrom_sizes = cnfg("wps_chrom_sizes", None)
wps_site_bed = cnfg("wps_site_bed", False)
wps_interval_size = cnfg("wps_interval_size", 5000)
wps_window_size = cnfg("wps_window_size", 120)
wps_min_len = cnfg("wps_min_len", 120)
wps_max_len = cnfg("wps_max_len", 180)
wps_mapq = cnfg("wps_mapq", 30)
wps_workers = cnfg("wps_workers", 1)

# adjust-wps
adjust_wps = cnfg("adjust_wps", False)
adjust_wps_interval_file = cnfg("adjust_wps_interval_file", False)
adjust_wps_chrom_sizes = cnfg("adjust_wps_chrom_sizes", False)
adjust_wps_interval_size = cnfg("adjust_wps_interval_size", 5000)
adjust_wps_median_window_size = cnfg("adjust_wps_median_window_size", 1000)
adjust_wps_savgol_window_size = cnfg("adjust_wps_savgol_window_size", 21)
adjust_wps_savgol_poly_deg = cnfg("adjust_wps_savgol_poly_deg", 2)
adjust_wps_exclude_savgol = cnfg("adjust_wps_exclude_savgol", True)
adjust_wps_workers = cnfg("adjust_wps_workers", 1)
adjust_wps_mean = cnfg("adjust_wps_mean", False)
adjust_wps_subtract_edges = cnfg("adjust_wps_subtract_edges", False)
adjust_wps_edge_size = cnfg("adjust_wps_edge_size", None) # Set default back to 500 when fixed

# delfi
delfi = cnfg("delfi", False)
delfi_blacklist_file = cnfg("delfi_blacklist_file", None)
delfi_gap_file = cnfg("delfi_gap_file", None)
delfi_gap_reference_genome = cnfg("delfi_gap_reference_genome", "hg19")
delfi_chrom_sizes = cnfg("delfi_chrom_sizes", False)
delfi_reference_file = cnfg("delfi_reference_file", False)
delfi_bins_file = cnfg("delfi_bins_file", False)
delfi_no_gc_correct = cnfg("delfi_no_gc_correct", True)
delfi_keep_nocov = cnfg("delfi_keep_nocov", True)
delfi_no_merge_bins = cnfg("delfi_no_merge_bins", True)
delfi_window_size = cnfg("delfi_window_size", 5000000)
delfi_quality_threshold = cnfg("delfi_quality_threshold", 30)
delfi_workers = cnfg("delfi_workers", 1)

# cleavage-profile
cleavage_profile = cnfg("cleavage_profile", False)
cleavage_profile_interval_file = cnfg("cleavage_profile_interval_file", None)
cleavage_profile_chrom_sizes = cnfg("cleavage_profile_chrom_sizes", None)
cleavage_profile_min_len = cnfg("cleavage_profile_min_len", None)
cleavage_profile_max_len = cnfg("cleavage_profile_max_len", None)
cleavage_profile_mapq = cnfg("cleavage_profile_mapq", 20)
cleavage_profile_left = cnfg("cleavage_profile_left", 0)
cleavage_profile_right = cnfg("cleavage_profile_right", 0)
cleavage_profile_workers = cnfg("cleavage_profile_workers", 1)

if (adjust_wps and not wps):
    raise SystemExit("wps is required to run adjust-wps")

using_finaletoolkit = frag_length_bins or frag_length_intervals or coverage or end_motifs or interval_end_motifs or wps or adjust_wps or delfi

bed_files = glob.glob(cmb(in_dir, "*.gz"))
bam_files = glob.glob(cmb(in_dir, "*.bam"))
cram_files = glob.glob(cmb(in_dir, "*.cram"))

sample_bed = [os.path.splitext(os.path.basename(f))[0] for f in bed_files]
sample_bam = [os.path.splitext(os.path.basename(f))[0] for f in bam_files]
sample_cram = [os.path.splitext(os.path.basename(f))[0] for f in cram_files]
rule all:
    input:
        # BED output
        io(["final.gz", "final.gz.tbi"], sample_bed, not using_finaletoolkit),
        # BAM output
        io(["final.bam", "final.bam.bai"], sample_bam, not using_finaletoolkit),
        # CRAM output
        io(["final.cram", "final.cram.crai"], sample_cram, not using_finaletoolkit),

        # Finale Toolkit outputs (separate bins and intervals)

        # (BED)
        io(["gz.frag_length_bins.tsv", "bed.gz.frag_length_bins.png"], sample_bed, frag_length_bins),
        io(["gz.frag_length_intervals.bed"], sample_bed, frag_length_intervals),
        io(["gz.coverage.bed"], sample_bed, coverage),
        io(["gz.end_motifs.tsv"], sample_bed, end_motifs),
        io(["gz.interval_end_motifs.tsv"], sample_bed, interval_end_motifs),
        io(["gz.wps.bw"], sample_bed, wps),
        io(["gz.adjust_wps.bw"], sample_bed, adjust_wps),
        io(["gz.delfi.bed"], sample_bed, delfi),
        io(["gz.cleavage_profile.bw"], sample_bed, cleavage_profile),

        # (BAM)
        io(["bam.frag_length_bins.tsv","bam.frag_length_bins.png"], sample_bam, frag_length_bins),
        io(["bam.frag_length_intervals.bed"], sample_bam, frag_length_intervals),
        io(["bam.coverage.bed"], sample_bam, coverage),
        io(["bam.end_motifs.tsv"], sample_bam, end_motifs),
        io(["bam.interval_end_motifs.tsv"], sample_bam, interval_end_motifs),
        io(["bam.wps.bw"], sample_bam, wps),
        io(["bam.adjust_wps.bw"], sample_bam, adjust_wps),
        io(["bam.delfi.bed"], sample_bam, delfi),
        io(["bam.cleavage_profile.bw"], sample_bam, cleavage_profile),

        # (CRAM)
        io(["cram.frag_length_bins.tsv","cram.frag_length_bins.png"], sample_cram, frag_length_bins),
        io(["cram.frag_length_intervals.bed"], sample_cram, frag_length_intervals),
        io(["cram.coverage.bed"], sample_cram, coverage),
        io(["cram.end_motifs.tsv"], sample_cram, end_motifs),
        io(["cram.interval_end_motifs.tsv"], sample_cram, interval_end_motifs),
        io(["cram.wps.bw"], sample_cram, wps),
        io(["cram.adjust_wps.bw"], sample_cram, adjust_wps),
        io(["cram.delfi.bed"], sample_cram, delfi),
        io(["cram.cleavage_profile.bw"], sample_cram, cleavage_profile),

# STEP 0: Performance improvements by compressing and indexing the blacklist file
if blacklist:
    blacklist_input_path = cmb(sup_dir, blacklist)

    if not os.path.exists(blacklist_input_path):
        raise FileNotFoundError(f"Blacklist file not found: {blacklist_input_path}. The blacklist file must be in the supplement directory.")

    if blacklist.endswith(".gz"):
        blacklist_output_gz_path = blacklist_input_path
        blacklist_output_tbi_path = blacklist_output_gz_path + ".tbi"
        if not os.path.exists(blacklist_output_tbi_path):
            shell(f"tabix -p bed '{blacklist_output_gz_path}'")
    else:
        blacklist_output_gz_path = cmb(sup_dir, blacklist + ".gz")
        blacklist_output_tbi_path = blacklist_output_gz_path + ".tbi"
        if not os.path.exists(blacklist_output_gz_path):
            print(blacklist_input_path,blacklist_output_gz_path)
            shell(f"bgzip -@ 4 -c {blacklist_input_path} > {blacklist_output_gz_path}")
        if not os.path.exists(blacklist_output_tbi_path):
            shell(f"tabix -p bed '{blacklist_output_gz_path}'")
    blacklist = blacklist_output_gz_path

# STEP 1: Run filtering based on mapq scores
rule filter_bed_mapq:
    input:
        raw_bed=cmb(in_dir, "{sample}.gz"),
        raw_tbi=cmb(in_dir, "{sample}.gz.tbi")
    output:
        bed=cmb(out_dir, "{sample}.gz"),
        tbi=cmb(out_dir, "{sample}.gz.tbi")
    threads: 4 
    run:
        if mapq > 0:
            shell(f"zcat {input.raw_bed} | awk -F '\\t' '$4 >= {mapq}' | bgzip -c -@ {threads} > {output.bed}")
            shell(f"tabix -p bed {output.bed}")
        else:
            shell(f"cp {input.raw_bed} {output.bed}")
            shell(f"cp {input.raw_tbi} {output.tbi}")

rule filter_bam_mapq:
    input:
        raw_bam=cmb(in_dir, "{sample}.bam"),
        raw_bai=cmb(in_dir, "{sample}.bam.bai")
    output:
        bam=cmb(out_dir, "{sample}.bam"),
        bai=cmb(out_dir, "{sample}.bam.bai")
    threads: 4
    run:
        if mapq > 0:
            shell(f"samtools view -b -q {mapq} -@ {threads} {input.raw_bam} > {output.bam}")
            shell(f"samtools index {output.bam}")
        else:
            shell(f"cp {input.raw_bam} {output.bam}")
            shell(f"cp {input.raw_bai} {output.bai}")

rule filter_cram_mapq:
    input:
        raw_cram=cmb(in_dir, "{sample}.cram"),
        raw_crai=cmb(in_dir, "{sample}.cram.crai")
    output:
        cram=cmb(out_dir, "{sample}.cram"),
        crai=cmb(out_dir, "{sample}.cram.crai")
    threads: 4
    run:
        if mapq > 0:
            shell(f"samtools view -b -q {mapq} -@ {threads} {input.raw_cram} > {output.cram}")
            shell(f"samtools index {output.cram}")
        else:
            shell(f"cp {input.raw_cram} {output.cram}")
            shell(f"cp {input.raw_crai} {output.crai}")

# STEP 2: Remove regions based off of blacklist file
rule blacklist_bed:
    input:
        bed=rules.filter_bed_mapq.output.bed,
        tbi=rules.filter_bed_mapq.output.tbi,
    output:
        bed=cmb(out_dir, "{sample}.final.gz"),
        tbi=cmb(out_dir, "{sample}.final.gz.tbi"),
    params:
        blacklist=blacklist
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
        bam=rules.filter_bam_mapq.output.bam,
        bai=rules.filter_bam_mapq.output.bai,
    output:
        bam=cmb(out_dir, "{sample}.final.bam"),
        bai=cmb(out_dir, "{sample}.final.bam.bai"),
    params:
        blacklist=blacklist
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
        cram=rules.filter_cram_mapq.output.cram,
        crai=rules.filter_cram_mapq.output.cram,
    output:
        cram=cmb(out_dir, "{sample}.final.cram"),
        crai=cmb(out_dir, "{sample}.final.cram.crai"),
    params:
        blacklist=blacklist
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
        lambda wildcards: cmb(out_dir, f"{wildcards.sample}.final.{wildcards.ext}")
    output:
        tsv=cmb(out_dir, "{sample}.{ext}.frag_length_bins.tsv"),
        png=cmb(out_dir, "{sample}.{ext}.frag_length_bins.png")
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
        data=lambda wildcards:  cmb(out_dir, f"{wildcards.sample}.final.{wildcards.ext}"),
        intervals=cmb(sup_dir, f"{frag_length_interval_file}")
    output:
        cmb(out_dir, "{sample}.{ext}.frag_length_intervals.bed")
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
        data=lambda wildcards:  cmb(out_dir, f"{wildcards.sample}.final.{wildcards.ext}"),
        intervals=cmb(sup_dir, f"{coverage_interval_file}")
    output:
        cmb(out_dir, "{sample}.{ext}.coverage.bed")
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
        data=lambda wildcards: cmb(out_dir, f"{wildcards.sample}.final.{wildcards.ext}"),
        refseq=cmb(sup_dir, f"{end_motifs_refseq_file}")
    output:
        cmb(out_dir, "{sample}.{ext}.end_motifs.tsv")
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
        command_parts.append(f"-o {output} -v")

        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")

rule finaletoolkit_interval_end_motifs:
    input:
        data=lambda wildcards:  cmb(out_dir, f"{wildcards.sample}.final.{wildcards.ext}"),
        refseq=cmb(sup_dir, f"{interval_end_motifs_refseq_file}"),
        intervals=cmb(sup_dir, f"{interval_end_motifs_interval_file}")
    output:
        cmb(out_dir, "{sample}.{ext}.interval_end_motifs.tsv")
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
        command_parts.append(f"-o {output} -v")

        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")

rule wps:
    input:
        data=lambda wildcards: cmb(out_dir, f"{wildcards.sample}.final.{wildcards.ext}"),
        site_bed=cmb(sup_dir, f"{wps_site_bed}")
    output:
        cmb(out_dir, "{sample}.{ext}.wps.bw")
    run:
        command_parts = ["finaletoolkit", "wps", input.data, input.site_bed]
        if wps_chrom_sizes is not None:
            command_parts.append(f"-c {cmb(sup_dir,wps_chrom_sizes)}")
        if wps_interval_size is not None:
            command_parts.append(f"-i {wps_interval_size}")
        if wps_window_size is not None:
            command_parts.append(f"-W {wps_window_size}")
        if wps_min_len is not None:
            command_parts.append(f"-min {wps_min_len}")
        if wps_max_len is not None:
            command_parts.append(f"-max {wps_max_len}")
        if wps_mapq is not None:
            command_parts.append(f"-q {wps_mapq}")
        if wps_workers is not None:
            command_parts.append(f"-w {wps_workers}")
        command_parts.append(f"-o {output} -v")
        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")

rule adjust_wps:
    input:
        wps=lambda wildcards: cmb(out_dir, f"{wildcards.sample}.{wildcards.ext}.wps.bw"),
        intervals=cmb(sup_dir, f"{adjust_wps_interval_file}"),
        chrom_sizes=cmb(sup_dir, f"{adjust_wps_chrom_sizes}")
    output:
        cmb(out_dir, "{sample}.{ext}.adjust_wps.bw")
    run:
        command_parts = [
            "finaletoolkit",
            "adjust-wps",
            input.wps,
            input.intervals,
            input.chrom_sizes,
        ]

        if adjust_wps_interval_size is not None:
            command_parts.append(f"-i {adjust_wps_interval_size}")
        if adjust_wps_median_window_size is not None:
            command_parts.append(f"-m {adjust_wps_median_window_size}")
        if adjust_wps_savgol_window_size is not None:
            command_parts.append(f"-s {adjust_wps_savgol_window_size}")
        if adjust_wps_savgol_poly_deg is not None:
            command_parts.append(f"-p {adjust_wps_savgol_poly_deg}")
        if adjust_wps_exclude_savgol:
            command_parts.append("-S")
        if adjust_wps_workers is not None:
            command_parts.append(f"-w {adjust_wps_workers}")
        if adjust_wps_mean:
            command_parts.append("--mean")
        if adjust_wps_subtract_edges:
            command_parts.append("--subtract-edges")
            if adjust_wps_edge_size is not None:
                command_parts.append(f"--edge-size {adjust_wps_edge_size}")
        command_parts.append(f"-o {output} -v")
        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")

rule delfi:
    input:
        input_file=lambda wildcards: cmb(out_dir, f"{wildcards.sample}.final.{wildcards.ext}"),
        chrom_sizes=cmb(sup_dir, f"{delfi_chrom_sizes}"),
        reference_file=cmb(sup_dir, f"{delfi_reference_file}"),
        bins_file=cmb(sup_dir, f"{delfi_bins_file}"),
    output:
        cmb(out_dir, "{sample}.{ext}.delfi.bed")
    run:
        input_file = input.input_file
        if input_file.endswith(".gz"):
            # Convince finaletoolkit that it is indeed a .frag.gz file
            input_file = input_file.replace(".final.", ".")
            shell(f"mv {input.input_file} {input_file} && mv {input.input_file}.tbi {input_file}.tbi")

        command_parts = [
            "finaletoolkit",
            "delfi",
            input_file, # Use the modified input_file here
            input.chrom_sizes,
            input.reference_file,
            input.bins_file,
        ]

        if delfi_blacklist_file:
            command_parts.append(f"-b {cmb(sup_dir, delfi_blacklist_file)}")

        if not delfi_gap_file:
            gap_bed_output = cmb(sup_dir, f"{delfi_gap_reference_genome}.gap.bed")
            shell(f"finaletoolkit gap-bed {delfi_gap_reference_genome} {gap_bed_output}")
            command_parts.append(f"-g {gap_bed_output}")
        else:
            command_parts.append(f"-g {cmb(sup_dir, delfi_gap_file)}")

        if delfi_no_gc_correct:
            command_parts.append("-G")
        if delfi_keep_nocov:
            command_parts.append("-R")
        if delfi_no_merge_bins:
            command_parts.append("-M")
        if delfi_window_size is not None:
            command_parts.append(f"-s {delfi_window_size}")
        if delfi_quality_threshold is not None:
            command_parts.append(f"-q {delfi_quality_threshold}")
        if delfi_workers is not None:
            command_parts.append(f"-w {delfi_workers}")

        command_parts.append(f"-o {output} -v")

        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")
        if input_file.endswith(".gz"): # Rename it again so that snakemake doesn't lose track
            shell(f"mv {input_file} {input.input_file} && mv {input_file}.tbi {input.input_file}.tbi")

rule cleavage_profile:
    input:
        data=lambda wildcards: cmb(out_dir, f"{wildcards.sample}.final.{wildcards.ext}"),
        intervals=cmb(sup_dir, f"{cleavage_profile_interval_file}")
    output:
        cmb(out_dir, "{sample}.{ext}.cleavage_profile.bw")
    run:
        command_parts = [
            "finaletoolkit",
            "cleavage-profile",
            input.data,
            input.intervals,
        ]

        if cleavage_profile_chrom_sizes is not None:
            command_parts.append(f"-c {cmb(sup_dir, cleavage_profile_chrom_sizes)}")
        if cleavage_profile_min_len is not None:
            command_parts.append(f"-min {cleavage_profile_min_len}")
        if cleavage_profile_max_len is not None:
            command_parts.append(f"-max {cleavage_profile_max_len}")
        if cleavage_profile_mapq is not None:
            command_parts.append(f"-q {cleavage_profile_mapq}")
        if cleavage_profile_left is not None:
            command_parts.append(f"-l {cleavage_profile_left}")
        if cleavage_profile_right is not None:
            command_parts.append(f"-r {cleavage_profile_right}")
        if cleavage_profile_workers is not None:
            command_parts.append(f"-w {cleavage_profile_workers}")
        command_parts.append(f"-o {output} -v")
        command = " ".join(command_parts)
        print("Running ", command)
        shell(f"{command}")

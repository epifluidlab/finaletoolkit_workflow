import os
import glob
import subprocess

cmb = os.path.join
spl = os.path.splitext
cnfg = config.get


def check_tools():
    tools = {
    "finaletoolkit": "finaletoolkit --help",
    "bedtools": "bedtools --help",
    "htslib": "htsfile --help",
    "samtools": "samtools --help",
    "bedmappabilityfilter": "bedMappabilityFilter --help"
    }
    missing_tools = []
    for tool, command in tools.items():
        try:
            subprocess.run(command, shell=True, check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing_tools.append(tool)
    if missing_tools:
        raise SystemExit(f"Error: The following tools are not installed: {', '.join(missing_tools)}. Please install them:\n\n pip install finaletoolkit && conda install bioconda::bedtools bioconda::samtools bioconda::htslib kudosbeluga::bedmappabilityfilter\n\n")

check_tools()

out_dir = cnfg("output_dir", "output")
in_dir = cnfg("input_dir", "input")
sup_dir = cnfg("supplement_dir", "supplement") # Should contain interval, blacklist, .chrom.sizes, .2bit, and other files that supplement the input files
blk = cnfg("filter_file_blacklist_file", "")
wht = cnfg("filter_file_whitelist_file", "")
mappability_file = cnfg("mappability_file", None)
mappability_threshold = cnfg("mappability_threshold", 0)

filter_file = cnfg("filter_file", False)
filter_file_mapq = cnfg("filter_file_mapq", 0)
filter_file_min_length = cnfg("filter_file_min_length", 0)
filter_file_max_length = cnfg("filter_file_max_length", (1 << 31) - 1)
filter_file_blacklist = cmb(sup_dir, blk)
filter_file_whitelist = cmb(sup_dir, wht)
filter_file_intersect_policy = cnfg("filter_file_intersect_policy", "midpoint")
filter_file_workers = cnfg("filter_file_workers", 8)

# Finale Toolkit parameters
# frag-length-bins
frag_length_bins = cnfg("frag_length_bins", False)
frag_length_bins_mapq = cnfg("frag_length_bins_mapq", None)
frag_length_bins_policy = cnfg("frag_length_bins_intersect_policy", None)
frag_length_bins_bin_size = cnfg("frag_length_bins_bin_size", None)
frag_length_bins_min_len = cnfg("frag_length_bins_min_length", None)
frag_length_bins_max_len = cnfg("frag_length_bins_max_length", None)
frag_length_bins_chrom = cnfg("frag_length_bins_contig", "")
frag_length_bins_start = cnfg("frag_length_bins_start", None)
frag_length_bins_end = cnfg("frag_length_bins_stop", None)

# frag-length-intervals
frag_length_intervals = cnfg("frag_length_intervals", False)
frag_length_interval_file = cnfg("frag_length_interval_file", "")
frag_length_intervals_mapq = cnfg("frag_length_intervals_mapq", None)
frag_length_intervals_policy = cnfg("frag_length_intervals_intersect_policy", None)
frag_length_intervals_min_len = cnfg("frag_length_intervals_min_length", None)
frag_length_intervals_max_len = cnfg("frag_length_intervals_max_length", None)
frag_length_intervals_workers = cnfg("frag_length_intervals_workers", 1)

# coverage
coverage = cnfg("coverage", False)
coverage_interval_file = cnfg("coverage_interval_file", "")
coverage_mapq = cnfg("coverage_mapq", None)
coverage_intersect_policy = cnfg("coverage_intersect_policy", None)
coverage_min_len = cnfg("coverage_min_length", None)
coverage_max_len = cnfg("coverage_max_length", None)
coverage_workers = cnfg("coverage_workers", 1)
coverage_normalize = cnfg("coverage_normalize", False)
coverage_scale_factor = cnfg("coverage_scale_factor", None)

# end-motifs
end_motifs = cnfg("end_motifs", False)
end_motifs_refseq_file = cnfg("end_motifs_refseq_file", None)
end_motifs_kmer_length = cnfg("end_motifs_k", 4)
end_motifs_min_len = cnfg("end_motifs_min_length", None)
end_motifs_max_len = cnfg("end_motifs_max_length", None)
end_motifs_single_strand = cnfg("end_motifs_no_both_strands", False)
end_motifs_negative_strand = cnfg("end_motifs_negative_strand", False)
end_motifs_mapq = cnfg("end_motifs_mapq", 20)
end_motifs_workers = cnfg("end_motifs_workers", 1)

# interval-end-motifs
interval_end_motifs = cnfg("interval_end_motifs", False)
interval_end_motifs_refseq_file = cnfg("interval_end_motifs_refseq_file", None)
interval_end_motifs_interval_file = cnfg("interval_end_motifs_interval_file", "")
interval_end_motifs_kmer_length = cnfg("interval_end_motifs_kmer_length", 4)
interval_end_motifs_min_len = cnfg("interval_end_motifs_min_length", None)
interval_end_motifs_max_len = cnfg("interval_end_motifs_max_length", None)
interval_end_motifs_single_strand = cnfg("interval_end_motifs_single_strand", False)
interval_end_motifs_negative_strand = cnfg("interval_end_motifs_negative_strand", False)
interval_end_motifs_mapq = cnfg("interval_end_motifs_mapq", 20)
interval_end_motifs_workers = cnfg("interval_end_motifs_workers", None)

# mds
mds = cnfg("mds", False)
mds_sep = cnfg("mds_sep", " ")
mds_header = cnfg("mds_header", 0)

# interval-mds
interval_mds = cnfg("interval_mds", False)
interval_mds_sep = cnfg("interval_mds_sep", " ")
interval_mds_header = cnfg("interval_mds_header", 0)

# wps
wps = cnfg("wps", False)
wps_chrom_sizes = cnfg("wps_chrom_sizes", None)
wps_site_bed = cnfg("wps_site_bed", "")
wps_interval_size = cnfg("wps_interval_size", 5000)
wps_window_size = cnfg("wps_window_size", 120)
wps_min_len = cnfg("wps_min_length", 120)
wps_max_len = cnfg("wps_max_length", 180)
wps_mapq = cnfg("wps_mapq", 30)
wps_workers = cnfg("wps_workers", 1)

# adjust-wps
adjust_wps = cnfg("adjust_wps", False)
adjust_wps_interval_file = cnfg("adjust_wps_interval_file", "")
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
delfi_bins_file = cnfg("delfi_bins_file", "")
delfi_no_gc_correct = cnfg("delfi_no_gc_correct", True)
delfi_keep_nocov = cnfg("delfi_keep_nocov", True)
delfi_no_merge_bins = cnfg("delfi_no_merge_bins", True)
delfi_window_size = cnfg("delfi_window_size", 5000000)
delfi_mapq = cnfg("delfi_mapq", 30)
delfi_workers = cnfg("delfi_workers", 1)

# cleavage-profile
cleavage_profile = cnfg("cleavage_profile", False)
cleavage_profile_interval_file = cnfg("cleavage_profile_interval_file", "")
cleavage_profile_chrom_sizes = cnfg("cleavage_profile_chrom_sizes", None)
cleavage_profile_min_len = cnfg("cleavage_profile_min_length", None)
cleavage_profile_max_len = cnfg("cleavage_profile_max_length", None)
cleavage_profile_mapq = cnfg("cleavage_profile_mapq", 20)
cleavage_profile_left = cnfg("cleavage_profile_left", 0)
cleavage_profile_right = cnfg("cleavage_profile_right", 0)
cleavage_profile_workers = cnfg("cleavage_profile_workers", 1)

# agg-bw
agg_bw = cnfg("agg_bw", False)
agg_bw_interval_file = cnfg("agg_bw_interval_file", "")
agg_bw_median_window_size = cnfg("agg_bw_median_window_size", None)
agg_bw_mean = cnfg("agg_bw_mean", False)

if (adjust_wps and not wps):
    raise SystemExit("wps is required to run adjust-wps.")

if (filter_file_max_length < filter_file_min_length):
    raise SystemExit("Minimum read length must be smaller than the maximum read length.")

if (mds and not end_motifs):
    raise SystemExit("end-motifs is required to run mds.")

if (interval_mds and not interval_end_motifs):
    raise SystemExit("interval-end-motifs is required to run interval-mds.")

if (agg_bw and not cleavage_profile):
    raise SystemExit("cleavage-profile is required to run agg-bw.")

if blk and not os.path.exists(filter_file_blacklist):
    raise SystemExit(f"Blacklist file not found: {filter_file_blacklist}. The blacklist file must be in the supplement directory.")

if wht and not os.path.exists(filter_file_whitelist):
    raise SystemExit(f"Whitelist file not found: {filter_file_whitelist}. The whitelist file must be in the supplement directory.")

using_finaletoolkit = frag_length_bins or frag_length_intervals or coverage or end_motifs or interval_end_motifs or mds or interval_mds or wps or adjust_wps or delfi or cleavage_profile or agg_bw

bed_files = glob.glob(cmb(in_dir, "*.bed.gz"))
bam_files = glob.glob(cmb(in_dir, "*.bam"))
cram_files = glob.glob(cmb(in_dir, "*.cram"))

sample_bed = [spl(spl(os.path.basename(f))[0])[0] for f in bed_files]
sample_bam = [spl(os.path.basename(f))[0] for f in bam_files]
sample_cram = [spl(os.path.basename(f))[0] for f in cram_files]

sample_interval = [frag_length_interval_file, coverage_interval_file, cleavage_profile_interval_file, wps_site_bed, adjust_wps_interval_file, delfi_bins_file,interval_end_motifs_interval_file]
sample_interval_bins = [spl(fp)[0] for fp in sample_interval if fp and fp.endswith(".bins")]
sample_interval_bed = [spl(fp)[0] for fp in sample_interval if fp and fp.endswith(".bed")]

def io(endings, samples, condition,dir=out_dir):
    return expand(cmb(dir, "{sample}.{ending}"), sample=samples, ending=endings) if condition else []
rule all:
    input:
        # BED output
        io(["final.bed.gz","final.bed.gz.tbi"], sample_bed, True),
        # BAM output
        io(["final.bam","final.bam.bai"], sample_bam, True),
        # CRAM output
        io(["final.cram","final.cram.crai"], sample_cram, True),

        # Blacklist/Whitelist compression
        filter_file_blacklist+".tbi" if blk and blk.endswith(".gz") and not os.path.exists(filter_file_blacklist + ".tbi") else [],
        expand(filter_file_blacklist+"{endings}", endings=['.gz','.gz.tbi']) if blk and not blk.endswith(".gz") else [],
        filter_file_whitelist+".tbi" if wht and wht.endswith(".gz") and not os.path.exists(filter_file_whitelist + ".tbi") else [],
        expand(filter_file_whitelist+"{endings}", endings=['.gz','.gz.tbi']) if wht and not wht.endswith(".gz") else [],

        # Interval file output (mappability)
        io(["filtered.bed"], sample_interval_bed, len(sample_interval) > 0, sup_dir),
        io(["filtered.bins"], sample_interval_bins, len(sample_interval) > 0, sup_dir),

        # Finale Toolkit outputs
        # (BED)
        io(["final.bed.gz.frag_length_bins.tsv", "final.bed.gz.frag_length_bins.png"], sample_bed, frag_length_bins),
        io(["final.bed.gz.frag_length_intervals.bed"], sample_bed, frag_length_intervals),
        io(["final.bed.gz.coverage.bed"], sample_bed, coverage),
        io(["final.bed.gz.end_motifs.tsv"], sample_bed, end_motifs),
        io(["final.bed.gz.interval_end_motifs.tsv"], sample_bed, interval_end_motifs),
        io(["final.bed.gz.mds.txt"], sample_bed, mds),
        io(["final.bed.gz.interval_mds.tsv"], sample_bed, interval_mds),
        io(["final.bed.gz.wps.bw"], sample_bed, wps),
        io(["final.bed.gz.adjust_wps.bw"], sample_bed, adjust_wps),
        io(["final.bed.gz.delfi.bed"], sample_bed, delfi),
        io(["final.bed.gz.cleavage_profile.bw"], sample_bed, cleavage_profile),
        io(["final.bed.gz.agg_bw.wig"], sample_bed, agg_bw),
        # (BAM)
        io(["final.bam.frag_length_bins.tsv","final.bam.frag_length_bins.png"], sample_bam, frag_length_bins),
        io(["final.bam.frag_length_intervals.bed"], sample_bam, frag_length_intervals),
        io(["final.bam.coverage.bed"], sample_bam, coverage),
        io(["final.bam.end_motifs.tsv"], sample_bam, end_motifs),
        io(["final.bam.interval_end_motifs.tsv"], sample_bam, interval_end_motifs),
        io(["final.bam.mds.txt"], sample_bam, mds),
        io(["final.bam.interval_mds.tsv"], sample_bam, interval_mds),
        io(["final.bam.wps.bw"], sample_bam, wps),
        io(["final.bam.adjust_wps.bw"], sample_bam, adjust_wps),
        io(["final.bam.delfi.bed"], sample_bam, delfi),
        io(["final.bam.cleavage_profile.bw"], sample_bam, cleavage_profile),
        io(["final.bam.agg_bw.wig"], sample_bam, agg_bw),

        # (CRAM)
        io(["final.cram.frag_length_bins.tsv","cram.frag_length_bins.png"], sample_cram, frag_length_bins),
        io(["final.cram.frag_length_intervals.bed"], sample_cram, frag_length_intervals),
        io(["final.cram.coverage.bed"], sample_cram, coverage),
        io(["final.cram.end_motifs.tsv"], sample_cram, end_motifs),
        io(["final.cram.interval_end_motifs.tsv"], sample_cram, interval_end_motifs),
        io(["final.cram.mds.txt"], sample_cram, mds),
        io(["final.cram.interval_mds.tsv"], sample_cram, interval_mds),
        io(["final.cram.wps.bw"], sample_cram, wps),
        io(["final.cram.adjust_wps.bw"], sample_cram, adjust_wps),
        io(["final.cram.delfi.bed"], sample_cram, delfi),
        io(["final.cram.cleavage_profile.bw"], sample_cram, cleavage_profile),
        io(["final.cram.agg_bw.wig"], sample_cram, agg_bw),

# STEP 1: Compress and index the blacklist/whitelist file, just in case
rule blacklist:
    input: 
        filter_file_blacklist
    output:
        gz=filter_file_blacklist+".gz",
        tbi=filter_file_blacklist+".gz.tbi"
    threads: 4
    shell:
        """
            bgzip -@ {threads} -c {input} > {output.gz}
            tabix -p bed {output.gz}
        """
rule blacklist_index_only:
    input: 
        filter_file_blacklist
    output:
        tbi=filter_file_blacklist+".tbi"
    shell:
        """
            tabix -p bed {input}
        """
rule whitelist:
    input: 
        filter_file_whitelist
    output:
        gz=filter_file_whitelist+".gz",
        tbi=filter_file_whitelist+".gz.tbi"
    threads: 4
    shell:
        """
            bgzip -@ {threads} -c {input} > {output.gz}
            tabix -p bed {output.gz}
        """
rule whitelist_index_only:
    input: 
        filter_file_whitelist
    output:
        tbi=filter_file_whitelist+".tbi"
    shell:
        """
            tabix -p bed {input}
        """

# STEP 2: Run filtering using filter-file
def filter_file_helper(input, output, threads, ending):
        if filter_file:
            command = f"""finaletoolkit filter-file \
                {" -W " if wht else ""}{"" if not wht else filter_file_whitelist + ".gz" if not filter_file_whitelist.endswith(".gz") else filter_file_whitelist} \
                {" -B " if blk else ""}{"" if not blk else filter_file_blacklist + ".gz" if not filter_file_blacklist.endswith(".gz") else filter_file_blacklist} \
                -o {output.main} -q {filter_file_mapq} -min {filter_file_min_length} -max {filter_file_max_length} \
                -p {filter_file_intersect_policy} -w {threads} {input}"""

            print(f"Running command: {command}")   
            shell(command)
        else:
            shell(f"cp {input} {output.main}")
            shell(f"cp {input}.{ending} {output.index}")

rule filter_file_bed:
    input:
        cmb(in_dir, "{sample}.bed.gz")
    output:
        main=cmb(out_dir, "{sample}.final.bed.gz"),
        index=cmb(out_dir, "{sample}.final.bed.gz.tbi")
    threads: filter_file_workers
    run:
        filter_file_helper(input, output, threads, "tbi")

rule filter_file_bam:
    input:
        lambda wildcards: cmb(in_dir, f"{wildcards.sample}.bam")
    output:
        main=cmb(out_dir, "{sample}.final.bam"),
        index=cmb(out_dir, "{sample}.final.bam.bai")
    threads: filter_file_workers
    run:
        filter_file_helper(input, output, threads, "bai")

rule filter_file_cram:
    input:
        lambda wildcards: cmb(in_dir, f"{wildcards.sample}.cram")
    output:
        main=cmb(out_dir, "{sample}.final.cram"),
        index=cmb(out_dir, "{sample}.final.cram.crai")
    threads: filter_file_workers
    run:
        filter_file_helper(input, output, threads, "crai")
# STEP 3: Filter interval files using bedMappabilityFilter

rule filter_mappability_bed:
    input:
        lambda wildcards: cmb(sup_dir, f"{wildcards.sample}.{wildcards.ext}")
    output:
        cmb(sup_dir, "{sample}.filtered.{ext}"),
    threads: 4
    run:
        if mappability_file is not None:
            shell(f"bedMappabilityFilter --bigwig {cmb(sup_dir,mappability_file)} --bed {input} --output {output} --minimum-mappability {mappability_threshold} --threads {threads}")
        else:
            shell("mv {input} {output}")

# STEP 4: Regular Finaletoolkit commands. Using "is not None" to avoid falsy values accidentaly bypassing flags

rule frag_length_bins:
    input:
        cmb(out_dir, "{sample}")
    output:
        tsv=cmb(out_dir, "{sample}.frag_length_bins.tsv"),
        png=cmb(out_dir, "{sample}.frag_length_bins.png")
    threads: 1
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
        print("Running: ", command)
        shell(f"{command}")

rule frag_length_intervals:
    input:
        data= cmb(out_dir, "{sample}"),
        intervals=cmb(sup_dir, f"{spl(frag_length_interval_file)[0]}.filtered{spl(frag_length_interval_file)[1]}")
    output:
        cmb(out_dir, "{sample}.frag_length_intervals.bed")
    threads: frag_length_intervals_workers
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
        print("Running: ", command)
        shell(f"{command}")

rule coverage:
    input:
        data= cmb(out_dir, "{sample}"),
        intervals=cmb(sup_dir, f"{spl(coverage_interval_file)[0]}.filtered{spl(coverage_interval_file)[1]}")
    output:
        cmb(out_dir, "{sample}.coverage.bed")
    threads: coverage_workers
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
        print("Running: ", command)
        shell(f"{command}")

rule end_motifs:
    input:
        data= cmb(out_dir, "{sample}"),
        refseq=cmb(sup_dir, f"{end_motifs_refseq_file}")
    output:
        cmb(out_dir, "{sample}.end_motifs.tsv")
    threads: end_motifs_workers
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

        if end_motifs_single_strand:
            command_parts.append("-B")
            if end_motifs_negative_strand:
                command_parts.append("-n")

        if end_motifs_mapq is not None:
            command_parts.append(f"-q {end_motifs_mapq}")
        if end_motifs_workers is not None:
            command_parts.append(f"-w {end_motifs_workers}")
        command_parts.append(f"-o {output} -v")

        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule interval_end_motifs:
    input:
        data= cmb(out_dir, "{sample}"),
        refseq=cmb(sup_dir, f"{interval_end_motifs_refseq_file}"),
        intervals=cmb(sup_dir, f"{spl(interval_end_motifs_interval_file)[0]}.filtered{spl(interval_end_motifs_interval_file)[1]}")
    output:
        cmb(out_dir, "{sample}.interval_end_motifs.tsv")
    threads: interval_end_motifs_workers
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

        if interval_end_motifs_single_strand:
            command_parts.append("-B")
            if interval_end_motifs_negative_strand:
                command_parts.append("-n")

        if interval_end_motifs_mapq is not None:
            command_parts.append(f"-q {interval_end_motifs_mapq}")
        if interval_end_motifs_workers is not None:
            command_parts.append(f"-w {interval_end_motifs_workers}")
        command_parts.append(f"-o {output} -v")

        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule mds:
    input:
        cmb(out_dir, "{sample}.end_motifs.tsv")
    output:
        cmb(out_dir, "{sample}.mds.txt")
    run:
        command_parts = [
            "finaletoolkit",
            "mds",
            f"{input}"
        ]
        if mds_sep is not None and mds_sep != " ":
            command_parts.append(f"-s {mds_sep}")
        if mds_header is not None:
            command_parts.append(f"--header {mds_header}")

        command_parts.append(f"> {output}")
        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule interval_mds:
    input:
        cmb(out_dir, "{sample}.interval_end_motifs.tsv")
    output:
        cmb(out_dir, "{sample}.interval_mds.tsv")
    run:
        command_parts = [
            "finaletoolkit",
            "interval-mds",
            f"{input}",
            f"{output}"
        ]
        if interval_mds_sep is not None and interval_mds_sep != " ":
            command_parts.append(f"-s {interval_mds_sep}")
        if interval_mds_header is not None:
            command_parts.append(f"--header {interval_mds_header}")

        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule wps:
    input:
        data=cmb(out_dir, "{sample}"),
        site_bed=cmb(sup_dir, f"{spl(wps_site_bed)[0]}.filtered{spl(wps_site_bed)[1]}")
    output:
        cmb(out_dir, "{sample}.wps.bw")
    threads: wps_workers
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
        print("Running: ", command)
        shell(f"{command}")

rule adjust_wps:
    input:
        wps=cmb(out_dir, "{sample}.wps.bw"),
        intervals=cmb(sup_dir, f"{spl(adjust_wps_interval_file)[0]}.filtered{spl(adjust_wps_interval_file)[1]}"),
        chrom_sizes=cmb(sup_dir, f"{adjust_wps_chrom_sizes}")
    output:
        cmb(out_dir, "{sample}.adjust_wps.bw")
    threads: adjust_wps_workers
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
        print("Running: ", command)
        shell(f"{command}")

rule delfi:
    input:
        input_file=cmb(out_dir, "{sample}"),
        chrom_sizes=cmb(sup_dir, f"{delfi_chrom_sizes}"),
        reference_file=cmb(sup_dir, f"{delfi_reference_file}"),
        bins_file=cmb(sup_dir, f"{spl(delfi_bins_file)[0]}.filtered{spl(delfi_bins_file)[1]}"),
    output:
        cmb(out_dir, "{sample}.delfi.bed")
    threads: delfi_workers
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
        if delfi_mapq is not None:
            command_parts.append(f"-q {delfi_mapq}")
        if delfi_workers is not None:
            command_parts.append(f"-w {delfi_workers}")

        command_parts.append(f"-o {output} -v")

        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")
        if input_file.endswith(".gz"): # Rename it again so that snakemake doesn't lose track
            shell(f"mv {input_file} {input.input_file} && mv {input_file}.tbi {input.input_file}.tbi")

rule cleavage_profile:
    input:
        data= cmb(out_dir, "{sample}"),
        intervals=cmb(sup_dir, f"{spl(cleavage_profile_interval_file)[0]}.filtered{spl(cleavage_profile_interval_file)[1]}")
    output:
        cmb(out_dir, "{sample}.cleavage_profile.bw")
    threads: cleavage_profile_workers
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
        print("Running: ", command)
        shell(f"{command}")
rule agg_bw:
    input:
        data= cmb(out_dir, "{sample}.cleavage_profile.bw"),
        intervals=cmb(sup_dir, f"{spl(cleavage_profile_interval_file)[0]}.filtered{spl(cleavage_profile_interval_file)[1]}")
    output:
        cmb(out_dir, "{sample}.agg_bw.wig")
    threads: 1
    run:
        command_parts = [
            "finaletoolkit",
            "agg-bw",
            input.data,
            input.intervals,
        ]

        if output:
            command_parts.append(f"-o {output}")
        if agg_bw_median_window_size is not None:
            command_parts.append(f"-m {agg_bw_median_window_size}")
        if agg_bw_mean:
            command_parts.append("-a")
        command_parts.append("-v")

        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

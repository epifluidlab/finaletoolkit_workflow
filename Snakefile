import os
import sys
import glob
import subprocess

# Check whether or not the relevent dependencies exist in the environment.
def check_tools():
    tools = {
    "finaletoolkit": "finaletoolkit --help",
    "bedtools": "bedtools --help",
    "htslib": "htsfile --help",
    "samtools": "samtools --help",
    }
    missing_tools = []
    for tool, command in tools.items():
        try:
            subprocess.run(command, shell=True, check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing_tools.append(tool)
    if missing_tools:
        raise SystemExit(f"Error: The following tools are not installed: {', '.join(missing_tools)}.")

check_tools()

# Function to check whether an argument exists:
def exists(arg):
    return arg in config

# Set some basic arguments and locations and which commands will be run

out_dir = config.get("output_dir", "output")
in_dir = config.get("input_dir", "input")
sup_dir = config.get("supplement_dir", "supplement")
file_format = config.get("file_format","bed.gz")

filter_file = config.get("filter_file", False)
frag_length_bins = config.get("frag_length_bins", False)
frag_length_intervals = config.get("frag_length_intervals", False)
coverage = config.get("coverage", False)
end_motifs = config.get("end_motifs", False)
interval_end_motifs = config.get("interval_end_motifs", False)
mds = config.get("mds", False)
interval_mds = config.get("interval_mds", False)
wps = config.get("wps", False)
adjust_wps = config.get("adjust_wps", False)
delfi = config.get("delfi", False)
cleavage_profile = config.get("cleavage_profile", False)
agg_bw = config.get("agg_bw", False)

if (adjust_wps and not wps):
    raise SystemExit("wps is required to run adjust-wps.")

if (mds and not end_motifs):
    raise SystemExit("end-motifs is required to run mds.")

if (interval_mds and not interval_end_motifs):
    raise SystemExit("interval-end-motifs is required to run interval-mds.")

if (file_format != "bam" and file_format != "cram" and file_format != "bed.gz" and file_format != "frag.gz"):
    raise SystemExit('File format must be of type "bed.gz," "bam," or "cram."')

if file_format == "bam":
    index = "bai"
elif file_format == "cram":
    index = "crai"
elif file_format == "bed.gz" or file_format == "frag.gz":
    index = "tbi"

using_finaletoolkit = frag_length_bins or frag_length_intervals or coverage or end_motifs or interval_end_motifs or mds or interval_mds or wps or adjust_wps or delfi or cleavage_profile or agg_bw

# Set the expand function for the files.
def io(endings, samples, condition, dir=out_dir):
    return expand(os.path.join(dir, "{sample}.{ending}"), sample=samples, ending=endings) if condition else []

# Grab all of the sample files.
sample_files = { 
    'bed.gz': [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in glob.glob(os.path.join(in_dir, "*.bed.gz"))],
    'frag.gz': [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in glob.glob(os.path.join(in_dir, "*.frag.gz"))],
    'bam': [os.path.splitext(os.path.basename(f))[0] for f in glob.glob(os.path.join(in_dir, "*.bam"))],
    'cram': [os.path.splitext(os.path.basename(f))[0] for f in glob.glob(os.path.join(in_dir, "*.cram"))]
}

#Define all expected outputs and files.
rule all:
    input:
        # Output files for the filter-file.
        io(["filtered.bed.gz","filtered.bed.gz.tbi"], sample_files['bed.gz'], file_format == "bed.gz"),
        io(["filtered.frag.gz","filtered.frag.gz.tbi"], sample_files['frag.gz'], file_format == "frag.gz"),
        io(["filtered.bam","filtered.bam.bai"], sample_files['bam'], file_format == "bam"),
        io(["filtered.cram","filtered.cram.crai"], sample_files['cram'], file_format == "cram"),

        # Relevant output file if filtering out the mappability file.
        io(["filtered"], [config.get('interval_file', "")], exists('interval_file'), sup_dir),
        
        # Relevant FinaleToolkit function output files.
        *[
            io(["frag_length_bins.tsv", "frag_length_bins.png"], value, frag_length_bins) +
            io(["frag_length_intervals.bed"], value, frag_length_intervals) +
            io(["coverage.bed"], value, coverage) +
            io(["end_motifs.tsv"], value, end_motifs) +
            io(["interval_end_motifs.tsv"], value, interval_end_motifs) +
            io(["mds.txt"], value, mds) +
            io(["interval_mds.tsv"], value, interval_mds) +
            io(["wps.bw"], value, wps) +
            io(["adjust_wps.bw"], value, adjust_wps) +
            io(["delfi.bed"], value, delfi) +
            io(["cleavage_profile.bw"], value, cleavage_profile) +
            io(["agg_bw.wig"], value, agg_bw)
            for key, value in sample_files.items() if key == file_format
        ]

# STEP 1: Run filtering using filter-file

rule filter_file:
    input:
        main=os.path.join(in_dir, "{sample}.") + file_format,
        index=os.path.join(in_dir, "{sample}.") + file_format + "." + index
    output:
        main=os.path.join(out_dir, "{sample}.filtered.") + file_format,
        index=os.path.join(out_dir, "{sample}.filtered.") + file_format + "." + index,
    params:
        filter_file_mapq = config.get("filter_file_mapq", 0),
        filter_file_min_length = config.get("filter_file_min_length", 0),
        filter_file_max_length = config.get("filter_file_max_length", sys.maxsize),
        filter_file_blacklist = os.path.join(sup_dir, config['filter_file_blacklist']) if 'filter_file_blacklist' in config else None,
        filter_file_whitelist = os.path.join(sup_dir, config['filter_file_whitelist']) if 'filter_file_whitelist' in config else None,
        filter_file_intersect_policy = config.get("filter_file_intersect_policy", "midpoint"),
        filter_file_workers = config.get("workers", 1)
    run:
        if filter_file:
            command = f"""finaletoolkit filter-file \
                {f" -W {params.filter_file_whitelist}" if params.filter_file_whitelist else ""} \
                {f" -B {params.filter_file_blacklist}" if params.filter_file_blacklist else ""} \
                -o {output.main} -q {params.filter_file_mapq} -min {params.filter_file_min_length} -max {params.filter_file_max_length} \
                -p {params.filter_file_intersect_policy} -w {params.filter_file_workers} {input.main}"""
            shell(command)
        else:
            shell(f"cp {input.main} {output.main}")
            shell(f"cp {input.index} {output.index}")

rule filter_interval_file:
    input:
        interval = os.path.join(sup_dir, f"{config.get('interval_file', '')}")
    output:
        filtered = os.path.join(sup_dir, f"{config.get('interval_file', '')}" + ".filtered")
    params:
        mappability_file = os.path.join(sup_dir, f"{config.get('mappability_file')}"),
        threshold = config.get("mappability_threshold", 0)
    script:
        "scripts/mappability_filter.py"


# STEP 4: Regular Finaletoolkit commands.
rule frag_length_bins:
    input:
        os.path.join(out_dir, "{sample}.filtered.") + file_format
    output:
        tsv=os.path.join(out_dir, "{sample}.frag_length_bins.tsv"),
        png=os.path.join(out_dir, "{sample}.frag_length_bins.png")
    threads: 1
    params:
        mapq = config.get("frag_length_bins_mapq"),
        bin_size = config.get("frag_length_bins_bin_size"),
        policy = config.get("frag_length_bins_policy"),
        min_len = config.get("frag_length_bins_min_len"),
        max_len = config.get("frag_length_bins_max_len"),
        chrom = config.get("frag_length_bins_chrom"),
        start = config.get("frag_length_bins_start"),
        end = config.get("frag_length_bins_end")
    run:
        command_parts = [f"finaletoolkit frag-length-bins {input}"]
        if params.mapq is not None:
            command_parts.append(f"-q {params.mapq}")
        if params.bin_size is not None:
            command_parts.append(f"--bin-size {params.bin_size}")
        if params.policy is not None:
            command_parts.append(f"-p {params.policy}")
        if params.min_len is not None:
            command_parts.append(f"-min {params.min_len}")
        if params.max_len is not None:
            command_parts.append(f"-max {params.max_len}")
        if params.chrom is not None and params.chrom != "":
            command_parts.append(f"-c {params.chrom}")
        if params.start is not None:
            command_parts.append(f"-S {params.start}")
        if params.end is not None:
            command_parts.append(f"-E {params.end}")
        command_parts.append(f"-o {output.tsv}")
        command_parts.append(f"--histogram-path {output.png} -v")
        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule frag_length_intervals:
    input:
        data = os.path.join(out_dir, "{sample}.filtered.") + file_format,
        intervals = os.path.join(sup_dir, f"{config.get('interval_file')}" + ".filtered")
    output:
        os.path.join(out_dir, "{sample}.frag_length_intervals.bed")
    threads: config.get("frag_length_intervals_workers")
    params:
        min_len = config.get("frag_length_intervals_min_len"),
        max_len = config.get("frag_length_intervals_max_len"),
        policy = config.get("frag_length_intervals_policy"),
        mapq = config.get("frag_length_intervals_mapq"),
        workers = config.get("frag_length_intervals_workers")
    run:
        command_parts = [
            "finaletoolkit",
            "frag-length-intervals",
            input.data,
            input.intervals,
        ]

        if params.min_len is not None:
            command_parts.append(f"-min {params.min_len}")
        if params.max_len is not None:
            command_parts.append(f"-max {params.max_len}")
        if params.policy is not None:
            command_parts.append(f"-p {params.policy}")
        if params.mapq is not None:
            command_parts.append(f"-q {params.mapq}")
        if params.workers is not None:
            command_parts.append(f"-w {params.workers}")

        command_parts.append(f"-o {output[0]} -v")  # Ensure output indexing

        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")


rule coverage:
    input:
        data= os.path.join(out_dir, "{sample}.filtered.") + file_format,
        intervals = os.path.join(sup_dir, f"{config.get('interval_file')}" + ".filtered")
    output:
        os.path.join(out_dir, "{sample}.coverage.bed")
    threads: config.get("coverage_workers")
    params:
        min_len = config.get("coverage_min_len"),
        max_len = config.get("coverage_max_len"),
        normalize = config.get("coverage_normalize", False),
        scale_factor = config.get("coverage_scale_factor"),
        intersect_policy = config.get("coverage_intersect_policy"),
        mapq = config.get("coverage_mapq"),
        workers = config.get("coverage_workers")
    run:
        command_parts = [
            "finaletoolkit",
            "coverage",
            input.data,
            input.intervals,
        ]

        if params.min_len is not None:
            command_parts.append(f"-min {params.min_len}")
        if params.max_len is not None:
            command_parts.append(f"-max {params.max_len}")

        if params.normalize:
            command_parts.append("-n")
            if params.scale_factor is not None:
                command_parts.append(f"-s {params.scale_factor}")

        if params.intersect_policy is not None:
            command_parts.append(f"-p {params.intersect_policy}")

        if params.mapq is not None:
            command_parts.append(f"-q {params.mapq}")

        if params.workers is not None:
            command_parts.append(f"-w {params.workers}")

        command_parts.append(f"-o {output} -v")

        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule end_motifs:
    input:
        data = os.path.join(out_dir, "{sample}.filtered.") + file_format,
        refseq = os.path.join(sup_dir, f"{config.get('end_motifs_refseq_file')}")
    output:
        os.path.join(out_dir, "{sample}.end_motifs.tsv")
    threads: config.get("end_motifs_workers")
    params:
        kmer_length = config.get("end_motifs_kmer_length"),
        min_len = config.get("end_motifs_min_len"),
        max_len = config.get("end_motifs_max_len"),
        single_strand = config.get("end_motifs_single_strand", False),
        negative_strand = config.get("end_motifs_negative_strand", False),
        mapq = config.get("end_motifs_mapq"),
        workers = config.get("end_motifs_workers")
    run:
        command_parts = [
            "finaletoolkit",
            "end-motifs",
            input.data,
            input.refseq,
        ]

        if params.kmer_length is not None:
            command_parts.append(f"-k {params.kmer_length}")
        if params.min_len is not None:
            command_parts.append(f"-min {params.min_len}")
        if params.max_len is not None:
            command_parts.append(f"-max {params.max_len}")

        if params.single_strand:
            command_parts.append("-B")
            if params.negative_strand:
                command_parts.append("-n")

        if params.mapq is not None:
            command_parts.append(f"-q {params.mapq}")
        if params.workers is not None:
            command_parts.append(f"-w {params.workers}")
        command_parts.append(f"-o {output} -v")

        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule interval_end_motifs:
    input:
        data = os.path.join(out_dir, "{sample}.filtered.") + file_format,
        refseq = os.path.join(sup_dir, f"{config.get('interval_end_motifs_refseq_file')}"),
        intervals = os.path.join(sup_dir, f"{config.get('interval_file')}.filtered")
    output:
        os.path.join(out_dir, "{sample}.interval_end_motifs.tsv")
    threads: config.get("interval_end_motifs_workers")
    params:
        kmer_length = config.get("interval_end_motifs_kmer_length"),
        min_len = config.get("interval_end_motifs_min_len"),
        max_len = config.get("interval_end_motifs_max_len"),
        single_strand = config.get("interval_end_motifs_single_strand", False),
        negative_strand = config.get("interval_end_motifs_negative_strand", False),
        mapq = config.get("interval_end_motifs_mapq"),
        workers = config.get("interval_end_motifs_workers")
    run:
        command_parts = [
            "finaletoolkit",
            "interval-end-motifs",
            input.data,
            input.refseq,
            input.intervals,
        ]

        if params.kmer_length is not None:
            command_parts.append(f"-k {params.kmer_length}")
        if params.min_len is not None:
            command_parts.append(f"-min {params.min_len}")
        if params.max_len is not None:
            command_parts.append(f"-max {params.max_len}")

        if params.single_strand:
            command_parts.append("-B")
            if params.negative_strand:
                command_parts.append("-n")

        if params.mapq is not None:
            command_parts.append(f"-q {params.mapq}")
        if params.workers is not None:
            command_parts.append(f"-w {params.workers}")
        command_parts.append(f"-o {output[0]} -v")

        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule mds:
    input:
        os.path.join(out_dir, "{sample}.end_motifs.tsv")
    output:
        os.path.join(out_dir, "{sample}.mds.txt")
    params:
        sep = config.get("mds_sep", " "),
        header = config.get("mds_header")
    run:
        command_parts = [
            "finaletoolkit",
            "mds",
            input[0],  # Accessing the input file properly
        ]
        if params.sep is not None and params.sep != " ":
            command_parts.append(f"-s {params.sep}")
        if params.header is not None:
            command_parts.append(f"--header {params.header}")

        # Redirect output to the file
        command = " ".join(command_parts) + f" {output[0]}"
        print("Running: ", command)
        shell(f"{command}")

rule interval_mds:
    input:
        os.path.join(out_dir, "{sample}.interval_end_motifs.tsv")
    output:
        os.path.join(out_dir, "{sample}.interval_mds.tsv")
    params:
        sep = config.get("interval_mds_sep", " "),
        header = config.get("interval_mds_header")
    run:
        command_parts = [
            "finaletoolkit",
            "interval-mds",
            input[0],  # Accessing the input file properly
        ]
        if params.sep is not None and params.sep != " ":
            command_parts.append(f"-s {params.sep}")
        if params.header is not None:
            command_parts.append(f"--header {params.header}")

        # Redirect output to the file
        command = " ".join(command_parts) + f" {output[0]}"
        print("Running: ", command)
        shell(f"{command}")

rule wps:
    input:
        data = os.path.join(out_dir, "{sample}.filtered.") + file_format,
        site_bed = os.path.join(sup_dir, f"{config.get('wps_site_bed')}")
    output:
        os.path.join(out_dir, "{sample}.wps.bw")
    threads: config.get("wps_workers")
    params:
        chrom_sizes = config.get("wps_chrom_sizes"),
        interval_size = config.get("wps_interval_size"),
        window_size = config.get("wps_window_size"),
        min_len = config.get("wps_min_len"),
        max_len = config.get("wps_max_len"),
        mapq = config.get("wps_mapq"),
        workers = config.get("wps_workers")
    run:
        command_parts = ["finaletoolkit", "wps", input.data, input.site_bed]
        if params.chrom_sizes is not None:
            command_parts.append(f"-c {os.path.join(sup_dir, params.chrom_sizes)}")
        if params.interval_size is not None:
            command_parts.append(f"-i {params.interval_size}")
        if params.window_size is not None:
            command_parts.append(f"-W {params.window_size}")
        if params.min_len is not None:
            command_parts.append(f"-min {params.min_len}")
        if params.max_len is not None:
            command_parts.append(f"-max {params.max_len}")
        if params.mapq is not None:
            command_parts.append(f"-q {params.mapq}")
        if params.workers is not None:
            command_parts.append(f"-w {params.workers}")
        command_parts.append(f"-o {output[0]} -v")
        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule adjust_wps:
    input:
        wps = os.path.join(out_dir, "{sample}.wps.bw"),
        intervals = os.path.join(sup_dir, f"{config.get('interval_file')}" + ".filtered"),
        chrom_sizes = os.path.join(sup_dir, f"{config.get('adjust_wps_chrom_sizes')}")
    output:
        os.path.join(out_dir, "{sample}.adjust_wps.bw")
    threads: config.get("adjust_wps_workers")
    params:
        interval_size = config.get("adjust_wps_interval_size"),
        median_window_size = config.get("adjust_wps_median_window_size"),
        savgol_window_size = config.get("adjust_wps_savgol_window_size"),
        savgol_poly_deg = config.get("adjust_wps_savgol_poly_deg"),
        exclude_savgol = config.get("adjust_wps_exclude_savgol", False),
        workers = config.get("adjust_wps_workers"),
        mean = config.get("adjust_wps_mean", False),
        subtract_edges = config.get("adjust_wps_subtract_edges", False),
        edge_size = config.get("adjust_wps_edge_size")
    run:
        command_parts = [
            "finaletoolkit",
            "adjust-wps",
            input.wps,
            input.intervals,
            input.chrom_sizes,
        ]

        if params.interval_size is not None:
            command_parts.append(f"-i {params.interval_size}")
        if params.median_window_size is not None:
            command_parts.append(f"-m {params.median_window_size}")
        if params.savgol_window_size is not None:
            command_parts.append(f"-s {params.savgol_window_size}")
        if params.savgol_poly_deg is not None:
            command_parts.append(f"-p {params.savgol_poly_deg}")
        if params.exclude_savgol:
            command_parts.append("-S")
        if params.workers is not None:
            command_parts.append(f"-w {params.workers}")
        if params.mean:
            command_parts.append("--mean")
        if params.subtract_edges:
            command_parts.append("--subtract-edges")
            if params.edge_size is not None:
                command_parts.append(f"--edge-size {params.edge_size}")
        command_parts.append(f"-o {output[0]} -v")
        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule delfi:
    input:
        input_file = os.path.join(out_dir, "{sample}.filtered.") + file_format,
        chrom_sizes = os.path.join(sup_dir, f"{config.get('delfi_chrom_sizes')}"),
        reference_file = os.path.join(sup_dir, f"{config.get('delfi_reference_file')}"),
        bins_file = os.path.join(sup_dir, f"{config.get('delfi_bins_file')}")
    output:
        os.path.join(out_dir, "{sample}.delfi.bed")
    threads: config.get("delfi_workers")
    params:
        blacklist_file = config.get("delfi_blacklist_file"),
        gap_file = config.get("delfi_gap_file"),
        gap_reference_genome = config.get("delfi_gap_reference_genome"),
        no_gc_correct = config.get("delfi_no_gc_correct", False),
        keep_nocov = config.get("delfi_keep_nocov", False),
        no_merge_bins = config.get("delfi_no_merge_bins", False),
        window_size = config.get("delfi_window_size"),
        mapq = config.get("delfi_mapq"),
        workers = config.get("delfi_workers")
    run:

        command_parts = [
            "finaletoolkit",
            "delfi",
            input.input_file,  # Use the modified input_file here
            input.chrom_sizes,
            input.reference_file,
            input.bins_file
        ]

        if params.blacklist_file:
            command_parts.append(f"-b {os.path.join(sup_dir, params.blacklist_file)}")

        if not params.gap_file:
            gap_bed_output = os.path.join(sup_dir, f"{params.gap_reference_genome}.gap.bed")
            shell(f"finaletoolkit gap-bed {params.gap_reference_genome} {gap_bed_output}")
            command_parts.append(f"-g {gap_bed_output}")
        else:
            command_parts.append(f"-g {os.path.join(sup_dir, params.gap_file)}")

        if params.no_gc_correct:
            command_parts.append("-G")
        if params.keep_nocov:
            command_parts.append("-R")
        if params.no_merge_bins:
            command_parts.append("-M")
        if params.window_size is not None:
            command_parts.append(f"-s {params.window_size}")
        if params.mapq is not None:
            command_parts.append(f"-q {params.mapq}")
        if params.workers is not None:
            command_parts.append(f"-w {params.workers}")

        command_parts.append(f"-o {output[0]} -v")

        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule cleavage_profile:
    input:
        data = os.path.join(out_dir, "{sample}.filtered.") + file_format,
        intervals = os.path.join(sup_dir, f"{config.get('interval_file')}.filtered")
    output:
        os.path.join(out_dir, "{sample}.cleavage_profile.bw")
    threads: config.get("cleavage_profile_workers")
    params:
        chrom_sizes = config.get("cleavage_profile_chrom_sizes"),
        min_len = config.get("cleavage_profile_min_len"),
        max_len = config.get("cleavage_profile_max_len"),
        mapq = config.get("cleavage_profile_mapq"),
        left = config.get("cleavage_profile_left"),
        right = config.get("cleavage_profile_right"),
        workers = config.get("cleavage_profile_workers")
    run:
        command_parts = [
            "finaletoolkit",
            "cleavage-profile",
            input.data,
            input.intervals
        ]

        if params.chrom_sizes is not None:
            command_parts.append(f"-c {os.path.join(sup_dir, params.chrom_sizes)}")
        if params.min_len is not None:
            command_parts.append(f"-min {params.min_len}")
        if params.max_len is not None:
            command_parts.append(f"-max {params.max_len}")
        if params.mapq is not None:
            command_parts.append(f"-q {params.mapq}")
        if params.left is not None:
            command_parts.append(f"-l {params.left}")
        if params.right is not None:
            command_parts.append(f"-r {params.right}")
        if params.workers is not None:
            command_parts.append(f"-w {params.workers}")

        command_parts.append(f"-o {output[0]} -v")
        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

rule agg_bw:
    input:
        data = os.path.join(out_dir, "{sample}.cleavage_profile.bw"),
        intervals = os.path.join(sup_dir, f"{config.get('interval_file')}.filtered")
    output:
        os.path.join(out_dir, "{sample}.agg_bw.wig")
    threads: 1
    params:
        median_window_size = config.get("agg_bw_median_window_size"),
        mean = config.get("agg_bw_mean", False)
    run:
        command_parts = [
            "finaletoolkit",
            "agg-bw",
            input.data,
            input.intervals
        ]

        if params.median_window_size is not None:
            command_parts.append(f"-m {params.median_window_size}")
        if params.mean:
            command_parts.append("-a")

        command_parts.append(f"-o {output[0]}")
        command_parts.append("-v")

        command = " ".join(command_parts)
        print("Running: ", command)
        shell(f"{command}")

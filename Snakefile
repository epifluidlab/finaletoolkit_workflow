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

using_finaletoolkit = frag_length_bins or frag_length_intervals or coverage or end_motifs or interval_end_motifs or mds or interval_mds or wps or adjust_wps or delfi or cleavage_profile or agg_bw

# Set the expand function for the files.
def io(endings, samples, condition, dir=out_dir):
    return expand(os.path.join(dir, "{sample}.{ending}"), sample=samples, ending=endings) if condition else []

# Grab all of the sample files.
sample_files = { 
    'bed': [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in glob.glob(os.path.join(in_dir, "*.bed.gz"))],
    'bam': [os.path.splitext(os.path.basename(f))[0] for f in glob.glob(os.path.join(in_dir, "*.bam"))],
    'cram': [os.path.splitext(os.path.basename(f))[0] for f in glob.glob(os.path.join(in_dir, "*.cram"))]
}

#Define all expected outputs and files.
rule all:
    input:
        # Output files for the filter-file.
        io(["filtered.bed.gz","filtered.bed.gz.tbi"], sample_files['bed'], True),
        io(["filtered.bam","filtered.bam.bai"], sample_files['bam'], True),
        io(["filtered.cram","filtered.cram.crai"], sample_files['cram'], True),

        # Relevant output file if filtering out the mappability file.
        
        io([".filtered"], [config.get('interval_file', "")], exists('interval_file'), sup_dir),

        # Relevant FinaleToolkit function output files.
        # for key, value in sample_files.items():
        #     io([".frag_length_bins.tsv", "final.bed.gz.frag_length_bins.png"], value, frag_length_bins),
        #     io([".frag_length_intervals.bed"], value, frag_length_intervals),
        #     io([".coverage.bed"], value, coverage),
        #     io([".end_motifs.tsv"], value, end_motifs),
        #     io([".interval_end_motifs.tsv"], value, interval_end_motifs),
        #     io([".mds.txt"], value, mds),
        #     io([".interval_mds.tsv"], value, interval_mds),
        #     io([".wps.bw"], value, wps),
        #     io([".adjust_wps.bw"], value, adjust_wps),
        #     io([".delfi.bed"], value, delfi),
        #     io([".cleavage_profile.bw"], value, cleavage_profile),
        #     io([".agg_bw.wig"], value, agg_bw),

# STEP 1: Run filtering using filter-file
def filter_file_helper(input, output, params, ending):
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

ruleorder: filter_file_bam > filter_file_cram > filter_file_bed

rule filter_file_bam:
    input:
        main = lambda wildcards: os.path.join(in_dir, f"{wildcards.sample}.bam"),
        index = lambda wildcards: os.path.join(in_dir, f"{wildcards.sample}.bam.bai")
    output:
        main = os.path.join(out_dir, "{sample}.final.bam"),
        index = os.path.join(out_dir, "{sample}.final.bam.bai")
    params:
        filter_file_mapq = config.get("filter_file_mapq", 0),
        filter_file_min_length = config.get("filter_file_min_length", 0),
        filter_file_max_length = config.get("filter_file_max_length", sys.maxsize),
        filter_file_blacklist = lambda wildcards: os.path.join(sup_dir, config['filter_file_blacklist']) if 'filter_file_blacklist' in config else None,
        filter_file_whitelist = lambda wildcards: os.path.join(sup_dir, config['filter_file_whitelist']) if 'filter_file_whitelist' in config else None,
        filter_file_intersect_policy = config.get("filter_file_intersect_policy", "midpoint"),
        filter_file_workers = config.get("workers", 1)
    run:
        filter_file_helper(input, output, params, "bai")

rule filter_file_cram:
    input:
        main = lambda wildcards: os.path.join(in_dir, f"{wildcards.sample}.cram"),
        index = lambda wildcards: os.path.join(in_dir, f"{wildcards.sample}.cram.crai")
    output:
        main = os.path.join(out_dir, "{sample}.final.cram"),
        index = os.path.join(out_dir, "{sample}.final.cram.crai")
    params:
        filter_file_mapq = config.get("filter_file_mapq", 0),
        filter_file_min_length = config.get("filter_file_min_length", 0),
        filter_file_max_length = config.get("filter_file_max_length", sys.maxsize),
        filter_file_blacklist = lambda wildcards: os.path.join(sup_dir, config['filter_file_blacklist']) if 'filter_file_blacklist' in config else None,
        filter_file_whitelist = lambda wildcards: os.path.join(sup_dir, config['filter_file_whitelist']) if 'filter_file_whitelist' in config else None,
        filter_file_intersect_policy = config.get("filter_file_intersect_policy", "midpoint"),
        filter_file_workers = config.get("workers", 1)
    run:
        filter_file_helper(input, output, params, "crai")

rule filter_file_bed:
    input:
        main = lambda wildcards: os.path.join(in_dir, f"{wildcards.sample}.bed.gz"),
        index = lambda wildcards: os.path.join(in_dir, f"{wildcards.sample}.bed.gz.tbi")
    output:
        main = os.path.join(out_dir, "{sample}.filtered.bed.gz"),
        index = os.path.join(out_dir, "{sample}.filtered.bed.gz.tbi")
    params:
        filter_file_mapq = config.get("filter_file_mapq", 0),
        filter_file_min_length = config.get("filter_file_min_length", 0),
        filter_file_max_length = config.get("filter_file_max_length", sys.maxsize),
        filter_file_blacklist = lambda wildcards: os.path.join(sup_dir, config['filter_file_blacklist']) if 'filter_file_blacklist' in config else None,
        filter_file_whitelist = lambda wildcards: os.path.join(sup_dir, config['filter_file_whitelist']) if 'filter_file_whitelist' in config else None,
        filter_file_intersect_policy = config.get("filter_file_intersect_policy", "midpoint"),
        filter_file_workers = config.get("workers", 1)
    run:
        filter_file_helper(input, output, params, "tbi")

rule filter_interval_file:
    input:
        interval = config.get('interval_file', "")
    output:
        filtered = os.path.join(sup_dir, f"{config.get('interval_file')}.filtered")
    params:
        mappability_file = config.get("mappability_file"),
        threshold = config.get("mappability_threshold", 0)
    script:
        "scripts/mappability_filter.py"

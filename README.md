# Finaletoolkit Workflow

## Extract epigenomic features from several files in a folder 
**Features:**
* Workflow-specific: Mappability filtering of interval files, parallelization, slurm, BED/BAM/CRAM compatability

* Finaletoolkit: frag-length-bins, frag-length-intervals, coverage, end-motifs, interval-end-motifs, mds, interval-mds, wps, adjust-wps, gap-bed, delfi, cleavage-profile, agg-bw, filter-file


## Example usage
`snakemake --configfile params.yaml -c8 -j 2`

This runs the snakemake pipeline using parameters specified in params.yaml on 8 cores, with up to 2 rules running at a time.

Input files from the specified input directory are processed into the specified output directory. Blacklist and secondary input files for finaletoolkit go in the supplementary directory.

`snakemake --profile slurm_profile > snakemake.log 2>&1 &`

This runs the snakemake pipeline using flags from slurm_profile/config.yaml. In the sample config.yaml given in this repository, this would create a slurm submission that allows for up to 2GB per snakemake job, up to 4 jobs running in parallel, and up to 8 cores per job. It runs in the background (`&`) with output going to snakemake.log.

# Documentation

## Basic I/O
The `Snakefile` in this repository and a YAML file specifying the parameters to process your data must be present when running the workflow, in a manner similar to the example usage section above. 

YAML parameters (mappings) are of the form: `key: value`, where strings must be wrapped in double quotes, and each parameter is on its own line. This workflow requires the keys `input_dir` and `output_dir` to specify the folders in the working directory of which genomic data should be taken from and where it should be processed into. 

`supplement_dir` is also necessary to denote the folder of files that are not processed themselves in a command, but are needed to process input data. For example, you would put your blacklist, whitelist, mappability, and interval files in the folder specified by `supplement_dir`. 

## Filtering by mappability
`mappability_file` and `mappability_threshold` are used to filter interval files located in `supplement_dir`. `mappability_file` should locate the bigWig file in `supplement_dir` that filters all interval files, and `mappability_threshold` should be a floating point value from 0.0 to 1.0 to denote the minimum average mappability value that intervals will be filtered by as per the mappability file. Finaletoolkit commands that take in these interval files will recieve the interval file filtered by mappability.

## Finaletoolkit features

Finaletoolkit CLI commands and flags directly correspond to parameters in this workflow, with hyphens turning into underscores, and flags being seperated from their command by an underscore. For example, `adjust-wps` would be `adjust_wps`, and the max length flag of this command would be set through the key `adjust_wps_max_length`

If you want to use a Finaletoolkit command, set the converted command name in YAML to `True`. For example, if you wanted to run `adjust-wps` on your input files, include the line `adjust_wps: True`. Flags may be set through the naming scheme as specified in the paragraph above (`adjust_wps_max_length: 250`), and flags do not exist for `input_file`, `output_file`, deperecated flags, or `verbose`, which is on by default. 

## Filter-file

The `filter-file` command can be used in a different manner than the other Finaletoolkit commands. If only `filter-file` is set to run, then it will be the output of the workflow. However, if other finaletoolkit commands are set to run, then they will take in the output of `filter-file`.

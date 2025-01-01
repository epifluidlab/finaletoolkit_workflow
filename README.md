# Finaletoolkit Workflow

## Extract epigenomic features from several files in a folder 


* Done: Blacklisting regions (midpoints), mapq score filtering, mappability filtering, BED/BAM/CRAM compatability, parallelization, Finaletoolkit frag-length-bins, frag-length-intervals, coverage, end-motifs, interval-end-motifs, mds, interval-mds, wps, adjust-wps, gap-bed, delfi, cleavage-profile, agg-bw slurm

* Todo: None?

## Example usage
`snakemake --configfile params.yaml -c8 -j 2`

This runs the snakemake pipeline using parameters specified in params.yaml on 8 cores, with up to 2 rules running at a time.

Input files from the specified input directory are processed into the specified output directory. Blacklist and secondary input files for finaletoolkit go in the supplementary directory.

`snakemake --profile slurm_profile > snakemake.log 2>&1 &`

This runs the snakemake pipeline using flags from slurm_profile/config.yaml. In the sample config.yaml given in this repository, this would create a slurm submission that allows for up to 2GB per snakemake job, up to 4 jobs running in parallel, and up to 8 cores per job. It runs in the background (`&`) with output going to snakemake.log.

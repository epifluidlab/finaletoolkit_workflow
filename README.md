# Finaletoolkit Workflow

## Extract epigenomic features from several files in a folder 


* Done: Blacklisting regions (midpoints), mapq score filtering, BED/BAM/CRAM compatability, Finaletoolkit frag-length-bins, parallelization

* Todo: Combine BAM and CRAM rules, mappability, most other finaletoolkit commands, slurm/batch

## Example usage
`snakemake --configfile config.yaml -c8 -j 2`

This runs the snakemake pipeline using parameters specified in config.yaml on 8 cores, with up to 2 rules running at a time.

Input files from the specified input directory are processed into the specified output directory. Blacklist and secondary input files for finaletoolkit go in the supplementary directory.

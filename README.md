# Finaletoolkit Workflow

## Extract epigenomic features from several files in a folder 


* Done: Blacklisting regions (midpoints), mapq score filtering, BED/BAM/CRAM compatability, Finaletoolkit frag-length-bins

* Todo: Combine BAM and CRAM rules, mappability, most other finaletoolkit commands, parallelism, slurm/batch

## Example usage
`snakemake --configfile config.yaml -c1`
This runs the snakemake pipeline using parameters specified in config.yaml on one core. 

Input files from the specified input directory are processed into the specified output directory.

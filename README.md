# Finaletoolkit Workflow

## Extract epigenomic features from several files in a folder 


* Done: Blacklisting regions (midpoints), mapq score filtering, BED/BAM/CRAM compatability, parallelization, Finaletoolkit frag-length-bins, frag-length-intervals, coverage, end-motifs, interval-end-motifs, mds, interval-mds, wps, adjust-wps, gap-bed, delfi, cleavage-profile

* Todo: Mappability, slurm/batch, agg-bw (filter-bam functionality will be integrated into workflow)

## Example usage
`snakemake --configfile config.yaml -c8 -j 2`

This runs the snakemake pipeline using parameters specified in config.yaml on 8 cores, with up to 2 rules running at a time.

Input files from the specified input directory are processed into the specified output directory. Blacklist and secondary input files for finaletoolkit go in the supplementary directory.

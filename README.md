
# Finaletoolkit Workflow

This Snakemake workflow automates the extraction of epigenomic features using Finaletoolkit, supporting parallel processing, SLURM, and common genomic file formats.

## Key Features

*   **Finaletoolkit Support:** Implements most Finaletoolkit CLI commands (excluding `gap-bed` and `delfi-gc-correct`).
*   **File Compatibility:** Works with BED, BAM, and CRAM files.
*   **Mappability Filtering:** Filters interval files based on mappability scores.
*   **Parallelization:** Supports multi-core processing.
*   **SLURM Integration:** Enables job submission to SLURM clusters.

## Dependencies

This workflow relies on the following tools being installed and accessible by your system PATH. We recommend that you install FinaleToolkit through `pip` and the other packages through `conda` in the Bioconda channel:

* `finaletoolkit`: A command-line tool for epigenomic analysis.
* `bedtools`: A suite of utilities for working with genomic intervals.
* `htslib`: A library that includes `bgzip`, necessary to GZIP uncompressed BED files. 
* `samtools`: A set of tools for manipulating and analyzing sequencing BAM/CRAM data

## Quick Start

1.  **Configuration:**  Create a `params.yaml` file defining your input, output, and processing options (reference below sections).
2.  **Basic Execution:** Run the workflow with `snakemake --configfile params.yaml -c <cores> -j <jobs>`.
   * `-c`: Number of CPU cores to use.
   * `-j`: Maximum number of concurrent jobs.
3.  **SLURM Execution:** Submit to SLURM to run in the background with `snakemake --profile slurm_profile > snakemake.log 2>&1 &` (see `slurm_profile/config.yaml` for default settings).

## Workflow Structure

*   **Input:**
    *   Genomic data files are located in the directory specified by `input_dir`.
    *   Supplemental files (blacklist, mappability, intervals) are located in the directory specified by `supplement_dir`.
*   **Output:** Processed files are written to the directory specified by `output_dir`.
*   **Configuration:** `params.yaml` dictates workflow parameters.

## YAML Parameters

*   **Required:**
    *   `input_dir`: Path to the input directory. Defaults to `input` if not specified
    *   `output_dir`: Path to the output directory. Defaults to `output` if not specified.
    *    `file_format`: `"bed.gz"`, `"bam"`, or `"cram"` indicating the format of the input files. Defaults to `bed.gz` if not specified.

*   **Optional:**
    *   `supplement_dir`: Path to supplemental files directory. Defaults to `supplement` if not specified. 
    *   `mappability_file`: Name of the bigWig mappability file in `supplement_dir`.
    *    `mappability_threshold`: Minimum average mappability score (0.0-1.0) for interval filtering.
    *  `interval_file`: Path to interval file in `supplement_dir`.
    *   `finaletoolkit_command: True/False`: Enables a specific Finaletoolkit command, using hyphens replaced by underscores (e.g., `adjust-wps` becomes `adjust_wps: True`).
    *   `finaletoolkit_command_flag: value`: Sets flags for a Finaletoolkit command (e.g., `adjust_wps_max_length: 250`). Flags that take input files, output files, or `verbose` flags do not exist here.  `mapping_quality` is shortened to `mapq` for flags (e.g., `coverage_mapping_quality` becomes `coverage_mapq`).

## Output File Naming

*   **Filtered Files:** Files are always are given a `.filtered` extension before the file format when passed into the output directory (e.g., `file.filtered.bed.gz`).
*   **Command-Processed Files:** Files processed by a Finaletoolkit command have the command name (with underscores) inserted before their format (e.g., `file.frag_length_bins.bed.gz`).
*   **Multiple Commands:** Input files will be processed for each enabled Finaletoolkit command.

## Mappability Filtering

*   Uses ``mappability_file`` and ``mappability_threshold`` (float from 0-1) to filter intervals specified by ``interval_file``.
*   Interval files used in Finaletoolkit commands are pre-filtered by mappability quallity to at least the threshold if set. 

## Using Finaletoolkit Commands

*   Finaletoolkit commands are specified in `params.yaml` with underscores instead of hyphens.
*   Set command flags by appending `_<flag_name>` to the converted command name.
*   **Dependencies:**  The workflow respects command dependencies.  For example, `adjust-wps` requires `wps` output, and `mds` needs `end-motifs`.

## `filter-file` Command

*   If only `filter-file` is set to `True`, the output of the workflow will be the filtered files.
*   If any other Finaletoolkit commands are set, they will use the output of `filter-file` as their input.

## Notes

*   This workflow uses `verbose` for Finaletoolkit commands by default.
*  Deprecated flags are not used in this workflow.

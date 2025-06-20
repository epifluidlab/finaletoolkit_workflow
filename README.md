
# FinaleToolkit Workflow

This Snakemake workflow automates the extraction of epigenomic features using Finaletoolkit, supporting parallel processing, SLURM, and common genomic file formats.

## Key Features

*   **Finaletoolkit Support:** Implements most Finaletoolkit CLI commands (excluding `gap-bed` and `delfi-gc-correct`).
*   **File Compatibility:** Works with BED, BAM, and CRAM files.
*   **Mappability Filtering:** Filters interval files based on mappability scores.
*   **Parallelization:** Supports multi-core processing.
*   **SLURM Integration:** Enables job submission to SLURM clusters.

## Installation

Reference the following commands for setup and execution:

```bash
$ git clone https://github.com/epifluidlab/finaletoolkit_workflow # Download the repository containing the workflow
$ cd finaletoolkit_workflow # Enter the repository folder
$ conda env create -f environment.yml # Create environment with relevant conda packages
$ conda activate finaletoolkit_workflow # Use environment for finaletoolkit-workflow
$ pip install finaletoolkit # Install finaletoolkit seperately through pip 
$ snakemake --configfile params.yaml --cores 4 --jobs 2 # Run with parameters set in params.yaml
```

## Dependencies

This workflow relies on the following tools being installed and accessible by your system PATH. FinaleToolkit must be installed through `pip` and the other packages through bioconda in `conda` (already installed if you activated the conda environment from `environment.yml`)

* `finaletoolkit`: A command-line tool for epigenomic feature extraction.
* `snakemake`: A workflow engine that determines which operations ("rules") to carry out on genomic files.
* `bedtools`: A suite of utilities for working with and manipulating genomic intervals.
* `htslib`: A library that includes `bgzip`, necessary to GZIP uncompressed BED files. 
* `samtools`: A set of tools for manipulating and analyzing sequencing BAM/CRAM data

## Quick Start

1.  **Configuration:**  Create a `params.yaml` file defining your input, output, and processing options (reference below sections).
2.  **Basic Execution:** Run the workflow in the directory with the `Snakefile` present through the following command:
```bash
cd finaletoolkit_workflow # Enter the folder with the workflow Snakefile

# --cores: Number of CPU cores to use.
# --jobs: Maximum number of concurrent jobs (limited by --cores).
snakemake --configfile params.yaml --cores <cores> --jobs <jobs>
```
3.  **SLURM Execution:** Submit to SLURM to run the workflow through the command below (see `slurm_profile/config.yaml` for default settings).
```bash
cd finaletoolkit_workflow # Enter the folder with the workflow Snakefile

# Runs the command through params specified in slurm_profile/config.yaml in the background (&),
# Redirects all command-related output to snakemake.log
snakemake  --profile slurm_profile > snakemake.log 2>&1 &
```

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
    *    `file_format`: `"bed.gz"`, `"frag.gz"`, `"bam"`, or `"cram"` indicating the format of the input files. Defaults to `bed.gz` if not specified.

*   **Optional:**
    *   `supplement_dir`: Path to supplemental files directory. Defaults to `supplement` if not specified. 
    *   `mappability_file`: Name of the bigWig mappability file in `supplement_dir`.
    *    `mappability_threshold`: Minimum average mappability score (0.0-1.0) for interval filtering.
    *  `interval_file`: Path to interval file in `supplement_dir`.
    *   `finaletoolkit_command: True/False`: Enables a specific Finaletoolkit command, using hyphens replaced by underscores (e.g., `adjust-wps` becomes `adjust_wps: True`).
    *   `finaletoolkit_command_flag: value`: Sets flags for a Finaletoolkit command (e.g., `adjust_wps_max_length: 250`). Flags that take input files, output files, or `verbose` flags do not exist here.  `mapping_quality` is shortened to `mapq` for flags (e.g., `coverage_mapping_quality` becomes `coverage_mapq`).

## Output File Naming

*   **Filtered Files:** Files are always given a `.filtered` extension before the file format when passed into the output directory (e.g., `file.filtered.bed.gz`).
*   **Command-Processed Files:** Files processed by a Finaletoolkit command have the command name (with underscores) inserted before their format (e.g., `file.frag_length_bins.bed.gz`).
*   **Multiple Commands:** Input files will be processed for each enabled Finaletoolkit command.

## Mappability Filtering

*   Uses ``mappability_file`` and ``mappability_threshold`` (float from 0-1) to filter intervals specified by ``interval_file``.
*   Interval files used in Finaletoolkit commands are pre-filtered by mappability quality to at least the threshold if set. 

## Using Finaletoolkit Commands

*   Finaletoolkit commands are specified in `params.yaml` with underscores instead of hyphens.
*   Set command flags by appending `_<flag_name>` to the converted command name.
*   This workflow respects command dependencies.  For example, `adjust-wps` requires `wps` output, and `mds` needs `end-motifs`.

## `filter-file` Command

*   If only `filter-file` is set to `True`, the output of the workflow will be only the filtered files.
*   If any other Finaletoolkit commands are set, they will use the output of `filter-file` as their input.

## Notes

*   This workflow uses `verbose` for Finaletoolkit commands by default.
*   Deprecated flags cannot be used in this workflow.


## Citation
Li J*, Bandaru R*, Liu Y (2024) FinaleToolkit: Accelerating Cell-Free DNA Fragmentation Analysis with a High-Speed Computational Toolkit. BioRxiv Preprint [![Static Badge](https://img.shields.io/badge/DOI-10.1101%2F2024.05.29.596414-blue?style=flat-square)](https://doi.org/10.1101/2024.05.29.596414)

## Contact
- Kundan Baliga: kundanbal2969@k12.ipsd.org
- Ravi Bandaru: ravi.bandaru@northwestern.edu
- Yaping Liu: yaping@northwestern.edu

## License
This project falls under an MIT license. See the included `LICENSE` file for details.



# FinaleToolkit Workflow

This Snakemake workflow automates the extraction of epigenomic features using Finaletoolkit, supporting parallel processing, SLURM, and common genomic file formats. It targets **FinaleToolkit 1.0.0**, which introduced a redesigned CLI (see the [Migration to FinaleToolkit 1.0.0](#migration-to-finaletoolkit-100) section).

## Key Features

*   **Finaletoolkit Support:** Implements most Finaletoolkit 1.0.0 CLI commands. (`delfi-gc-correct` was removed in 1.0.0 — GC correction now runs automatically inside `delfi`; `gap-bed` is not exposed directly but is used internally by the `delfi` step.)
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
$ snakemake --configfile params.yaml --cores 4 --jobs 2 # Run with parameters set in params.yaml
```

## Dependencies

This workflow relies on the following tools being installed and accessible by your system PATH. FinaleToolkit and other
packages may be installed through bioconda in `conda` (already installed if you activated the conda environment from
`environment.yml`)

* `finaletoolkit`: A command-line tool for epigenomic feature extraction.
* `snakemake`: A workflow engine that determines which operations ("rules") to carry out on genomic files.
* `bedtools`: A suite of utilities for working with and manipulating genomic intervals.
* `htslib`: A library that includes `bgzip`, necessary to GZIP uncompressed BED files. 
* `samtools`: A set of tools for manipulating and analyzing sequencing BAM/CRAM data

## Quick Start

1.  **Configuration:**  Create a `params.yaml` file defining your input, output, and processing options (reference below
sections).
2.  **Basic Execution:** Run the workflow in the directory with the `Snakefile` present through the following command:
```bash
cd finaletoolkit_workflow # Enter the folder with the workflow Snakefile

# --cores: Number of CPU cores to use.
# --jobs: Maximum number of concurrent jobs (limited by --cores).
snakemake --configfile params.yaml --cores <cores> --jobs <jobs>
```
3.  **SLURM Execution:** Before using this workflow with SLURM, first install
the snakemake executor plugin: slurm:

```bash
conda activate finaletoolkit_workflow # activate the conda environment if not already activated
conda install bioconda::snakemake-executor-plugin-slurm # install executor plugin from bioconda
```

Submit to SLURM to run the workflow through the command below (see `slurm_profile/config.yaml`
for default settings).
```bash
cd finaletoolkit_workflow # Enter the folder with the workflow Snakefile

# Runs the command through params specified in slurm_profile/config.yaml in the background (&),
# Redirects all command-related output to snakemake.log
snakemake  --profile slurm_profile > snakemake.log 2>&1 &
```

## Workflow Structure

*   **Input:**
    *   Genomic data files are located in the directory specified by `input_dir`.
    *   Supplemental files (blacklist, mappability, intervals) are located in the directory specified by
    `supplement_dir`.
*   **Output:** Processed files are written to the directory specified by `output_dir`.
*   **Configuration:** `params.yaml` dictates workflow parameters.

## YAML Parameters

*   **Required:**
    *   `input_dir`: Path to the input directory. Defaults to `input` if not specified
    *   `output_dir`: Path to the output directory. Defaults to `output` if not specified.
    *    `file_format`: `"bed.gz"`, `"frag.gz"`, `"bam"`, or `"cram"` indicating the format of the input files. Defaults
    to `bed.gz` if not specified.

*   **Optional:**
    *   `supplement_dir`: Path to supplemental files directory. Defaults to `supplement` if not specified. 
    *   `mappability_file`: Name of the bigWig mappability file in `supplement_dir`.
    *    `mappability_threshold`: Minimum average mappability score (0.0-1.0) for interval filtering.
    *  `interval_file`: Path to interval file in `supplement_dir`.
    *   `finaletoolkit_command: True/False`: Enables a specific Finaletoolkit command, using hyphens replaced by
    underscores (e.g., `adjust-wps` becomes `adjust_wps: True`).
    *   `finaletoolkit_command_flag: value`: Sets flags for a Finaletoolkit command (e.g.,
    `adjust_wps_max_length: 250`). Flags that take input files, output files, or `verbose` flags do not exist here.  `mapping_quality` is shortened to `mapq` for flags (e.g., `coverage_mapping_quality` becomes `coverage_mapq`). Config keys follow FinaleToolkit 1.0.0 terminology — worker counts use `_threads`, fragment bounds use `_min_length`/`_max_length`, motif strand selection uses a single `_strand` key (`both`/`forward`/`reverse`), and cleavage padding uses `_pad_left`/`_pad_right`.

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
*   This workflow respects command dependencies.  For example, `adjust-wps` requires `wps` output, `mds` needs `end-motifs`, and `regional_mds` (renamed from `interval-mds` in 1.0.0) needs `interval-end-motifs`.

## `filter-file` Command

*   If only `filter-file` is set to `True`, the output of the workflow will be only the filtered files.
*   If any other Finaletoolkit commands are set, they will use the output of `filter-file` as their input.

## Migration to FinaleToolkit 1.0.0

FinaleToolkit 1.0.0 redesigned its CLI, so several `params.yaml` config keys were renamed to match. If you are upgrading an existing `params.yaml`, update the following keys:

| Old key | New key |
| --- | --- |
| `<command>_workers`, `workers` | `<command>_threads`, `threads` |
| `<command>_min_len` / `<command>_max_len` | `<command>_min_length` / `<command>_max_length` |
| `<motif>_single_strand` + `<motif>_negative_strand` | `<motif>_strand` (`both` / `forward` / `reverse`) |
| `cleavage_profile_left` / `cleavage_profile_right` | `cleavage_profile_pad_left` / `cleavage_profile_pad_right` |
| `frag_length_bins_short_fraction` | `frag_length_bins_short_threshold` |
| `frag_length_intervals_short_reads` | `frag_length_intervals_short_threshold` |
| `adjust_wps_exclude_savgol: True` | `adjust_wps_savgol: False` (meaning inverted; savgol is on by default) |
| `delfi_window_size` | `delfi_merge_size` |
| `delfi_keep_nocov: True` | `delfi_remove_nocov: False` (meaning inverted; removal is on by default) |
| `delfi_no_merge_bins: True` | `delfi_merge_bins: False` (meaning inverted; merging is on by default) |
| `interval_mds` (and `interval_mds_sep` / `interval_mds_header`) | `regional_mds` (and `regional_mds_sep` / `regional_mds_header`) |

Motif strand selection: use `_strand: "forward"` for the old `single_strand: True` (positive strand only) and `_strand: "reverse"` for `single_strand: True` + `negative_strand: True`. The default `both` matches the old default.

Output note: the `regional_mds` outputs land in `output/regional_mds/*.regional_mds.tsv` (previously `interval_mds`).

## Notes

*   This workflow uses `verbose` for Finaletoolkit commands by default.
*   Deprecated flags cannot be used in this workflow.


## Citation
Li JW, Bandaru R, Baliga K, Liu Y. FinaleToolkit: accelerating cell-free DNA fragmentation analysis with a high-speed computational toolkit. Bioinform Adv. 2025;5(1):vbaf236. Published 2025 Sep 26. doi:10.1093/bioadv/vbaf236

## Contact
- Kundan Baliga: kundanbal2969@k12.ipsd.org
- Ravi Bandaru: ravi.bandaru@northwestern.edu
- Yaping Liu: yaping@northwestern.edu

## License
This project falls under an MIT license. See the included `LICENSE` file for details.


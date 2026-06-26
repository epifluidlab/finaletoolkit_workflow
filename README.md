<!-- Swap the placeholder src below for a kawaii irasutoya image (DNA / test tube),
     matching the title style of the other ravibandaru-lab repos. -->
# <img alt="dna with letters FT" src="https://github.com/epifluidlab/FinaleToolkit/blob/b99b38e22a3b07ee9b7e0fa44a488a1eb5442efe/docs/_static/finaletoolkit_logo_rounded.png?raw=true" height="60"> finaletoolkit_workflow

A Snakemake workflow that automates cell-free DNA fragmentation feature extraction with
[FinaleToolkit](https://github.com/epifluidlab/FinaleToolkit) — supporting **hg38** and
**T2T-CHM13**, parallel processing, SLURM, and BED/BAM/CRAM inputs.

[![Docs](https://img.shields.io/badge/docs-epifluidlab.github.io-836eaa?style=flat-square)](https://epifluidlab.github.io/FinaleToolkit/documentation/workflow/)

### Contents

|     Installation     |     Reference Setup     |     Usage     |     Parameters     |
|:--------------------:|:-----------------------:|:-------------:|:------------------:|
| [Install](#installation) | [Reference setup](#reference--supplement-setup) | [Usage](#usage) | [YAML parameters](#yaml-parameters) |
| [Genomes](#genome-support) | [Mappability](#mappability) | [SLURM](#slurm-execution) | [Output naming](#output-file-naming) |

### Genome support

Two genomes are supported out of the box, each with a runnable example config and a one-command
reference setup:

| Genome | Example config | Setup command |
|:--|:--|:--|
| hg38 | [`params.hg38.yaml`](./params.hg38.yaml) | `scripts/setup_reference.sh hg38 supplement 500` |
| T2T-CHM13 (hs1) | [`params.t2t-chm13.yaml`](./params.t2t-chm13.yaml) | `scripts/setup_reference.sh t2t-chm13 supplement 500` |

### Installation

```bash
git clone https://github.com/epifluidlab/finaletoolkit_workflow
cd finaletoolkit_workflow
conda env create -f environment.yml
conda activate finaletoolkit_workflow
```

Core tools (installed by the environment): `finaletoolkit`, `snakemake`, `bedtools`, `htslib`,
`samtools`, `pybigwig`.

### Reference & supplement setup

The workflow needs per-genome supplement files (chrom sizes, `.2bit`, interval bins, blacklist, gap, and
a mappability bigWig). `scripts/setup_reference.sh` builds all of them, with the exact filenames the
example configs expect:

```bash
# T2T-CHM13: UCSC hs1 2bit/sizes + BEDbase excluderanges blacklist + Zenodo mappability track
scripts/setup_reference.sh t2t-chm13 supplement 500

# hg38: UCSC 2bit/sizes + Boyle-Lab Blacklist + finaletoolkit gap-bed + Zenodo mappability track
scripts/setup_reference.sh hg38 supplement 500
```

This writes into `supplement/`:

```
<g>.chrom.sizes   <g>.2bit   <g>.<N>kb.bins   <g>.delfi.chrom.sizes   <g>.blacklist.bed   <g>.45mer.mappability.bw
                  ( + hg38.gap.bed for hg38 only — T2T-CHM13 is gap-free )
```

(`g` = `hg38` or `chm13`). `<g>.delfi.chrom.sizes` is chr1–22, X, Y (DELFI requires
centromere-bearing contigs; the general bins keep chrM/chrY). Everything except the mappability bigWig is built live from public sources
(UCSC, Boyle-Lab Blacklist, BEDbase). The `*.bins.filtered` interval file used by DELFI is generated
automatically during the run by the mappability filter.

### Usage

1. Put input fragment/alignment files in `input/` (or set `input_dir`).
2. Pick a config and run:

```bash
snakemake --configfile params.t2t-chm13.yaml --cores <N> \
  --rerun-incomplete --default-resources "tmpdir='./tmp'"
```

`--cores` sets CPU cores; `--jobs` caps concurrent jobs. Use `-n` for a dry run to preview the DAG.

### SLURM execution

The SLURM executor plugin is already in `environment.yml`. Set your account/partition in
[`slurm_profile/config.yaml`](./slurm_profile/config.yaml), then submit:

```bash
./ftk_exc.sh params.t2t-chm13.yaml ./tmp     # runs in the background -> snakemake.log
```

### Mappability

Interval bins are kept only if their **mean** mappability over the bin is ≥ `mappability_threshold`
(filtered by `scripts/mappability_filter.py`). Because the filter uses the per-bin mean, the *kind* of
track matters:

- **Continuous track** (recommended) — each position = `1 / (k-mer occurrences)`, so the bin mean is the
  **average mappability** (a position that maps twice contributes 0.5). The shipped **hg38** and
  **T2T-CHM13** tracks are both continuous 45-mer GenMap tracks (`-E 0`), hosted on Zenodo
  ([record 20724659](https://zenodo.org/records/20724659)) and fetched automatically by
  `setup_reference.sh`. To build one for another genome / k-mer:

  ```bash
  conda create -n genmap -c bioconda -c conda-forge genmap ucsc-bedgraphtobigwig
  conda activate genmap
  scripts/genmap_mappability.sh genome.fa supplement/<g>.chrom.sizes \
      supplement/<g>.45mer.mappability.bw 45 0   # memory/time-heavy: use a compute node
  ```

T2T-CHM13 specifics: the assembly is gap-free and has no finaletoolkit gap preset, so DELFI runs
**without a gap file** (centromeres are already removed by the mappability + blacklist filtering); the
blacklist is the [`excluderanges`](https://dozmorovlab.github.io/excluderanges/) `T2T.excluderanges` set.
DELFI uses `<g>.delfi.chrom.sizes` (chr1–22, X, Y) because finaletoolkit crashes on contigs without a
centromere (e.g. chrM) when a gap file is supplied.

### Data availability

The continuous 45-mer mappability bigWigs (hg38, T2T-CHM13) are deposited in the **FinaleToolkit
Dataset** on Zenodo — DOI [10.5281/zenodo.20724659](https://doi.org/10.5281/zenodo.20724659)
(concept DOI [10.5281/zenodo.14284132](https://doi.org/10.5281/zenodo.14284132) always resolves to the
latest version). All other supplement files are built live from public sources (UCSC, Boyle-Lab
Blacklist, BEDbase) by `setup_reference.sh`.

### YAML parameters

* **Required:** `input_dir`, `output_dir`, `file_format` (`bed.gz`, `frag.gz`, `bam`, or `cram`).
* **Optional:** `supplement_dir`, `interval_file`, `mappability_file`, `mappability_threshold`, and any
  FinaleToolkit command.
* **Commands:** enable a FinaleToolkit command with underscores instead of hyphens (e.g. `adjust-wps`
  → `adjust_wps: True`); set flags by appending `_<flag>` (e.g. `coverage_mapq: 30`). The workflow
  respects command dependencies (e.g. `mds` needs `end_motifs`, `regional_mds` needs
  `interval_end_motifs`). The per-region motif-diversity step is `regional_mds`. See
  [`params.yaml`](./params.yaml) for the fully annotated list of every option.

### Output file naming

* **Filtered files** get `.filtered` before the format (e.g. `file.filtered.bed.gz`).
* **Command outputs** insert the command name (e.g. `file.frag_length_intervals.bed`).
* Each input is processed for every enabled command.

### Citation
Li JW, Bandaru R, Baliga K, Liu Y (2025). *FinaleToolkit: accelerating cell-free DNA fragmentation
analysis with a high-speed computational toolkit.* **Bioinformatics Advances**.
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioadv%2Fvbaf236-836eaa?style=flat-square)](https://doi.org/10.1093/bioadv/vbaf236)

### Contact
- Ravi Bandaru: ravi.bandaru@northwestern.edu
- James Li: james.li3@northwestern.edu
- Kundan Baliga: kundanbal2969@k12.ipsd.org
- Yaping Liu: yaping@northwestern.edu

### License
See [`LICENSE`](./LICENSE).

#!/bin/bash
snakemake --rerun-incomplete --profile slurm_profile > snakemake.log 2>&1 &
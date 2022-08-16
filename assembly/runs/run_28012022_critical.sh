#!/bin/bash
#SBATCH --job-name tr14uHMW_critical
#SBATCH --ntasks 8
#SBATCH --nodes 1
#SBATCH --mem 150GB
#SBATCH --time 20:00:00
#SBATCH --output slurm-tr14uHMW_critical.txt
#SBATCH -p critical

set -x
set -e
snakemake --rerun-incomplete -F --use-conda --cores 8  --snakefile assembly.smk --config  metadata=runs/metadata_tr14_28012022_test200GB.txt

#!/bin/bash
#SBATCH --job-name tr14uHMW_critical
#SBATCH --ntasks 8
#SBATCH --nodes 1
#SBATCH --mem 500GB
#SBATCH --time 20:00:00
#SBATCH --output slurm-tr14uHMW.txt

set -x
set -e
snakemake --rerun-incomplete -F --use-conda --cores 8  --snakefile assembly.smk --config  metadata=runs/metadata_tr14_28012022.txt

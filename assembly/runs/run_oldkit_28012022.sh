#!/bin/bash
#SBATCH --job-name=tr14HMW_oldkit
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH -p long
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-tr14HMW_oldkit.txt

set -x
set -e
snakemake --rerun-incomplete -F --use-conda --cores 8  --snakefile assembly.smk --config metadata=runs/metadata_tr14_oldkit_28012022.txt

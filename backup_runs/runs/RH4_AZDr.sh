#!/bin/bash
#SBATCH --job-name=RH4_AZDr
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH -p long
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-RH4_AZDr.txt

set -x

date
hostname

snakemake -F --jobs 8 --max-jobs-per-second 10 --keep-going --rerun-incomplete --printshellcmds --config metafile=runs/RH4_AZDr.txt tmp_dir=/fast/users/giurgium_c/scratch/RH4_AZDr working_dir=/fast/users/giurgium_c/CircleSeq/Nanopore project_name=Celllines --unlock

snakemake -F --jobs 8 --max-jobs-per-second 10 --keep-going --rerun-incomplete --printshellcmds --config metafile=runs/RH4_AZDr.txt tmp_dir=/fast/users/giurgium_c/scratch/RH4_AZDr working_dir=/fast/users/giurgium_c/CircleSeq/Nanopore project_name=Celllines 

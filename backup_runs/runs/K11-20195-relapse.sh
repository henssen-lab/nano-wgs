#!/bin/bash
#SBATCH --job-name=K11-20195-relapse
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH -p long
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-K11-20195-relapse.txt

set -x

date
hostname

snakemake -F --jobs 8 --max-jobs-per-second 10 --keep-going --rerun-incomplete --printshellcmds --config metafile=runs/K11-20195-relapse.txt tmp_dir=/fast/users/giurgium_c/scratch/K11-20195-relapse working_dir=/fast/users/giurgium_c/CircleSeq/Nanopore project_name=Evolution --unlock

snakemake -F --jobs 8 --max-jobs-per-second 10 --keep-going --rerun-incomplete --printshellcmds --config metafile=runs/K11-20195-relapse.txt tmp_dir=/fast/users/giurgium_c/scratch/K11-20195-relapse working_dir=/fast/users/giurgium_c/CircleSeq/Nanopore project_name=Evolution 

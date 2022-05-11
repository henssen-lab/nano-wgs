#!/bin/bash
#SBATCH --job-name=CB1008-initial
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH -p long
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-CB1008-initial.txt

set -x

date
hostname

snakemake -F --jobs 8 --max-jobs-per-second 10 --keep-going --rerun-incomplete --printshellcmds --config metafile=runs/CB1008-initial.txt tmp_dir=/fast/users/giurgium_c/scratch/CB1008-initial working_dir=/fast/users/giurgium_c/CircleSeq/Nanopore project_name=Evolution --unlock

snakemake -F --jobs 8 --max-jobs-per-second 10 --keep-going --rerun-incomplete --printshellcmds --config metafile=runs/CB1008-initial.txt tmp_dir=/fast/users/giurgium_c/scratch/CB1008-initial working_dir=/fast/users/giurgium_c/CircleSeq/Nanopore project_name=Evolution 

#!/bin/bash
#SBATCH --job-name=K11-20195-relapse
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH -p long
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-K11-20195-relapse_filter.txt

set -x

date
hostname

snakemake --snakefile workflows/postprocess.smk -F --jobs 8 --max-jobs-per-second 10 --keep-going --printshellcmds --config metafile=runs/K11-20195-relapse_filter.txt tmp_dir=/fast/users/giurgium_c/scratch/K11-20195-relapse_filter working_dir=/fast/users/giurgium_c/CircleSeq/Nanopore/ project_name=Evolution --unlock

snakemake --snakefile workflows/postprocess.smk -F --jobs 8 --max-jobs-per-second 10 --keep-going --printshellcmds --config metafile=runs/K11-20195-relapse_filter.txt tmp_dir=/fast/users/giurgium_c/scratch/K11-20195-relapse_filter working_dir=/fast/users/giurgium_c/CircleSeq/Nanopore/ project_name=Evolution 
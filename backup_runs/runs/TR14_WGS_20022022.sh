#!/bin/bash
#SBATCH --job-name=TR14_WGS_20022022
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH -p long
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-TR14_WGS_20022022.txt

set -x

date
hostname

meta=../runs/TR14_WGS_20022022.txt
tmp=/fast/users/giurgium_c/scratch/TR14_WGS_20022022
working_dir=/fast/users/giurgium_c/Nanopore/Celllines
project_name=Process

# go to the root of snakemake pipeline
cd nano-wgs

snakemake \
	--jobs 8 \
	--max-jobs-per-second 10 \
	--keep-going --rerun-incomplete --printshellcmds \
	--config metafile=${meta} \
	tmp_dir=${tmp} \
	working_dir=${working_dir} \
	project_name=${project_name} --unlock

snakemake \
	--jobs 8 \
        --max-jobs-per-second 10 \
        --keep-going --rerun-incomplete --printshellcmds \
        --config metafile=${meta} \
        tmp_dir=${tmp} \
        working_dir=${working_dir} \
        project_name=${project_name}


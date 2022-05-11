#!/bin/bash
#SBATCH --job-name=CHP212-DNA
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH -p long
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-CHP212-DNA.txt

set -x

date
hostname

meta=../runs/CHP212_WGS_08032022.txt
tmp=/fast/users/giurgium_c/scratch/CHP212-DNA
working_dir=/fast/users/giurgium_c/Nanopore
project_name=Celllines

# go to the root of snakemake pipeline
cd nano-wgs

snakemake -F \
	--jobs 8 \
	--max-jobs-per-second 10 \
	--keep-going --rerun-incomplete --printshellcmds \
	--config metafile=${meta} \
	tmp_dir=${tmp} \
	working_dir=${working_dir} \
	project_name=${project_name} --unlock

snakemake -F \
	--jobs 8 \
        --max-jobs-per-second 10 \
        --keep-going --rerun-incomplete --printshellcmds \
        --config metafile=${meta} \
        tmp_dir=${tmp} \
        working_dir=${working_dir} \
        project_name=${project_name}


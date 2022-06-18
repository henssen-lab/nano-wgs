#!/bin/bash
#SBATCH --job-name=TR14-uHMW-DNA
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-TR14-uHMW-DNA.txt

set -x

date
hostname

meta=../runs/TR14-uHMW-DNA.txt
tmp=/fast/users/giurgium_c/scratch/TR14-uHMW-DNA
working_dir=/fast/users/giurgium_c/Nanopore
project_name=Celllines
refs="hg38"

# go to the root of snakemake pipeline
cd nano-wgs

snakemake -F --cores 32 \
	--keep-going --rerun-incomplete --printshellcmds \
	--config metafile=${meta} \
	tmp_dir=${tmp} \
	threads=64 \
	read_length=5000 \
	working_dir=${working_dir} \
	project_name=${project_name} --unlock

snakemake -F --cores 32 \
        --keep-going --rerun-incomplete --printshellcmds \
        --config metafile=${meta} \
	refs=${refs} \
	threads=64 \
	read_length=5000 \
	tmp_dir=${tmp} \
        working_dir=${working_dir} \
        project_name=${project_name}


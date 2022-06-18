#!/bin/bash
#SBATCH --job-name=TR14-uHMW-DNA-hg19
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-TR14-uHMW-DNA-hg19.txt

set -x

date
hostname

meta=../runs/TR14-uHMW-DNA_hg19.txt
tmp=/fast/users/giurgium_c/scratch/TR14-uHMW-DNA-hg19
working_dir=/fast/users/giurgium_c/Nanopore/Celllines
project_name=Process
refs="hg19"

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


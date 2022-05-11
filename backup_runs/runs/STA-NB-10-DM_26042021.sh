#!/bin/bash
#SBATCH --job-name=STA-10-DM
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-STA-NB-10-DM.txt

set -x

date
hostname

meta=../runs/STA-NB-10-DM_26042021.txt
tmp=/fast/users/giurgium_c/scratch/STA-NB-10-DM
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
	read_length=300 \
	working_dir=${working_dir} \
	project_name=${project_name} --unlock

snakemake -F --cores 32 \
        --keep-going --rerun-incomplete --printshellcmds \
        --config metafile=${meta} \
	refs=${refs} \
	threads=64 \
	read_length=300 \
	tmp_dir=${tmp} \
        working_dir=${working_dir} \
        project_name=${project_name}


#!/bin/bash
#SBATCH --job-name=Rocio-hg38
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-Rocio-Multiplexed-TR14-CHP212-hg38.txt

set -x

date
hostname

meta=../runs/Rocio_Multiplexed.txt
mux=../runs/Rocio_Multiplexed_Barcodes.txt
tmp=/fast/users/giurgium_c/scratch/Rocio-Multiplexed-TR14-CHP212-hg38
working_dir=/fast/users/giurgium_c/Nanopore
project_name=Celllines
refs="hg38"

# go to the root of snakemake pipeline
cd nano-wgs

snakemake -F --cores 4 \
	--keep-going --rerun-incomplete --printshellcmds \
	--config metafile=${meta} \
	muxfile=${mux} \
	tmp_dir=${tmp} \
	threads=8 \
	read_length=300 \
	working_dir=${working_dir} \
	project_name=${project_name} --unlock

snakemake -F --cores 4 \
        --keep-going --rerun-incomplete --printshellcmds \
        --config metafile=${meta} \
	muxfile=${mux} \
	refs=${refs} \
	threads=8 \
	read_length=300 \
	tmp_dir=${tmp} \
        working_dir=${working_dir} \
        project_name=${project_name}


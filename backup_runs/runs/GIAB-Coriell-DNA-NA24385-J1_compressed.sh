#!/bin/bash
#SBATCH --job-name=GIAB-Coriell-DNA-NA24385-J1_compressed
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH -p long
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-GIAB-Coriell-DNA-NA24385-J1_compressed.txt

set -x

date
hostname

snakemake -F --jobs 8 --max-jobs-per-second 10 --keep-going --rerun-incomplete --printshellcmds --config metafile=runs/GIAB-Coriell-DNA-NA24385-J1_compressed.txt tmp_dir=/fast/users/giurgium_c/scratch/GIAB-Coriell-DNA-NA24385-J1_compressed working_dir=/fast/users/giurgium_c/scratch project_name=Process --unlock

snakemake -F --jobs 8 --max-jobs-per-second 10 --keep-going --rerun-incomplete --printshellcmds --config metafile=runs/GIAB-Coriell-DNA-NA24385-J1_compressed.txt tmp_dir=/fast/users/giurgium_c/scratch/GIAB-Coriell-DNA-NA24385-J1_compressed working_dir=/fast/users/giurgium_c/scratch project_name=Process 

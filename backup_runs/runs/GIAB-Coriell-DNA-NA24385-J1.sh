#!/bin/bash
#SBATCH --job-name=GIAB-Coriell-DNA-NA24385-J1
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH -p long
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-GIAB-Coriell-DNA-NA24385-J1.txt

set -x

date
hostname

snakemake -F --jobs 8 --max-jobs-per-second 10 --keep-going --rerun-incomplete --printshellcmds --config metafile=runs/GIAB-Coriell-DNA-NA24385-J1.txt tmp_dir=/fast/users/giurgium_c/scratch/GIAB-Coriell-DNA-NA24385-J1 working_dir=/fast/users/giurgium_c/scratch project_name=Process --unlock

snakemake -F --jobs 8 --max-jobs-per-second 10 --keep-going --rerun-incomplete --printshellcmds --config metafile=runs/GIAB-Coriell-DNA-NA24385-J1.txt tmp_dir=/fast/users/giurgium_c/scratch/GIAB-Coriell-DNA-NA24385-J1 working_dir=/fast/users/giurgium_c/scratch project_name=Process 

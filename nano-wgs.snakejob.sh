#!/bin/bash
#SBATCH --job-name=nanowgs
#SBATCH --output=log.txt
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --time=04-00
#SBATCH --mem-per-cpu=100M
#SBATCH -p medium

rm slurm-* # horrible quick fix
set -x

date
hostname
>&2 echo "Running Snakemake..."
snakemake \
    --use-conda \
    --jobs 20 \
    --max-jobs-per-second 10 \
    --cluster-config config_slurm.yaml \
    --drmaa " --partition={cluster.partition} --time={cluster.time} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu}" \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds

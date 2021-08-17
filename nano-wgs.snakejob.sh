#!/bin/bash
#SBATCH --job-name=nanowgs
#SBATCH --output=slurm.log
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --time=04-00
#SBATCH --mem-per-cpu=1000M
#SBATCH -p medium
#SBATCH --output=slurm-log.txt

mkdir -p slurm-out && rm -r slurm-out && mkdir -p slurm-out
set -x

date
hostname
>&2 echo "Running Snakemake..."
snakemake \
    --use-conda \
    --jobs 32 \
    --max-jobs-per-second 10 \
    --cluster-config config_slurm.yaml \
    --drmaa " --partition={cluster.partition} --time={cluster.time} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} --output=slurm-out/slurm-%j.out" \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds

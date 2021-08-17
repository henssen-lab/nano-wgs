#!/bin/bash
#SBATCH --job-name=nanowgs
#SBATCH --output=slurm.log
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000M
#SBATCH -p long
#SBATCH --output=slurm-log.txt

mkdir -p slurm-out && rm -r slurm-out && mkdir -p slurm-out
set -x

date
hostname
>&2 echo "Running Snakemake..."
snakemake \
    --use-conda \
    --jobs 4 \
    --max-jobs-per-second 10 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds

# --cluster-config config_slurm.yaml \
# --drmaa " --partition={cluster.partition} --time={cluster.time} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} --output=slurm-out/slurm-%j.out" \


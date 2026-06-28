#!/bin/bash

# metadata_patients_rerun_28062026.txt
meta=$1
working_dir=$2
project_name=$3
header=`head -1 $meta`

mkdir -p runs/

tail -n+2 $meta | awk '{print $3}' | sort | uniq | while read sample
do

echo "$header" > runs/${sample}.txt
grep $sample $meta >> runs/${sample}.txt


cat runs/${sample}.txt


echo """#!/bin/bash
#SBATCH --job-name=81mapping_$sample
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=150GB
#SBATCH --partition=highmem
#SBATCH --time 7-00:00:00
#SBATCH --output=logs/slurm-${sample}.txt

set -x

date
hostname
""" > runs/${sample}.sh

echo """snakemake \
-F \
--cores 8 \
--max-jobs-per-second 10 \
--keep-going \
--rerun-incomplete \
--printshellcmds \
--config threads=8 read_length=300 refs=hg38 metafile=runs/${sample}.txt tmp_dir=/data/cephfs-1/scratch/projects/henssen-nanopore/madalina/${sample} working_dir=${working_dir} project_name=${project_name} \
--unlock

snakemake \
-F \
--cores 8 \
--max-jobs-per-second 10 \
--keep-going \
--rerun-incomplete \
--printshellcmds \
--config threads=8 read_length=300 refs=hg38 metafile=runs/${sample}.txt tmp_dir=/data/cephfs-1/scratch/projects/henssen-nanopore/madalina/${sample} working_dir=${working_dir} project_name=${project_name} """ >> runs/${sample}.sh

sbatch runs/${sample}.sh

done





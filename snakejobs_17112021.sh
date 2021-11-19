#!/bin/bash

# metadata_12102021.txt
meta=$1
working_dir=$2
project_name=$3
header=`head -1 $meta`

mkdir -p runs/

tail -n+2 $meta | awk '{print $3}' | sort | uniq | while read sample
do

echo "$header" > runs/${sample}.txt
grep $sample $meta >> runs/${sample}.txt

echo """#!/bin/bash
#SBATCH --job-name=$sample
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH -p long
#SBATCH --time 7-00:00:00
#SBATCH --output=slurm-${sample}.txt

set -x

date
hostname
""" > runs/${sample}.sh

echo """snakemake \
-F \
--jobs 8 \
--max-jobs-per-second 10 \
--keep-going \
--rerun-incomplete \
--printshellcmds \
--config metafile=runs/${sample}.txt tmp_dir=/fast/users/giurgium_c/scratch/${sample} working_dir=${working_dir} project_name=${project_name} \
--unlock

snakemake \
-F \
--jobs 8 \
--max-jobs-per-second 10 \
--keep-going \
--rerun-incomplete \
--printshellcmds \
--config metafile=runs/${sample}.txt tmp_dir=/fast/users/giurgium_c/scratch/${sample} working_dir=${working_dir} project_name=${project_name} """ >> runs/${sample}.sh

sbatch runs/${sample}.sh

done





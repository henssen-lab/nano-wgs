#!/bin/bash
#SBATCH --nodes 1
#SBATCH --mem 200G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --time 80:00:00

set -x
set -e

root_dir=$1
fast5_dir=$2

cd ${root_dir}
compress_fast5 --input_path ${fast5_dir} --in_place --compression vbz --recursive --threads 32


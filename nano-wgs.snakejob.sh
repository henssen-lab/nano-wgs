#!/bin/bash
# Join stderr and stdout log files into stdout log file
#$ -j y
# Keep current environment variables, keeps virtualenv
#$ -V
# Use current working directory as working directory of the job.  This
# requires you to cd/popd into the directory where the snake_job.sh lies.
#$ -cwd
# Send a mail upon job completion and error
#$ -m eas
#$ -M max.mustermann(at)charite.de
#$ -P medium
#$ -l h_rt=96:00:00
#
#$ -o sge_log

set -x

# Enforce existence of TMPDIR -----------------------------------------------

export TMPDIR=/fast/users/${USER}/scratch/tmp
mkdir -p ${TMPDIR}

# Create one log directory per Snakemake run --------------------------------

touch sge_log
mkdir -p sge_log_dir
test -z "${JOB_ID}" && JOB_ID=$(date +%Y-%m-%d_%H-%M)
LOGDIR=sge_log_dir/${JOB_ID}
mkdir -p ${LOGDIR}

# Activate bash cmd printing, debug info ------------------------------------

set -x
>&2 hostname
>&2 date

# Kick off Snakemake --------------------------------------------------------

MAX_JOBS=30
snakemake \
    --use-conda \
    -j ${MAX_JOBS} \
    --max-jobs-per-second 10 \
    --cluster-config config_sge.yaml \
    --drmaa " -V -l h_vmem={cluster.h_vmem} -l h_rt={cluster.h_rt} -pe {cluster.pe} -P {cluster.project} -j yes -o ${PWD}/${LOGDIR}" \
    -k \
    -p \
    $*

# Print date after finishing, for good measure ------------------------------

>&2 date
>&2 echo "All done. Have a nice day."

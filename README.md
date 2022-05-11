### Install dependencies


### Run pipeline


#### 1. Data structure

The input of the pipeline are `fastq` files. 
If you have only the `fast5` and not the `fastq` you need to perform the basecalling.
Usually for the nanopore projects the basecalling is already performed (currently using Guppy5).

Set up your data as follows:

```
<Project Directory>
├── metadata.txt
└── RawData/
    ├── <Run Directory 1>
        └── ...
            └── fast5/
                ├── <fast5 file 1>.fast5
                ├── ...
                └── <fast5 file n>.fast5
            └── fastq/
                ├── <fastq file 1>.fastq
                ├── ...
                └── <fastq file n>.fastq
    ├── ...
    └── <Run Directory n>
        └── ...
```

#### 2. Input data

To run the pipeline you need to include a `metadata.txt` file.<br/>
If you have a multiplexed run then you need to include also a `barcodes.txt` file.<br/>
For examples please check `backup_runs` folder.<br/>

#### 3. Configuration file

Please check the configurable parameters of the snakemake pipeline in `configs/config_run.yaml`.<br/>
You dont need to change this file because snakemake allows to override the configuration parameters during runtime.<br/>
For example, if you want to run the pipeline with only `hg38` then use:

```
snakemake <other_paramters> --config refs=hg38
```

The default minimal read length is `300 bp`.<br/>

#### 4. Run snakemake pipeline

Run pipeline using the structure of the following script:

```
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

# metadata.txt - configuration file
meta=backup_runs/runs/Rocio_Multiplexed.txt 

# barcodes.txt - barcodes specification in case of multiplexing
mux=backup_runs/runs/Rocio_Multiplexed_Barcodes.txt

# set temp folder
tmp=/fast/users/giurgium_c/scratch/Rocio-Multiplexed-TR14-CHP212-hg38
working_dir=/fast/users/giurgium_c/Nanopore
project_name=Celllines
refs="hg38"

# go to the root of snakemake pipeline

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

```

### Install dependencies

1. Make sure you have conda installed.

2. Install guppy manually.

3. Then, download nano-wgs directory and install dependencies using conda.

```
git -clone https://github.com/k0n5/nano-wgs.git
cd nano-wgs
conda env create -f nano-wgs-env.yaml
```

### Run pipeline

Set up your data as follows:

```
<Project Directory>
├── metadata.csv
└── RawData/
    ├── <Run Directory 1>
        └── ...
            └── fast5/
                ├── <fast5 file 1>.fast5
                ├── ...
                └── <fast5 file n>.fast5
    ├── ...
    └── <Run Directory n>
        └── ...
```

`metadata.csv` is a colon-separated file specifying the different Nanopore runs, library kits, flowcell versions, barcodes and samples. Entries in the `Run` column have to match directory names in `RawData/`. There can be only one `Kit` and `Flowcell` per run. `Kit` and `Flowcell` have  Make sure If no barcoding was used, use `none`. Samples with identical names will be merged across runs.

Example:
```
Run;Kit;Flowcell;Barcode;Sample
RMS_1;SQK-RBK004;FLO-MIN106;barcode01;Rh30
RMS_1;SQK-RBK004;FLO-MIN106;barcode03;Rh36
RMS_1;SQK-RBK004;FLO-MIN106;barcode04;TE441T
RMS_1;SQK-RPK004;FLO-MIN106;barcode05;T174
HMW_IMR575_run1;SQK-LSK109;FLO-MIN106;none;IMR-5-75
HMW_IMR575_run2;SQK-LSK109;FLO-MIN106;none;IMR-5-75
```

Then, edit the header of the `<PATH TO>/nano-wgs/Snakefile` to match your setup:
```
# PROJECT INFORMATION
PROJECT_NAME = "<CHOOSE A PROJECT NAME HERE>"
WORKING_DIR = "<PATH TO PROJECT DIRECTORY>"

# GENERAL
GUPPY_BASECALLER = "<PATH TO>/ont-guppy-cpu/bin/guppy_basecaller"
HG19 = "<PATH TO>/hg19.fa"
```

Once this is set up, you are ready to go:

```
cd <PATH TO>/nano-wgs/
conda activate nano-wgs-env
qsub nano-wgs.snakejob.sh
```

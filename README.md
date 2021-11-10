### Install dependencies


### Run pipeline

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


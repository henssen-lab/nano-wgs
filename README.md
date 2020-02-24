### Setup

1. Make sure you have conda installed.

2. Install guppy manually.

3. Then, download nano-wgs directory and install dependencies using conda.
```
git -clone https://github.com/k0n5/nano-wgs.git
cd nano-wgs
```

### Run pipeline

Set up your data as follows:
```
<PROJECT DIRECTORY NAME>
├── metadata.csv
└── RawData/
    ├── <FIRST RUN DIRECTORY>
    ├── <...>
    └── <LAST RUN DIRECTORY>
```

metadata.csv is a file containing

```
conda activate nano-wgs-env
qsub nano-wgs.snakejob.sh
```

# Perform CNV calling

Assumptions. The data was mapped to hg19/hs37d5 using ngmlr.

First, clone this repo locally:

```

git clone https://github.com/henssen-lab/nano-wgs.git

```

Then go to `cnv` folder and clone the `smurfseq` project:

```

cd nano-wgs/cnv
git clone https://github.com/smithlabcode/smurfseq_scripts.git

```

To run the CNV calling:

```
sbatch run_cnv.sh ${folder_to_bamfile} ${samplename}
```
This will create a `cnv` folder in the `${folder_to_bamfile}` where the CNV profile and the raw segmental data are stored. 

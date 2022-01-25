import os
import pandas as pd
import sys
from snakemake.utils import min_version

min_version("6.4.1")

configfile: "../configs/config_assembly.yaml"

OUTDIR = config["outputdir"]
MINLEN = config["filtering"]["assembly"]["readlen"]
THREADS = config["threads"]
METADATA = pd.read_csv(config["metadata"],header=0,sep="\t")
REF = config["hg38"]["genome"]
GTF = config["hg38"]["gtf"]
BWA_ALING = config["params"]["mapping"]

rule all:
    input:
        expand(["{outdir}/{sample}/assembly/Assembly.fasta",
                "{outdir}/{sample}/assembly/assembly_bwa.sam"],outdir=OUTDIR,sample=METADATA.Sample.drop_duplicates().tolist())

rule merge_fastq:
    input:
        fastq=lambda wildcards: ["{run}".format(run=row.Run) for index, row in
                                 METADATA[METADATA.Sample == wildcards.sample].iterrows()],
    output:
        allfastq=temp("{outdir}/{sample}/assembly/all.fastq")
    shell:
        """mkdir -p {wildcards.outdir}/{wildcards.sample}/assembly && find {input.fastq} -name '*.fastq' | xargs cat > {output.allfastq}"""

rule assembly:
    input:
        fastq="{outdir}/{sample}/assembly/all.fastq",
    output:
        dir=directory("{outdir}/{sample}/assembly"),
        f="{outdir}/{sample}/assembly/Assembly.fasta"
    threads: THREADS
    params:
        minlen=MINLEN
    conda:
        "../envs/assembly-env.yaml"
    log:
        "{outdir}/{sample}/logs/assembly.log"
    shell:
        """
        mkdir -p {output.dir}
        shasta \
            --threads {threads}\
            --input {input.fastq} \
            --config Nanopore-Oct2021 \
            --Reads.minReadLength {params.minlen} \
            --assemblyDirectory {output.dir} &> {log}
        """

rule assembly_qc:
    input:
        assembly="{outdir}/{sample}/assembly/Assembly.fasta",
        ref=REF,
        gtf=GTF
    output:
        directory("{outdir}/{sample}/assembly_qc")
    threads: THREADS
    conda:
        "../envs/assembly-env.yaml"
    log:
        "{outdir}/logs/{sample}/assembly_qc.log"
    shell:
        """
        quast {input.assembly} \
            -r {input.ref} \
            -g {input.gtf} \
            -o {output} &> {log}
        """

rule assembly_against_ref:
    input:
        assembly="{outdir}/{sample}/assembly/Assembly.fasta",
        ref=REF,
        gtf=GTF
    output:
        sam="{outdir}/{sample}/assembly/assembly_bwa.sam"
    threads: THREADS
    params:
        BWA_ALING
    conda:
        "../envs/assembly-env.yaml"
    shell:
        """
        bwa mem {params} -t {threads} {input.ref} {input.assembly} | samtools sort -o {output.sam} - 
        """

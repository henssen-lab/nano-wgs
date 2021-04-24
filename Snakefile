# Installation:
# install guppy manually
# conda install pandas numpy
# conda install -c bioconda samtools=1.9 nanoplotter minimap2 ngmlr sniffles deeptools qcat bioconductor-qdnaseq bioconductor-qdnaseq.hg19
# conda env export --from-history > nano-wgs-env.yaml
# samtools must be of version 1.9! There are issues with long reads and 1.10

import os
import sys
import pandas as pd

# PROJECT INFORMATION
PROJECT_NAME = "BER05_NanoporeCircleSeq"
WORKING_DIR = "/fast/groups/ag_henssen/work/Kons_Nanopore_circleseq_BER05pool/"

# GENERAL
GUPPY_BASECALLER = "/fast/users/helmsauk_c/work/ont-guppy-cpu/bin/guppy_basecaller"
HG19 = "/fast/users/helmsauk_c/work/resources/hg19_bwa/hg19.fa" # should but need not be a writable directory

# TODO: implement test that each run has exactly one kit
# TODO: implement test that each run has exactly one flowcell
# TODO: implement test that each barcode has one sample for each run

metadata = pd.read_csv(WORKING_DIR + "metadata.csv", sep=";", comment="#")
runs = list(set(metadata.Run.tolist()))
samples = list(set(metadata.Sample.tolist()))
kit_dict = pd.Series(metadata.Kit.values,index=metadata.Run).to_dict()
flowcell_dict = pd.Series(metadata.Flowcell.values,index=metadata.Run).to_dict()

rule all:
    input:
        expand(WORKING_DIR + "RunsQC/{run}-QC/{run}-NanoPlot-report.html", run = runs),
        expand(WORKING_DIR + "RunsQC/{run}-QC/{run}-demux.log", run = runs),
        expand(WORKING_DIR + "Samples/{sample}.fastq.gz", sample = samples),
        expand(WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bam", sample = samples),
        expand(WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bam.bai", sample = samples),
        expand(WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bw", sample = samples),
        expand(WORKING_DIR + "Samples/{sample}.ngmlr_hg19.stats.txt", sample = samples),
        expand(WORKING_DIR + "Samples/{sample}.ngmlr_hg19.sniffles.vcf", sample = samples),
        expand(WORKING_DIR + "Samples/{sample}.ngmlr_hg19.svim.vcf", sample = samples),
        #expand(WORKING_DIR + "Samples/{sample}.ngmlr_hg19.copynumber.pdf", sample = samples),
        WORKING_DIR + "SamplesQC/" + PROJECT_NAME + "-NanoComp-report.html",

rule basecalling:
    output:
        temp(WORKING_DIR + "Runs/{run}/{run}.fastq"),
        temp(WORKING_DIR + "Runs/{run}/sequencing_summary.txt"),
        temp(WORKING_DIR + "Runs/{run}/sequencing_telemetry.js")
    params:
        working_dir = WORKING_DIR,
        output_dir = WORKING_DIR + "Runs/{run}/",
        run = "{run}",
        kit = lambda wildcards: kit_dict[wildcards.run],
        flowcell = lambda wildcards: flowcell_dict[wildcards.run]
    shell:
       GUPPY_BASECALLER + " --num_callers 24 --input_path {params.working_dir}RawData/{params.run}/ --save_path {params.output_dir} --flowcell {params.flowcell} --kit {params.kit} --recursive -- && cat {params.output_dir}*.fastq > {params.output_dir}{params.run}.fastq && rm {params.output_dir}fastq_runid*.fastq && rm {params.output_dir}*.log"

rule nanoplot:
    input:
        WORKING_DIR + "Runs/{run}/sequencing_summary.txt"
    output:
        WORKING_DIR + "RunsQC/{run}-QC/{run}-NanoPlot-report.html"
    params:
        run = "{run}",
        output_dir = WORKING_DIR + "RunsQC/{run}-QC",
    shell:
        "NanoPlot --summary {input} --loglength --outdir {params.output_dir} --format pdf --N50 --title {params.run} --prefix {params.run}-"

rule demultiplex:
    input:
        WORKING_DIR + "Runs/{run}/{run}.fastq"
    output:
        log_file = WORKING_DIR + "RunsQC/{run}-QC/{run}-demux.log"
    params:
        run = "{run}",
        output_dir = WORKING_DIR + "Runs/{run}/"
    run:
        shell("qcat --trim --detect-middle --kit Auto --fastq {input} --barcode_dir {params.output_dir} 2> {output.log_file}")
        for index, row in metadata[metadata.Run == params.run].iterrows():
            shell("gzip -c {params.output_dir}" + str(row["Barcode"]) + ".fastq > {params.output_dir}" + str(row["Sample"]) + ".fastq.gz")
        shell("rm -f {params.output_dir}barcode*.fastq")
        shell("rm -f {params.output_dir}none.fastq")

rule merge_samples:
    input:
        lambda wildcards: [WORKING_DIR + "RunsQC/{run}-QC/{run}-demux.log".format(run=row.Run) for index, row in metadata[metadata.Sample == wildcards.sample].iterrows()]
    output:
        WORKING_DIR + "Samples/{sample}.fastq.gz"
    params:
        files = lambda wildcards: [WORKING_DIR + "Runs/{run}/{sample}.fastq.gz".format(run=row.Run, sample=row.Sample) for index, row in metadata[metadata.Sample == wildcards.sample].iterrows()]
    shell:
        "cat {params.files} > {output} && rm {params.files}"

rule ngmlr:
    input:
        WORKING_DIR + "Samples/{sample}.fastq.gz"
    output:
        bam = WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bam",
        bai = WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bam.bai"
    params:
        sam = WORKING_DIR + "Samples/{sample}.ngmlr_hg19.sam",
        reference = HG19,
        threads = 24
    shell:
        """
        ngmlr --bam-fix --threads {params.threads} --reference {params.reference} --query {input} --output {params.sam} --presets ont
        samtools sort --threads {params.threads} -o {output.bam} {params.sam}
        samtools index {output.bam}
        rm {params.sam}
        """

rule stats:
    input:
        WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bam"
    output:
        WORKING_DIR + "Samples/{sample}.ngmlr_hg19.stats.txt"
    shell:
        "samtools stats {input} > {output}"

rule nanocomp_samples:
    input:
        expand(WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bam", sample = samples)
    output:
        WORKING_DIR + "SamplesQC/" + PROJECT_NAME + "-NanoComp-report.html"
    params:
        name = PROJECT_NAME,
        output_dir = WORKING_DIR + "SamplesQC/",
        samples = expand("{s}", s=samples)
    shell:
        "NanoComp --bam {input} --names {params.samples} --outdir {params.output_dir} --format pdf --title {params.name} --prefix {params.name}-"

rule sniffles:
    input:
        WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bam"
    output:
        WORKING_DIR + "Samples/{sample}.ngmlr_hg19.sniffles.vcf"
    shell:
         "sniffles -m {input} -v {output}"

rule svim:
    input:
        WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bam"
    output:
        WORKING_DIR + "Samples/{sample}.ngmlr_hg19.svim.vcf"
    params:
        reference = HG19,
        svim_output_dir = WORKING_DIR + "Samples/{sample}-svim/",
        sample = "{sample}"
    shell:
         "svim alignment --read_names --insertion_sequences --sample {params.sample} {params.svim_output_dir} {input} {params.reference} && mv {params.svim_output_dir}variants.vcf {output}"

rule coverage:
    input:
        WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bam"
    output:
        WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bw"
    params:
        threads = 12
    shell:
        "bamCoverage --bam {input} -o {output} --binSize 20 -p {params.threads}"

rule copynumber:
    input:
        WORKING_DIR + "Samples/{sample}.ngmlr_hg19.bam"
    output:
        pdf = WORKING_DIR + "Samples/{sample}.ngmlr_hg19.copynumber.pdf",
        rdata = WORKING_DIR + "Samples/{sample}.ngmlr_hg19.copynumber.Rdata"
    shell:
        "Rscript cnprofile.R {input} {output.pdf} {output.rdata}"

onsuccess:
    shell("rm -rf " + WORKING_DIR + "Runs/")

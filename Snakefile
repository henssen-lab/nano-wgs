import os
import sys
import pandas as pd

configfile: "config_run.yaml"

PROJECT_NAME = config['PROJECT_NAME']
WORKING_DIR = config['WORKING_DIR']

# GENERAL
# GUPPY_BASECALLER = "/fast/users/helmsauk_c/work/ont-guppy-cpu/bin/guppy_basecaller"
HG19 = config['HG19']
HG38 = config['HG38']

MULTISAMPLE = config['MULTISAMPLE']

if MULTISAMPLE == True:
    metadata = pd.read_csv(WORKING_DIR + "metadata.csv", sep=";", comment="#")
    runs = list(set(metadata.Run.tolist()))
    samples = list(set(metadata.Sample.tolist()))
    kit_dict = pd.Series(metadata.Kit.values,index=metadata.Run).to_dict()
    flowcell_dict = pd.Series(metadata.Flowcell.values,index=metadata.Run).to_dict()
else:
    runs = [config['RUN']]
    samples = [config['SAMPLE']]
    summary = [config['SEQ_SUMMARY']]

os.chdir(os.path.join(WORKING_DIR, PROJECT_NAME, "Process"))

rule all:
    input:
        expand("{run}", run=runs)
        expand("{sample}/QC/NanoPlot-report.html"), sample = samples),
        expand("{sample}/ngmlr_hg19.bam"), sample = samples),
        expand("{sample}/ngmlr_hg19.bam.bai"), sample = samples),
        expand("{sample}/ngmlr_hg19.sniffles.vcf"), sample = samples),
        expand("{sample}/ngmlr_hg19.svim.vcf"), sample = samples),
        expand("{sample}/ngmlr_hg38.bam"), sample = samples),
        expand("{sample}/ngmlr_hg38.bam.bai"), sample = samples),
        expand("{sample}/ngmlr_hg38.sniffles.vcf"), sample = samples),
        expand("{sample}/ngmlr_hg38.svim.vcf"), sample = samples)
        # expand(os.path.join(WORKING_DIR, PROJECT_NAME, "Process", "{sample}", "QC", "NanoPlot-report.html"), sample = samples),


rule nanoplot:
    input:
        "{run}/sequencing_summary.txt"
    output:
        "{sample}/QC/NanoPlot-report.html"
    params:
        sample = "{sample}",
        output_dir = WORKING_DIR + "/" + PROJECT_NAME + "/Proces/" +  "{sample}" +  "/QC",
    shell:
        "NanoPlot --summary {input} --outdir {params.output_dir} --N50 --title {params.sample}"

rule nanfilt:
    input:
        fastqin="{run}/all.fastq"
    output:
        fastqout = temp("all.fastq")
    shell:
        "NanoFilt -l 1000 --headcrop 50 --tailcrop 50 --readtype 1D {input.fastqin} > {output.fastqout}"

rule ngmlr:
    input:
        fastq="{run}/all.fastq"
        reference=expand("{ref}", ref=[HG19, HG38])
    output:
        bam = expand("{sample}/ngmlr_{refid}.bam", refid=['hg19', 'hg38']),
        bai = expand("{sample}/ngmlr_{refid}.bam.bai", refid=['hg19', 'hg38']),
        sam = temp(expand("{sample}/ngmlr_{refid}.sam", refid=['hg19', 'hg38']))
    params:
        threads = 24
    shell:
        """
        ngmlr --bam-fix --threads {params.threads} --reference {input.reference} --query {input.fastq} --output {output.sam} --presets ont
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
         "sniffles -m {input} -v {output} --min_length 15 --genotype --min-support 3 --report-seq"

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

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

# if MULTISAMPLE == True:
#     metadata = pd.read_csv(WORKING_DIR + "metadata.csv", sep=";", comment="#")
#     runs = list(set(metadata.Run.tolist()))
#     samples = list(set(metadata.Sample.tolist()))
#     kit_dict = pd.Series(metadata.Kit.values,index=metadata.Run).to_dict()
#     flowcell_dict = pd.Series(metadata.Flowcell.values,index=metadata.Run).to_dict()
# else:
runs = [config['RUN']]
samples = [config['SAMPLE']]
# summary = [config['SEQ_SUMMARY']]

os.chdir(os.path.join(WORKING_DIR, PROJECT_NAME, "Process"))

rule all:
    input:
        expand("{run}", run=runs)
        expand("{sample}/QC/NanoPlot-report.html"), sample = samples),
        expand("{sample}/ngmlr_{refid}.bam"), sample = samples, refid=['hg19', 'hg38']),
        expand("{sample}/ngmlr_{refid}.bam.bai"), sample = samples, refid=['hg19', 'hg38']),
        expand("{sample}/ngmlr_{refid}.sniffles.vcf"), sample = samples, refid=['hg19', 'hg38']),
        expand("{sample}/ngmlr_{refid}.svim.vcf"), sample = samples, refid=['hg19', 'hg38']),
        expand("{sample}/ngmlr_{refid}.copynumber.Rdata"), sample = samples, refid=['hg19', 'hg38']),
        # expand("{sample}/ngmlr_hg38.bam"), sample = samples),
        # expand("{sample}/ngmlr_hg38.bam.bai"), sample = samples),
        # expand("{sample}/ngmlr_hg38.sniffles.vcf"), sample = samples),
        # expand("{sample}/ngmlr_hg38.svim.vcf"), sample = samples)

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
        fastqin = "{run}/all.fastq"
    output:
        fastqout = "all.fastq"
    params:
        rlen = config['read_length']
    shell:
        "NanoFilt -l {params.rlen} --headcrop 50 --tailcrop 50 --readtype 1D {input.fastqin} > {output.fastqout}"

rule ngmlr:
    input:
        fastq="all.fastq"
        reference=expand("{ref}", ref=[HG19, HG38])
    output:
        bam = expand("{sample}/ngmlr_{refid}.bam", refid=['hg19', 'hg38']),
        bai = expand("{sample}/ngmlr_{refid}.bam.bai", refid=['hg19', 'hg38']),
        sam = expand("{sample}/ngmlr_{refid}.sam", refid=['hg19', 'hg38'])
    params:
        threads = 32
    shell:
        """
        ngmlr --bam-fix --threads {params.threads} --reference {input.reference} --query {input.fastq} --output {output.sam} --presets ont
        samtools sort -O BAM -o {output.bam} {params.sam}
        samtools index {output.bam}
        """

rule stats:
    input:
        bam="{sample}/ngmlr_{refid}.bam"
        bai="{sample}/ngmlr_{refid}.bam.bai"
    output:
        "{sample}/ngmlr_{refid}.stats.txt"
    shell:
        "samtools stats {input.bam} > {output}"

rule sniffles:
    input:
        bam="{sample}/ngmlr_{refid}.bam"
        bai="{sample}/ngmlr_{refid}.bam.bai"
    output:
        "{sample}/ngmlr_{refid}.sniffles.vcf"
    params:
        threads = 32
    shell:
         "sniffles -t {params.threads} -m {input.bam} -v {output} --min_homo_af 0.7 --min_het_af 0.1 --min_length 500 --cluster --genotype --min-support 4 --report-seq"

rule svim:
    input:
        bam="{sample}/ngmlr_{refid}.bam"
        bai="{sample}/ngmlr_{refid}.bam.bai"
        reference="{ref}"
    output:
        "{sample}/ngmlr_{refid}.svim.vcf"
    params:
        svim_output_dir = "{sample}-svim/",
        sample = "{sample}"
    shell:
         "svim alignment --read_names --insertion_sequences --sample {params.sample} {params.svim_output_dir} {input.bam} {input.reference} && mv {params.svim_output_dir}variants.vcf {output}"

rule coverage:
    input:
        bam="{sample}/ngmlr_{refid}.bam"
        bai="{sample}/ngmlr_{refid}.bam.bai"
    output:
        "{sample}.ngmlr_{refid}.bw"
    params:
        threads = 32
    shell:
        "bamCoverage --bam {input.bam} -o {output} --binSize 20 -p {params.threads}"

rule copynumber:
    input:
        bam="{sample}/ngmlr_{refid}.bam"
        bai="{sample}/ngmlr_{refid}.bam.bai"
    output:
        pdf="{sample}/ngmlr_{refid}.copynumber.pdf",
        rdata = "{sample}/ngmlr_{refid}.copynumber.Rdata"
    shell:
        "Rscript cnprofile.R {input.bam} {output.pdf} {output.rdata}"

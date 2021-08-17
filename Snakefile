import os
import sys
import pandas as pd
from snakemake.utils import min_version 

min_version("6.4.1")

configfile: "config_run.yaml"

HG19 = 'hg19'
HG38 = 'hg38'
REFERENCE = 'reference'

PROJECT_NAME = config['project_name']
WORKING_DIR = config['working_dir']
METADATA_FILE = config['metafile']
TMP_DIR = config['tmp_dir']
HEADER_SUM = os.path.join(os.getcwd(), "envs/header_summary.txt")

metadata = pd.read_csv(METADATA_FILE, sep="\t", header=0)
runs = list(set(metadata.Run.tolist()))
samples = list(set(metadata.Sample.tolist()))
kit_dict = pd.Series(metadata.Kit.values,index=metadata.Run).to_dict()
flowcell_dict = pd.Series(metadata.Flowcell.values,index=metadata.Run).to_dict()

os.chdir(os.path.join(WORKING_DIR, PROJECT_NAME))

def get_reference(wildcards):
    """
    Get reference genome
    """
    if wildcards.refid == HG19:
        return config[REFERENCE][HG19]
    return config[REFERENCE][HG38]

rule all:
    input:
        expand(["Process/{sample}/QC/NanoPlot-report.html",
                "Process/{sample}/filt.fastq",
                "Process/{sample}/{refid}/ngmlr_{refid}.bam",
                "Process/{sample}/{refid}/ngmlr_{refid}.bam.bai",
                "Process/{sample}/{refid}/ngmlr_{refid}.stats.txt",
                "Process/{sample}/{refid}/ngmlr_{refid}.sniffles.vcf",
                "Process/{sample}/{refid}/ngmlr_{refid}.svim.vcf",
                "Process/{sample}/{refid}/coverage_{refid}.bw"
                ], sample=samples, refid=[HG19, HG38]),
        # expand("Process/{sample}/QC/NanoPlot-report.html", sample=samples)
        # expand("{run}/alltemp.fastq", run=[run]),
        # sample + "/QC/NanoPlot-report.html",
        # expand("{sample}/ngmlr_{refid}.bam", sample = [sample], refid=[HG19, HG38]),
        # expand("{sample}/ngmlr_{refid}.bam.bai", sample = [sample], refid=[HG19, HG38]),
        # expand("{sample}/ngmlr_{refid}.sniffles.vcf", sample = samples, refid=[HG19, HG38]),
        # expand("{sample}/ngmlr_{refid}.svim.vcf", sample = samples, refid=[HG19, HG38]),
        # expand("{sample}/ngmlr_{refid}.copynumber.Rdata", sample = samples, refid=[HG19, HG38]),
        # expand("{sample}/ngmlr_hg38.bam"), sample = samples),
        # expand("{sample}/ngmlr_hg38.bam.bai"), sample = samples),
        # expand("{sample}/ngmlr_hg38.sniffles.vcf"), sample = samples),
        # expand("{sample}/ngmlr_hg38.svim.vcf"), sample = samples)

        # expand(os.path.join(WORKING_DIR, PROJECT_NAME, "Process", "{sample}", "QC", "NanoPlot-report.html"), sample = samples),

rule merge_fastq:
    input:
        fastq = lambda wildcards: ["{run}/pass".format(run=row.Run) for index, row in metadata[metadata.Sample == wildcards.sample].iterrows()],
    output:
        allfastq = temp("Process/{sample}/all.fastq"),
    resources:
        tmpdir = TMP_DIR
    shell:
        """find {input.fastq} -name '*.fastq' | xargs cat > {output.allfastq}"""

rule merge_summary:
    input:
        sumfiles = lambda wildcards: ["{run}/sequencing_summary.txt".format(run=row.Run) for index, row in metadata[metadata.Sample == wildcards.sample].iterrows()],
    output:
        seqsum = temp("Process/{sample}/sequencing_summary.txt")
    resources:
        tmpdir = TMP_DIR
    params:
        shead = HEADER_SUM
    shell:
        """cat {params.shead} > {output.seqsum} && \
           cat {input.sumfiles} | grep -v "passes_filtering" >> {output.seqsum}"""

rule nanoplot:
    input:
        "Process/{sample}/sequencing_summary.txt"
    output:
        report = "Process/{sample}/QC/NanoPlot-report.html"
    resources:
        tmpdir = TMP_DIR
    conda:
        "envs/qc-env.yaml"
    params:
        sample = "{sample}",
        output_dir = "Process/{sample}/QC",
    shell:
        "NanoPlot --summary {input} --outdir {params.output_dir} --N50 --title {params.sample}"

rule nanfilt:
    input:
        fastqin = "Process/{sample}/all.fastq"
    output:
        fastqout = temp("Process/{sample}/filt.fastq")
    resources:
        tmpdir=TMP_DIR
    conda:
        "envs/mapping-env.yaml"
    params:
        rlen = config['read_length']
    shell:
        "NanoFilt -l {params.rlen} --headcrop 50 --tailcrop 50 --readtype 1D {input.fastqin} > {output.fastqout}"

rule ngmlr_mock:
    input:
        "Process/{sample}/filt.fastq"
    output:
        # outf = directory(expand("Process/{sample}/{refid}/", sample=["{sample}"], refid=[HG19, HG38])),
        outtemp = expand("Process/{sample}/{refid}/temp", sample=["{sample}"], refid=[HG19, HG38])
    resources:
        tmpdir = TMP_DIR
    conda:
        "envs/mapping-env.yaml"
    shell:
        "touch {output.outtemp}"

rule ngmlr:
    input:
        fastq = "Process/{sample}/filt.fastq",
        reference = get_reference,
        ftemp = "Process/{sample}/{refid}/temp"
    output:
        bam = protected("Process/{sample}/{refid}/ngmlr_{refid}.bam"),
        bai = protected("Process/{sample}/{refid}/ngmlr_{refid}.bam.bai"),
        sam = temp("Process/{sample}/{refid}/ngmlr_{refid}.sam")
    resources:
        tmpdir = TMP_DIR
    conda:
        "envs/mapping-env.yaml"
    params:
        threads = 32
    shell:
        """
        ngmlr --bam-fix --threads {params.threads} --reference {input.reference} --query {input.fastq} --output {output.sam} --presets ont
        samtools sort -O BAM -o {output.bam} {output.sam}
        samtools index {output.bam}
        """

rule coverage:
    input:
        bam="Process/{sample}/{refid}/ngmlr_{refid}.bam",
        bai="Process/{sample}/{refid}/ngmlr_{refid}.bam.bai"
    output:
        protected("Process/{sample}/{refid}/coverage_{refid}.bw")
    resources:
        tmpdir = TMP_DIR
    conda:
        "envs/mapping-env.yaml"
    params:
        threads = 32
    shell:
        "bamCoverage --bam {input.bam} -o {output} --binSize 50 -p {params.threads}"

rule stats:
    input:
        bam="Process/{sample}/{refid}/ngmlr_{refid}.bam",
        bai="Process/{sample}/{refid}/ngmlr_{refid}.bam.bai"
    output:
        "Process/{sample}/{refid}/ngmlr_{refid}.stats.txt"
    resources:
        tmpdir = TMP_DIR
    conda:
        "envs/mapping-env.yaml"
    shell:
        "samtools stats {input.bam} > {output}"

rule sniffles:
    input:
        bam="Process/{sample}/{refid}/ngmlr_{refid}.bam",
        bai="Process/{sample}/{refid}/ngmlr_{refid}.bam.bai"
    output:
        protected("Process/{sample}/{refid}/ngmlr_{refid}.sniffles.vcf")
    resources:
        tmpdir = TMP_DIR
    conda:
        "envs/snv-env.yaml"
    params:
        threads = 32
    shell:
         "sniffles -t {params.threads} -m {input.bam} -v {output} --min_homo_af 0.7 --min_het_af 0.1 --min_length 500 --cluster --genotype --min-support 4 --report-seq"

rule svim:
    input:
        bam="Process/{sample}/{refid}/ngmlr_{refid}.bam",
        bai="Process/{sample}/{refid}/ngmlr_{refid}.bam.bai",
        reference=get_reference
    output:
        protected("Process/{sample}/{refid}/ngmlr_{refid}.svim.vcf")
    resources:
        tmpdir = TMP_DIR
    conda:
        "envs/snv-env.yaml"
    params:
        svim_output_dir = "Process/{sample}/{refid}/svim/",
        sample = "{sample}"
    shell:
         "svim alignment --read_names --insertion_sequences --sample {params.sample} {params.svim_output_dir} {input.bam} {input.reference} && mv {params.svim_output_dir}variants.vcf {output}"


# rule copynumber:
#     input:
#         bam="{sample}/ngmlr_{refid}.bam",
#         bai="{sample}/ngmlr_{refid}.bam.bai"
#     output:
#         pdf="{sample}/ngmlr_{refid}.copynumber.pdf",
#         rdata = "{sample}/ngmlr_{refid}.copynumber.Rdata"
#     shell:
#         "Rscript cnprofile.R {input.bam} {output.pdf} {output.rdata}"

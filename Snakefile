import os
import sys
import pandas as pd

configfile: "config_run.yaml"

HG19 = 'hg19'
HG38 = 'hg38'
REFERENCE = 'reference'

PROJECT_NAME = config['project_name']
WORKING_DIR = config['working_dir']
METADATA_FILE = config['metafile']
TMP_DIR = config['tmp_dir']

metadata = pd.read_csv(METADATA_FILE, sep="\t", header=0)
runs = list(set(metadata.Run.tolist()))
samples = list(set(metadata.Sample.tolist()))
kit_dict = pd.Series(metadata.Kit.values,index=metadata.Run).to_dict()
flowcell_dict = pd.Series(metadata.Flowcell.values,index=metadata.Run).to_dict()

os.chdir(os.path.join(WORKING_DIR, PROJECT_NAME))
print(os.getcwd())


rule all:
    input:
        expand(["Process/{sample}/QC/NanoPlot-report.html",
                "Process/{sample}/filt.fastq"], sample=samples),
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
        allfastq = "Process/{sample}/all.fastq",
    resources:
        tmpdir = TMP_DIR
    shell:
        """find {input.fastq} -name '*.fastq' | xargs cat > {output.allfastq}"""

rule merge_summary:
    input:
        sumfiles = lambda wildcards: ["{run}/sequencing_summary.txt".format(run=row.Run) for index, row in metadata[metadata.Sample == wildcards.sample].iterrows()],
    output:
        seqsum = "Process/{sample}/sequencing_summary.txt"
    resources:
        tmpdir = TMP_DIR
    shell:
        """head -1 {input.sumfiles} | head -1 > {output.seqsum} && \
           tail -n+2 {input.sumfiles} >> {output.seqsum}"""

rule nanoplot:
    input:
        "Process/{sample}/sequencing_summary.txt"
    output:
        report = "Process/{sample}/QC/NanoPlot-report.html"
    resources:
        tmpdir = TMP_DIR
    params:
        sample = "{sample}",
        output_dir = "Process/{sample}/QC",
    shell:
        "NanoPlot --summary {input} --outdir {params.output_dir} --N50 --title {params.sample}"

rule nanfilt:
    input:
        fastqin = "Process/{sample}/all.fastq"
    output:
        fastqout = "Process/{sample}/filt.fastq"
    resources:
        tmpdir=TMP_DIR
    params:
        rlen = config['read_length']
    shell:
        "NanoFilt -l {params.rlen} --headcrop 50 --tailcrop 50 --readtype 1D {input.fastqin} > {output.fastqout}"

# rule ngmlr:
#     input:
#         fastq=run + "/alltemp.fastq",
#         reference=[config[REFERENCE][HG19], config[REFERENCE][HG38]],
#         refid=lambda wildcards: [HG19, HG38],
#         sample=lambda wildcards: sample
#     output:
#         bam = protected("{sample}/ngmlr_{refid}.bam"),
#         bai = protected("{sample}/ngmlr_{refid}.bam.bai"),
#         sam = temp("{sample}/ngmlr_{refid}.sam")
#     params:
#         threads = 32
#     shell:
#         """
#         ngmlr --bam-fix --threads {params.threads} --reference {input.reference} --query {input.fastq} --output {output.sam} --presets ont
#         samtools sort -O BAM -o {output.bam} {params.sam}
#         samtools index {output.bam}
#         """

# rule stats:
#     input:
#         bam="{sample}/ngmlr_{refid}.bam",
#         bai="{sample}/ngmlr_{refid}.bam.bai"
#     output:
#         "{sample}/ngmlr_{refid}.stats.txt"
#     shell:
#         "samtools stats {input.bam} > {output}"

# rule sniffles:
#     input:
#         bam="{sample}/ngmlr_{refid}.bam",
#         bai="{sample}/ngmlr_{refid}.bam.bai"
#     output:
#         "{sample}/ngmlr_{refid}.sniffles.vcf"
#     params:
#         threads = 32
#     shell:
#          "sniffles -t {params.threads} -m {input.bam} -v {output} --min_homo_af 0.7 --min_het_af 0.1 --min_length 500 --cluster --genotype --min-support 4 --report-seq"

# rule svim:
#     input:
#         bam="{sample}/ngmlr_{refid}.bam",
#         bai="{sample}/ngmlr_{refid}.bam.bai",
#         reference="{ref}"
#     output:
#         "{sample}/ngmlr_{refid}.svim.vcf"
#     params:
#         svim_output_dir = "{sample}-svim/",
#         sample = "{sample}"
#     shell:
#          "svim alignment --read_names --insertion_sequences --sample {params.sample} {params.svim_output_dir} {input.bam} {input.reference} && mv {params.svim_output_dir}variants.vcf {output}"

# rule coverage:
#     input:
#         bam="{sample}/ngmlr_{refid}.bam",
#         bai="{sample}/ngmlr_{refid}.bam.bai"
#     output:
#         "{sample}.ngmlr_{refid}.bw"
#     params:
#         threads = 32
#     shell:
#         "bamCoverage --bam {input.bam} -o {output} --binSize 20 -p {params.threads}"

# rule copynumber:
#     input:
#         bam="{sample}/ngmlr_{refid}.bam",
#         bai="{sample}/ngmlr_{refid}.bam.bai"
#     output:
#         pdf="{sample}/ngmlr_{refid}.copynumber.pdf",
#         rdata = "{sample}/ngmlr_{refid}.copynumber.Rdata"
#     shell:
#         "Rscript cnprofile.R {input.bam} {output.pdf} {output.rdata}"

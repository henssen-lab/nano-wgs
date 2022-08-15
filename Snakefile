import os

import pandas as pd
from snakemake.utils import min_version

min_version("6.4.1")

configfile: "configs/config_run.yaml"

HG19 = 'hg19'
HG38 = 'hg38'
T2T = 'T2T-CHM13'

REFERENCE = 'reference'
THREADS = config["threads"]
MUX = False

REFS = config['refs'].split(",")  # hg19,hg38
PROJECT_NAME = config['project_name']
WORKING_DIR = config['working_dir']
METADATA_FILE = config['metafile']
DEMUX_FILE = config['muxfile']
TMP_DIR = config['tmp_dir']
HEADER_SUM = os.path.join(os.getcwd(),"envs/header_summary.txt")

metadata = pd.read_csv(METADATA_FILE,sep="\t",header=0)

if len(DEMUX_FILE) > 0:
	# sample barcode
	demux = pd.read_csv(DEMUX_FILE,sep="\t",header=0)

	# generate all filenames based on the borcodes
	metadata = metadata.merge(demux,how='cross')
	metadata['RunRoot'] = metadata.Run.str.extract('(.*)\/demux',expand=True)
	metadata.Run = metadata.Run.str.cat(metadata.Barcode,sep='/')
	metadata.Sample = metadata.Sample.str.cat(metadata.Sample_ID,sep='/')
	metadata.Patient = metadata.Sample.str.cat(metadata.Barcode,sep='/')
	MUX = True

runs = list(set(metadata.Run.tolist()))
print(runs)
samples = list(set(metadata.Sample.tolist()))
print(samples)
kit_dict = pd.Series(metadata.Kit.values,index=metadata.Run).to_dict()
flowcell_dict = pd.Series(metadata.Flowcell.values,index=metadata.Run).to_dict()

os.chdir(os.path.join(WORKING_DIR,PROJECT_NAME))
print("current directory: ",os.getcwd())


def get_reference(wildcards):
	"""
	Get reference genome
	"""
	if wildcards.refid == HG19:
		return config[REFERENCE][HG19]
	elif wildcards.refid == HG38:
		return config[REFERENCE][HG38]
	elif wildcards.refid == T2T:
		return config[REFERENCE][T2T]
	else:
		return None


# include modules
include: "rules/sv.smk"

rule all:
	input:
		expand(["Process/{sample}/sequencing_summary.txt",
		        "Process/{sample}/QC/NanoPlot-report.html",
		        "Process/{sample}/filt.fastq",
		        "Process/{sample}/{refid}/ngmlr_{refid}.bam",
		        "Process/{sample}/{refid}/ngmlr_{refid}.bam.bai",
		        "Process/{sample}/{refid}/ngmlr_{refid}.stats.txt",
		        "Process/{sample}/{refid}/ngmlr_{refid}.sniffles.vcf",
		        "Process/{sample}/{refid}/ngmlr_{refid}.sniffles.conf.vcf",
		        "Process/{sample}/{refid}/ngmlr_{refid}.svim.vcf",
		        "Process/{sample}/{refid}/coverage_{refid}.bw",
		        "Process/{sample}/{refid}/ngmlr_{refid}.sniffles.bedpe",
		        "Process/{sample}/{refid}/ngmlr_{refid}.sniffles.conf.bedpe",
		        ],sample=samples,refid=REFS)


rule merge_fastq:
	input:
		fastq=lambda wildcards: ["{run}".format(run=row.Run) for index, row in
		                         metadata[metadata.Sample == wildcards.sample].iterrows()],
	output:
		allfastq=temp("Process/{sample}/all.fastq"),
	resources:
		tmpdir=TMP_DIR
	shell:
		"""find {input.fastq} -name '*.fastq' | xargs cat > {output.allfastq}"""

# single sample
if MUX == False:
	rule merge_summary:
		input:
			sumfiles=lambda wildcards: ["{run}/sequencing_summary.txt".format(run=row.Run) for index, row in
			                            metadata[metadata.Sample == wildcards.sample].iterrows()],
		output:
			seqsum="Process/{sample}/sequencing_summary.txt"
		resources:
			tmpdir=TMP_DIR
		params:
			shead=HEADER_SUM
		shell:
			"""
			cat {params.shead} > {output.seqsum} && \
			cat {input.sumfiles} | grep -v "passes_filtering" >> {output.seqsum}
			"""
# multiplexed samples
else:
	rule merge_summary:
		input:
			sumfiles=lambda wildcards: ["{run}/fastq/sequencing_summary.txt".format(run=row.RunRoot) for index, row in
			                            metadata[metadata.Sample == wildcards.sample].iterrows()],
			barcodesfile=lambda wildcards: ["{run}/demux/barcoding_summary.txt".format(run=row.RunRoot) for index, row
			                                in metadata[metadata.Sample == wildcards.sample].iterrows()]
		output:
			seqsum="Process/{sample}/sequencing_summary.txt"
		resources:
			tmpdir=TMP_DIR
		params:
			shead=HEADER_SUM,
			barcode=lambda wildcards:
			[row.Barcode for index, row in metadata[metadata.Sample == wildcards.sample].iterrows()][0]
		shell:
			"""
			cat {params.shead} > {output.seqsum} && \
			mkdir -p {resources.tmpdir}/{params.barcode} && \
			cat {input.barcodesfile} > '{resources.tmpdir}/{params.barcode}/input.txt' && \
			cat {input.sumfiles} > '{resources.tmpdir}/{params.barcode}/sumfiles.txt' && \
			grep '{params.barcode}' '{resources.tmpdir}/{params.barcode}/input.txt' | awk '{{print $1}}' | grep -f- '{resources.tmpdir}/{params.barcode}/sumfiles.txt' >> {output.seqsum}
			"""

rule nanoplot:
	input:
		"Process/{sample}/sequencing_summary.txt"
	output:
		report="Process/{sample}/QC/NanoPlot-report.html"
	resources:
		tmpdir=TMP_DIR
	conda:
		"envs/qc-env.yaml"
	params:
		sample="{sample}",
		output_dir="Process/{sample}/QC",
	shell:
		"NanoPlot --summary {input} --outdir {params.output_dir} --N50 --title {params.sample}"

rule nanfilt:
	input:
		fastqin="Process/{sample}/all.fastq"
	output:
		fastqout=temp("Process/{sample}/filt.fastq")
	resources:
		tmpdir=TMP_DIR
	conda:
		"envs/mapping-env.yaml"
	params:
		rlen=config['read_length']
	shell:
		"NanoFilt -l {params.rlen} --headcrop 50 --tailcrop 50 --readtype 1D {input.fastqin} > {output.fastqout}"

rule ngmlr_mock:
	input:
		fastq="Process/{sample}/filt.fastq"
	output:
		# outf = directory(expand("Process/{sample}/{refid}/", sample=["{sample}"], refid=[HG19, HG38])),
		outtemp=expand("Process/{sample}/{refid}/temp",sample=["{sample}"],refid=REFS)
	resources:
		tmpdir=TMP_DIR
	conda:
		"envs/mapping-env.yaml"
	shell:
		"touch {output.outtemp}"

rule ngmlr:
	input:
		fastq="Process/{sample}/filt.fastq",
		reference=get_reference,
		ftemp="Process/{sample}/{refid}/temp"
	output:
		bam="Process/{sample}/{refid}/ngmlr_{refid}.bam",
		bai="Process/{sample}/{refid}/ngmlr_{refid}.bam.bai",
		sam=temp("Process/{sample}/{refid}/ngmlr_{refid}.sam")
	resources:
		tmpdir=TMP_DIR
	conda:
		"envs/mapping-env.yaml"
	params:
		threads=THREADS
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
		"Process/{sample}/{refid}/coverage_{refid}.bw"
	resources:
		tmpdir=TMP_DIR
	conda:
		"envs/mapping-env.yaml"
	params:
		threads=16
	shell:
		"bamCoverage --bam {input.bam} -o {output} --binSize 50 -p {params.threads}"

rule stats:
	input:
		bam="Process/{sample}/{refid}/ngmlr_{refid}.bam",
		bai="Process/{sample}/{refid}/ngmlr_{refid}.bam.bai"
	output:
		"Process/{sample}/{refid}/ngmlr_{refid}.stats.txt"
	resources:
		tmpdir=TMP_DIR
	conda:
		"envs/mapping-env.yaml"
	shell:
		"samtools stats {input.bam} > {output}"

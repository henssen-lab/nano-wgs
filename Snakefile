import os

import pandas as pd
from snakemake.utils import min_version

# 7.32.4
#min_version("6.4.1")

configfile: "configs/config_run.yaml"

HG19 = 'hg19'
HG38 = 'hg38'
T2T = 'T2T-CHM13'
MM39 = 'mm39'
MM39C = 'mm39+construct'

genomes = [HG19,HG38,T2T,MM39,MM39C]

REFERENCE = 'reference'
THREADS = config["threads"]
MUX = False

REFS = config['refs'].split(",")  # hg19,hg38
PROJECT_NAME = config['project_name']
WORKING_DIR = config['working_dir']
METADATA_FILE = config['metafile']
DEMUX_FILE = config['muxfile']
TMP_DIR = config['tmp_dir']
NAME_FASTQ_FOLDER = config['name_fastq_folder']
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
samples = list(set(metadata.Sample.tolist()))
kit_dict = pd.Series(metadata.Kit.values,index=metadata.Run).to_dict()
flowcell_dict = pd.Series(metadata.Flowcell.values,index=metadata.Run).to_dict()

#os.chdir(os.path.join(WORKING_DIR,PROJECT_NAME))
#print("current directory: ",os.getcwd())


def get_reference(wildcards):
	"""
	Get reference genome
	"""
	if wildcards.refid in genomes:
		return config[REFERENCE][wildcards.refid]
	else:
		return None

# include modules
include: "rules/sv.smk"

rule all:
	input:
		expand(["{workdir}/{project}/Process/{sample}/sequencing_summary.txt",
		        "{workdir}/{project}/Process/{sample}/QC/NanoPlot-report.html",
		        "{workdir}/{project}/Process/{sample}/filt.fastq",
		        "{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam",
		        "{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam.bai",
		        "{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.stats.txt",
		        "{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sniffles.vcf",
		        "{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sniffles.conf.vcf",
		        "{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.svim.vcf",
		        "{workdir}/{project}/Process/{sample}/{refid}/coverage_{refid}.bw",
		        "{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sniffles.bedpe",
		        "{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sniffles.conf.bedpe",
		        ],sample=samples,refid=REFS,workdir=WORKING_DIR,project=PROJECT_NAME)


rule merge_fastq:
	input:
		fastq=lambda wildcards: ["{workdir}/{project}/{run}/{subdir}".format(workdir=wildcards.workdir,project=wildcards.project,run=row.Run,subdir=NAME_FASTQ_FOLDER)
					for index, row in metadata[metadata.Sample == wildcards.sample].iterrows()],
	output:
		allfastq=temp("{workdir}/{project}/Process/{sample}/all.fastq"),
	resources:
		tmpdir=TMP_DIR
	shell:
		"""
		find {input.fastq} -name '*.fastq' -print0 | xargs -0 -r cat > {output.allfastq} && \
		find {input.fastq} -name '*.fastq.gz' -print0 | xargs -0 -r zcat >> {output.allfastq}		   
		"""

# single sample
if MUX == False:
	rule merge_summary:
		input:
			sumfiles=lambda wildcards: ["{workdir}/{project}/{run}/sequencing_summary.txt".format(workdir=wildcards.workdir,project=wildcards.project,run=row.Run) for index, row in
			                            metadata[metadata.Sample == wildcards.sample].iterrows()],
		output:
			seqsum="{workdir}/{project}/Process/{sample}/sequencing_summary.txt"
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
			sumfiles=lambda wildcards: ["{workdir}/{project}/{run}/sequencing_summary.txt".format(workdir=wildcards.workdir,project=wildcards.project,run=row.Run) for index, row in
			                            metadata[metadata.Sample == wildcards.sample].iterrows()],
			barcodesfile=lambda wildcards: ["{workdir}/{project}/{run}/barcoding_summary.txt".format(run=row.RunRoot) for index, row
			                                in metadata[metadata.Sample == wildcards.sample].iterrows()]
		output:
			seqsum="{workdir}/{project}/Process/{sample}/sequencing_summary.txt"
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
		"{workdir}/{project}/Process/{sample}/sequencing_summary.txt"
	output:
		report="{workdir}/{project}/Process/{sample}/QC/NanoPlot-report.html"
	resources:
		tmpdir=TMP_DIR
	conda:
		"envs/qc-env.yaml"
	params:
		sample="{sample}",
		output_dir="{workdir}/{project}/Process/{sample}/QC",
	shell:
		"NanoPlot --summary {input} --outdir {params.output_dir} --N50 --title {params.sample}"

rule nanfilt:
	input:
		fastqin="{workdir}/{project}/Process/{sample}/all.fastq"
	output:
		fastqout=temp("{workdir}/{project}/Process/{sample}/filt.fastq")
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
		fastq="{workdir}/{project}/Process/{sample}/filt.fastq"
	output:
		# outf = directory(expand("{workdir}/{project}/Process/{sample}/{refid}/", sample=["{sample}"], refid=[HG19, HG38])),
		outtemp=expand("{workdir}/{project}/Process/{sample}/{refid}/temp",sample=["{sample}"],workdir=["{workdir}"],project=["{project}"],refid=REFS)
	resources:
		tmpdir=TMP_DIR
	conda:
		"envs/mapping-env.yaml"
	shell:
		"touch {output.outtemp}"

rule ngmlr:
	input:
		fastq="{workdir}/{project}/Process/{sample}/filt.fastq",
		reference=get_reference,
		ftemp="{workdir}/{project}/Process/{sample}/{refid}/temp"
	output:
		sam="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sam"
	resources:
		tmpdir=TMP_DIR
	conda:
		"envs/mapping-env.yaml"
	params:
		threads=THREADS
	shell:
		"""
		ngmlr --bam-fix --threads {params.threads} --reference {input.reference} --query {input.fastq} --output {output.sam} --presets ont
		"""

rule ngmlr_sort:
	input:
		sam="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sam"
	output:
		bam="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam",
                bai="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam.bai",
	params:
		threads=8
	shell:
		"""
		samtools sort -@ {params.threads} -O BAM -o {output.bam} {input.sam}
                samtools index {output.bam}
                """

rule coverage:
	input:
		bam="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam",
		bai="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam.bai"
	output:
		"{workdir}/{project}/Process/{sample}/{refid}/coverage_{refid}.bw"
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
		bam="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam",
		bai="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam.bai"
	output:
		"{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.stats.txt"
	resources:
		tmpdir=TMP_DIR
	conda:
		"envs/mapping-env.yaml"
	shell:
		"samtools stats {input.bam} > {output}"

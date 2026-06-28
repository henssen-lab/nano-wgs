"""
SV calling.
"""
rule sniffles:
    input:
        bam="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam",
        bai="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam.bai"
    output:
        sniff="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sniffles.vcf",
	sniffconf="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sniffles.conf.vcf"
    resources:
        tmpdir = TMP_DIR
    conda:
        "../envs/snv-env.yaml"
    params:
        threads = 8
    shell:
         """
	sniffles -t {params.threads} -m {input.bam} -v {output.sniff} --min_homo_af 0.7 --min_het_af 0.1 --min_length 50 --cluster --genotype --min_support 4 --report-seq
	sniffles -t {params.threads} -m {input.bam} -v {output.sniffconf} --min_homo_af 0.7 --min_het_af 0.1 --min_length 500 --cluster --genotype --min_support 10 --report-seq
	"""

rule svim:
    input:
        bam="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam",
        bai="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.bam.bai",
        reference=get_reference
    output:
        "{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.svim.vcf"
    resources:
        tmpdir = TMP_DIR
    conda:
        "../envs/snv-env.yaml"
    params:
        svim_output_dir = "{workdir}/{project}/Process/{sample}/{refid}/svim/",
        sample = "{sample}"
    shell:
         "svim alignment --read_names --insertion_sequences --sample {params.sample} {params.svim_output_dir} {input.bam} {input.reference} && mv {params.svim_output_dir}variants.vcf {output}"

rule survivor:
    input:
        unfilt="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sniffles.vcf",
	filt="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sniffles.conf.vcf"
    output:
        unfilt="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sniffles.bedpe",
	filt="{workdir}/{project}/Process/{sample}/{refid}/ngmlr_{refid}.sniffles.conf.bedpe"
    conda:
        "../envs/sv.yaml"
    log:
        "{workdir}/{project}/Process/{sample}/{refid}/logs/survivor.log"
    shell:
        """
        SURVIVOR vcftobed {input.unfilt} -1 -1 {output.unfilt} &> {log}
	SURVIVOR vcftobed {input.filt} -1 -1 {output.filt} &>> {log}
        """

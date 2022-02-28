"""
SV calling.
"""
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
        threads = 16
    shell:
         "sniffles -t {params.threads} -m {input.bam} -v {output} --min_homo_af 0.7 --min_het_af 0.1 --min_length 50 --cluster --genotype --min_support 4 --report-seq"

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

rule survivor:
    input:
        unfilt="Process/{sample}/{refid}/ngmlr_{refid}.sniffles.vcf"
    output:
        unfilt="Process/{sample}/{refid}/ngmlr_{refid}.sniffles.bedpe"
    conda:
        "../envs/sv.yaml"
    log:
        "Process/{sample}/logs/survivor.log"
    shell:
        """
        SURVIVOR vcftobed {input.unfilt} -1 -1 {output.unfilt} &> {log}
        """

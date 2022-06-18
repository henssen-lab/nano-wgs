#!/bin/bash

fin=$1
region=$2

samtools depth -a -b $regions ${root}/${sample}/hg38/ngmlr_hg38.bam > ${sample}.cov
cat ${sample}.cov | awk '{print $1"\t"$2"\t"$2+2"\t"$3}' > ${sample}.bed
bedtools intersect -wa -wb -a $regions -b ${sample}.bed | awk '{print $5"\t"$6"\t"$7"\t"$4"\t"$8}' | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,5 -o distinct,mean | awk -v a="$sample" '{print $0"\t"a}' > ${sample}_merge.bed


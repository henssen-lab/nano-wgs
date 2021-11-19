#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH -p long
#SBATCH --time 7-00:00:00

set -x
set -e

#datadir=/fast/users/giurgium_c/CircleSeq/Nanopore/Evolution/Process/CB1008-initial_01012020/hg19
datadir=$1
fin=${datadir}/ngmlr_hg19.bam
#sample=CB1008-initial
sample=$2

datatemp="/fast/users/giurgium_c/scratch/CNV/${sample}"
mkdir -p $datatemp

scripts_root=/fast/users/giurgium_c/work/nano-wgs/cnv/smurfseq_scripts/scripts/
scripts_data=/fast/users/giurgium_c/work/nano-wgs/cnv/smurfseq_scripts/data
scripts_cnv=/fast/users/giurgium_c/work/nano-wgs/cnv/smurfseq_scripts/scripts/cnvAnalysis.R

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

echo "$fin"
echo "check bam and bed for consistency in chromosome annotation!!"

mkdir -p ${datadir}/cnv_50k 
cd ${scripts_root}
./filterAlnScoreAndQual.py -i $fin  -o $datatemp/unambig_smurf_frags.sam -s 120 -q 1

samtools view -H $fin  |\
   sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' | \
   sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' | \
   sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' | \
   sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' | \
   sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | \
   sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' | \
   sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' | \
   sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' | \
   sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' | \
   sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' | \
   sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' | \
   sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' | \
   sed -e 's/SN:MT/SN:chrM/' > h.sam

samtools view -@16 -h -O BAM $datatemp/unambig_smurf_frags.sam > $datatemp/unambig_smurf_frags.bam
samtools reheader h.sam $datatemp/unambig_smurf_frags.bam | samtools view -h -O SAM  - > $datatemp/unambig_smurf_frags_chr.sam
#samtools reheader h.sam $fin | samtools view -h -O SAM  - > $datatemp/unambig_smurf_frags_chr.sam

./getBinCounts.py -v -i $datatemp/unambig_smurf_frags_chr.sam  -c ../data/hg19.chrom.sizes -b ../data/bins_50k_hg19.bed -o $datatemp/bin_counts.bed -s ../data/bin_stats.txt

cd $datadir/cnv_50k/
Rscript $scripts_cnv $datatemp/bin_counts.bed $sample ${scripts_data}/bins_50k_hg19_gc.txt ${scripts_data}/bins_50k_hg19_exclude.txt

rm $datatemp/unambig_smurf_frags*

import os
import sys
import pysam as ps
import numpy as np
import matplotlib
if sys.platform == "linux" or sys.platform == "linux2":
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set(style="whitegrid")

sam_fname = sys.argv[1]
paf_fname = sys.argv[2]
savepath = sam_fname + "-samanalysis/"
os.system("mkdir " + savepath)

text_output_fname = savepath + "stats.txt"
orig_stdout = sys.stdout
text_output_handle = open(text_output_fname, 'w')
sys.stdout = text_output_handle

samfile = ps.AlignmentFile(sam_fname, "r")
ref_lengths_dict = {}

read_lengths = [] # How long are the reads?
align_lengths = [] # How many bases align?
mapping_qualities = []
is_unmapped = []
is_secondary = []
is_supplementary = []
query_name = []

for read in samfile.fetch():
        query_name.append(read.query_name)
        is_unmapped.append(read.is_unmapped)
        is_secondary.append(read.is_secondary)
        is_supplementary.append(read.is_supplementary)
        mapping_qualities.append(read.mapping_quality)
        read_lengths.append(read.infer_read_length())
        if read.is_unmapped:
            align_lengths.append(float('nan'))
            #cigars.append("NA")
        else:
            align_lengths.append(read.get_overlap(read.reference_start, read.reference_end))
            #cigars.append(read.cigarstring)
sam_df = pd.DataFrame({'ReadName':query_name, 'isUnmapped':is_unmapped, 'isSecondary':is_secondary, 'isSupplementary':is_supplementary, 'ReadLength':read_lengths, 'AlignLength':align_lengths, 'MappingQuality':mapping_qualities})

paf_df = pd.read_csv(paf_fname, sep="\t", header=None, usecols=list(range(0,12)), names=["QuerySequenceName", "QuerySequenceLength", "QueryStart", "QueryEnd", "RelativeStrand", "TargetSequenceName", "TargetSequenceLength", "TargetStartOnOriginalStrand", "TargetEndOnOriginalStrand", "NumberOfResidueMatches", "AlignmentBlockLength", "MappingQuality"])
paf_df["AlignmentIdentity"] = paf_df["NumberOfResidueMatches"] / paf_df["AlignmentBlockLength"]

n_reads = sam_df[sam_df["isUnmapped"]].shape[0] + sam_df[~sam_df.isUnmapped & ~sam_df.isSecondary  & ~sam_df.isSupplementary].shape[0]
print("Number of Reads: ", n_reads, "\n")
print("Number of entries in sam (Mappings + Unmapped Reads): ", sam_df.shape[0], "\n")
print("Number of mappings: ", sam_df[~sam_df.isUnmapped].shape[0], "\n")
print("Number of unmapped reads: ", sam_df[sam_df["isUnmapped"]].shape[0], "\n")
print("Percent of unmapped reads: ", (100*(sam_df[sam_df["isUnmapped"]].shape[0] / n_reads)), "%\n")
print("Number of primary mappings: ", sam_df[~sam_df.isUnmapped & ~sam_df.isSecondary  & ~sam_df.isSupplementary].shape[0], "\n")
#print("Number of unique mappings: ", unique_mappings.shape[0])

## READ LENGTH DISTRIBUTION
print("Read Length Mean:\t", sam_df.ReadLength.mean(skipna=True), "\n")
print("Read Length STD:\t", sam_df.ReadLength.std(skipna=True), "\n")
v = np.array(sam_df.ReadLength)
v = v[~np.isnan(v)]
sns.distplot(v, kde=True).get_figure().savefig(savepath + "ReadLengthDist.png")
plt.clf()

## MAPPING QUALITY DIST
print("Mapping Quality Mean:\t", sam_df.MappingQuality.mean(), "\n")
print("Mapping Quality STD:\t", sam_df.MappingQuality.std(), "\n")
v = np.array(sam_df.MappingQuality)
v = v[~np.isnan(v)]
sns.distplot(v, kde=False).get_figure().savefig(savepath + "MappingQualityDist.png")
plt.clf()

## READ LENGTH VS. MAPPING Quality
sns.jointplot(x='ReadLength', y='MappingQuality', data=sam_df)
plt.savefig(savepath + "ReadLengthVSMappingQuality.png")
plt.clf()

## READ LENGTHS VS. ALIGN lengths
df = pd.DataFrame({'read_lengths':read_lengths, 'align_lengths':align_lengths})
sns.jointplot(x='align_lengths', y='read_lengths', data=df, alpha=0.1)
plt.savefig(savepath + "ReadLengthVSAlignLength.png")
plt.clf()

## ALIGNMENT IDENTITY DIST
print("Alignment Identity Mean:\t", np.mean(paf_df["AlignmentIdentity"]), "\n")
print("Alignment Identity STD:\t", np.std(paf_df["AlignmentIdentity"]), "\n")
sns.distplot(paf_df["AlignmentIdentity"], kde=True).get_figure().savefig(savepath + "AlignmentIdentityDist.png")
plt.clf()

## AlignLengths
print("Alignment Length Mean:\t", df.align_lengths.mean(skipna=True), "\n")
print("Alignment Length STD:\t", df.align_lengths.std(skipna=True), "\n")
v = np.array(df.align_lengths)
v = v[~np.isnan(v)]
sns.distplot(v).get_figure().savefig(savepath + "AlignmentLengths.png")
plt.clf()

## WRAP
sys.stdout = orig_stdout
text_output_handle.close()

#!/usr/bin/env python
# coding: utf-8

# Packages
import sys
import argparse
import glob
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

## Data csv loaded in from command-line
parser = argparse.ArgumentParser()
parser.add_argument("csv_table")
parser.add_argument("exp_fragments")
args = parser.parse_args()

## Load data csv into pandas dataframe
hg38_meta = pd.read_csv(args.csv_table, sep=',', header=0)

## Make new perctenage alignment columns
hg38_meta['hg38_alignment_rate'] = hg38_meta.loc[:, ('bt2_total_aligned')] / hg38_meta.loc[:, ('total_reads')] *100
hg38_meta['ecoli_alignment_rate'] = hg38_meta.loc[:, ('bt2_spike_total_aligned')] / hg38_meta.loc[:, ('total_reads')] *100
# print(hg38_meta)

## Table containing just the data used for this plot; columns total_reads, bt2_total_aligned, hg38_alignment_rate, ecoli_alignment_rate
seq_summary_table = hg38_meta.loc[:, ('sample_id', 'experiment', 'total_reads', 'bt2_total_aligned', 'hg38_alignment_rate', 'ecoli_alignment_rate')]
seq_summary_table.to_csv('./hg38_seq_summary_table.csv', index=False)
# print(seq_summary_table)

# Apply the default sns theme
# sns.set_theme()

# ---------- Plot 1 - Alignment Summary --------- #
## Construct quad plot
fig, seq_summary = plt.subplots(2,2)
fig.suptitle("Sequencing and Alignment Summary")

# Seq depth
sns.boxplot(data=hg38_meta, x='experiment', y='total_reads', order=['h3k27me3', 'h3k4me3', 'igg'], ax=seq_summary[0,0])
seq_summary[0,0].set_title("Sequencing Depth")
seq_summary[0,0].set_ylabel("Total Reads")

# Alignable fragments
sns.boxplot(data=hg38_meta, x='experiment', y='bt2_total_aligned', order=['h3k27me3', 'h3k4me3', 'igg'], ax=seq_summary[0,1])
seq_summary[0,1].set_title("Alignable Fragments")
seq_summary[0,1].set_ylabel("Total Aligned Reads")

# Alignment rate hg38
sns.boxplot(data=hg38_meta, x='experiment', y='hg38_alignment_rate', order=['h3k27me3', 'h3k4me3', 'igg'], ax=seq_summary[1,0])
seq_summary[1,0].set_title("Alignment Rate (hg38)")
seq_summary[1,0].set_ylabel("Percent of Fragments Aligned")

# Alignment rate e.coli
sns.boxplot(data=hg38_meta, x='experiment', y='ecoli_alignment_rate', order=['h3k27me3', 'h3k4me3', 'igg'], ax=seq_summary[1,1])
seq_summary[1,1].set_title("Alignment Rate (e.coli)")
seq_summary[1,1].set_ylabel("Percent of Fragments Aligned")

plt.subplots_adjust(wspace=0.5, hspace=0.5)
# plt.show()
fig.savefig('hg38_seq_summary_seaborn.png')

# ---------- Plots 2&3 - Mapped Fragment Distribution --------- #

# Parse './' to fragments argument as we expect all data to be in the same directory. 
fragments = args.exp_fragments

# Create list of deeptools raw fragment files
dt_frag_list = glob.glob(fragments + '*raw.csv')

df_list = list()
for i in list(range(0,len(dt_frag_list))):
    # create dataframe from csv file for each file and save to a list
    dt_frag_i = pd.read_csv(dt_frag_list[i], sep='\t', skiprows=[0], header=0)
    df_list.append( dt_frag_i )

    # create long forms of fragment histograms
    dt_frag_i_long = np.repeat(dt_frag_i['Size'].values, dt_frag_i['Occurrences'].values)
    dt_sample_i_long = np.repeat(dt_frag_i['Sample'][0], len(dt_frag_i_long))

    if i==0:
        frags_arr = dt_frag_i_long
        sample_arr = dt_sample_i_long
        og_frag_df = dt_frag_i
    else:
        frags_arr = np.append(frags_arr, dt_frag_i_long)
        sample_arr = np.append(sample_arr, dt_sample_i_long)
        og_frag_df = og_frag_df.append(dt_frag_i)

# create hue array using regex pattern matching
for i in list(range(0,len(sample_arr))):
    sample_i = sample_arr[i]
    sample_exp = re.findall("^[^_]*", sample_i)

    if i==0:
        sample_exp_arr = np.array(sample_exp[0])
    else:
        sample_exp_arr = np.append(sample_exp_arr, sample_exp[0])

df_long = pd.DataFrame( { "fragment_size" : frags_arr, "sample" : sample_arr , "histone_mark": sample_exp_arr}, index = np.arange(len(frags_arr)))
df_long.to_csv('./fragmanet_distribution_violin.csv', index=False)
og_frag_df.to_csv('./fragmanet_distribution_line.csv', index=False)
# print(df_long)

plt.clf()
ax1 = sns.violinplot(data=df_long, x="sample", y="fragment_size", hue="histone_mark")
# plt.show()
fig1 = ax1.get_figure()
fig1.savefig('fragmanet_distribution_violin.png')

plt.clf()
ax2 = sns.lineplot(data=og_frag_df, x="Size", y="Occurrences", hue="Sample")
# plt.show()
fig2 = ax2.get_figure()
fig2.savefig('fragmanet_distribution_line.png')

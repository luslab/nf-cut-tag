#!/usr/bin/env python
# coding: utf-8

# Packages
import sys
import argparse
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## Data csv loaded in from command-line
parser = argparse.ArgumentParser()
parser.add_argument("csv_table")
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

############
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
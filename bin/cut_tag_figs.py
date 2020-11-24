#!/usr/bin/env python
# coding: utf-8

# Packages
import numpy as np
import pandas as pd

# Load data csv
# chr20_meta = pd.read_csv("/Users/westc/dev/python_scripts/cut-tag-charting/meta_table.csv", sep=',', header=0)
chr20_meta = pd.read_csv("./meta_table.csv", sep=',', header=0)
print(chr20_meta)

# Subset dataframe and save to file
meta_seq_depth = chr20_meta[['sample_id', 'experiment', 'total_reads']]
print(meta_seq_depth)

# Save dataframe to file
# meta_seq_depth.to_csv('/Users/westc/dev/python_scripts/cut-tag-charting/meta_seq_depth.csv', header=True, index=False)
meta_seq_depth.to_csv('./meta_seq_depth.csv', header=True, index=False)
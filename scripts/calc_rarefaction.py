"""Find abundance of hits in blast search
"""

import pandas as pd
import os
import sys
import argparse
import re
import random

parser = argparse.ArgumentParser(description=
        'Calcualte rarefaction points')
parser.add_argument('--pc_table',
        help=('Tab seperate OTU style table with PC abundances'),
        type=str,
        required=True)
parser.add_argument('--iterations',
        help=('how many times to sample the data randomly'),
        type=int,
        required=True)

args = parser.parse_args()

pc_df = pd.read_table(args.pc_table, header=0, index_col=0, na_values='0')
output_handle = open('pc_rarefaction.txt', 'w')

n = 0
while n < args.iterations:
    n = n + 1
    set1 = set()
    sample_count = 0
    list_samples = pc_df.columns.tolist()
    random.shuffle(list_samples)
    while list_samples:
        sample = list_samples.pop()
        sample_count += 1
        sample_series = pc_df[sample].dropna()
        pcs = sample_series.index.values
        set1.update(pcs)
        output = '%s\t%s\t%s\n' % (n, sample_count, len(set1))
        output_handle.write(output)

output_handle.close()
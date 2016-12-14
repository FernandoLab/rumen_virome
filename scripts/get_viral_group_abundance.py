#! /usr/bin/env python

"""
"""

import os
import sys
import argparse
import re

group_list = []
with open('blast2network.txt-Partition-Contig-Sizes-filter.tab', 'rU') as group_file:
    lines = group_file.read().splitlines()
    for line in lines:
        group_list.append(line.split()[0])

contig_group_dict = {}
with open('blast2network.txt-Partition-Membership.tab', 'rU') as mem_file:
    lines = mem_file.read().splitlines()
    for line in lines:
        contig_group_dict[line.split()[0]] = line.split()[1]

files = []
for filename in os.listdir('vmg_ion_illumina_contigs_bowtie'):
    if filename.split('.', 1)[1] == 'contig_raw_abundance.txt':
        filename = os.path.join('vmg_ion_illumina_contigs_bowtie', filename)
        files.append(filename)

for file in files:
    with open(file, 'rU') as ab_file:
        group_ab_dict = {}
        lines = ab_file.read().splitlines()
        for line in lines:
            if (line.startswith('NODE_') and line.split()[0] in contig_group_dict):
                group = contig_group_dict[line.split()[0]]
                if group in group_list:
                    if group in group_ab_dict:
                        group_ab_dict[group] += int(line.split()[1])
                    else:
                        group_ab_dict[group] = int(line.split()[1])
        output_file = file.split('.', 1)[0] + '.group_raw_abundance.txt'
        output_handle = open(output_file, 'w')
        sample = file.split('/')[1].split('.', 1)[0]
        output_handle.write('OTU' + '\t' + sample + '\n')
        for group, ab in group_ab_dict.iteritems():
            if ab > 0:
                output = '%s\t%d\n' % (group, ab)
                output_handle.write(output)
        output_handle.close()


#! /usr/bin/env python

"""

"""

import argparse
import sys
import re

parser = argparse.ArgumentParser(description='filter viral KOs so a different ORF from the same contig has a viral hit')
parser.add_argument('-k', '--kegg_db_blast', help='blast file from ORFs vs. KEGG database', type=str, required=True)
parser.add_argument('-v', '--viral_db_blast', help='blast file from ORFs vs. viral database', type=str, required=True)
args = parser.parse_args()

viral_contig_dict = {}
viral_orf_dict = {}
in_viral = open(args.viral_db_blast, 'rU')
for line in in_viral:
    stripline = line.rstrip()
    blast_tabs = stripline.split('\t')
    match = re.search(r'NODE\S+\_cov', stripline)
    viral_contig_dict[match.group()] = blast_tabs[0]
    viral_orf_dict[blast_tabs[0]] = line,

in_viral.close()

contig_dict = {}
in_kegg = open(args.kegg_db_blast, 'rU')
for line in in_kegg:
    stripline = line.rstrip()
    blast_tabs = stripline.split('\t')
    match = re.search(r'NODE\S+\_cov', stripline)
    if match.group() in viral_contig_dict:
        if blast_tabs[0] not in viral_orf_dict:
            contig_dict[match.group()] = blast_tabs[0]

in_kegg.close()
out_kegg = open('viral_filtered_kegg_blast.txt', 'wb')

in_kegg = open(args.kegg_db_blast, 'rU')
for line in in_kegg:
    stripline = line.rstrip()
    blast_tabs = stripline.split('\t')
    match = re.search(r'NODE\S+\_cov', stripline)
    if match.group() in contig_dict:
        out_kegg.write(line,)

in_kegg.close()
out_kegg.close()


#! /usr/bin/env python

"""
"""

import os
import sys
import argparse
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(description =
        ('Extract contigs and ORFs from fasta files idnetified by VirSorter'))
parser.add_argument('--virsorter_csv',
        help=('viroster csv output'),
        type=str,
        required=True)
parser.add_argument('--orf_fasta',
        help=('ORF fasta file to extract ORFs from'),
        type=str,
        required=True)
parser.add_argument('--contig_fasta',
        help=('Contig fastsa file to extract contigs from'),
        type=str,
        required=True)
args = parser.parse_args()

vir_dict = {}
vir_count = 0
with open(args.virsorter_csv, 'rU') as vir_csv:
    for row in vir_csv:
        if row.strip().startswith('## 4 - Prophages'):
            break

        if not row.strip().startswith('#'):
            split_row = row.strip().split(',')
            vir_count += 1
            contig = split_row[0].split('_', 1)[1]
            if '-' in contig:
                contig = contig.rsplit('-', 1)[0]
            vir_dict[contig] = split_row[4]

contig_found = 0
contig_rec = []
for record in SeqIO.parse(open(args.contig_fasta, 'rU'), 'fasta'):
    id = record.id.replace('.', '_')
    if id in vir_dict:
        contig_found += 1
        contig_rec.append(record)

orf_found = 0
orf_rec = []
orf_contig_dict = {}
for record in SeqIO.parse(open(args.orf_fasta, 'rU'), 'fasta'):
    contig = record.id.rsplit('_', 1)[0].replace('.', '_')
    if contig in vir_dict:
        orf_contig_dict[contig] = record.id
        orf_found += 1
        orf_rec.append(record)

output_handle_contig = open('virsorter_contig.fasta', 'w')
contig_out = FastaIO.FastaWriter(output_handle_contig, wrap=None)
contig_out.write_file(contig_rec)

output_handle_orf = open('virsorter_orf.fasta', 'w')
orf_out = FastaIO.FastaWriter(output_handle_orf, wrap=None)
orf_out.write_file(orf_rec)

output_handle_contig.close()
output_handle_orf.close()
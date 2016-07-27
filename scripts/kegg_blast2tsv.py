#! /usr/bin/env python

"""Calculate KEGG KO abundance from BLAST and Bowtie inputs.
Use '-h' for parameter help.
"""

import os
import re
import argparse
import blast

parser = argparse.ArgumentParser(description = 
        ('Calculate KO abundance for each sample from KEGG database'
        'BLAST result and Bowtie mapping information for ORFs.'))
parser.add_argument('--blast', 
        help=('BLAST output file from search against KEGG db. Output '
        'format should be "-m 8".'),
        type=str,
        required=True)
parser.add_argument('--bowtie_dir', 
        help=('Directory containing output from bowtie in which reads '
        'were mapped to ORFs. Filenames should end in ".psl".'),
        type=str,
        required=True)
parser.add_argument('--gene_ko_mapping', 
        help=('File containing KEGG gene to KO mapping information. '
       'Usually found at kegg/genes/ko/ko_genes.list.'),
        type=str,
        required=True)
args = parser.parse_args()


orf_gene_dict = blast.get_orf_gene_dict(args.blast)
gene_ko_dict = blast.get_gene_ko_dict(args.gene_ko_mapping)
bowtie_files = [os.path.join(args.bowtie_dir, f) 
        for f in os.listdir(args.bowtie_dir) if f.endswith('.psl')]
for bowtie_file in bowtie_files:
    match = re.search(r'([\w.]+/)([\w.]+)([.]psl)', bowtie_file)
    sample = match.group(2)

    orf_raw_ab_file = sample + '.orf_raw_abundance.txt'
    orf_norm_ab_file = sample + '.orf_norm_abundance.txt'
    ko_raw_ab_file = sample + '.ko_raw_abundance.txt'
    ko_norm_ab_file = sample + '.ko_norm_abundance.txt'

    read_count, mapped_reads = blast.count_sample_reads_bowtie(bowtie_file)
    orf_abundance_dict = blast.get_orf_abundance_dict(bowtie_file)
    ko_abundance_dict = blast.get_ko_abundance_dict(orf_gene_dict,
            gene_ko_dict, orf_abundance_dict)

    
    orf_raw_ab_output = open(orf_raw_ab_file, 'w')
    orf_raw_ab_output.write('ID\t%s\n' % sample) 
    for orf, ab in orf_abundance_dict.items():
        orf_raw_ab_output.write('%s\t%d\n' % (orf, ab))
    orf_raw_ab_output.close()
    
    orf_norm_ab_output = open(orf_norm_ab_file, 'w')
    orf_norm_ab_output.write('ID\t%s\n' % sample)     
    for orf, ab in orf_abundance_dict.items():
        norm_ab = float(ab) / float(read_count)
        orf_norm_ab_output.write('%s\t%e\n' % (orf, norm_ab))
    orf_norm_ab_output.close()

    ko_count = 0
    ko_raw_ab_output = open(ko_raw_ab_file, 'w')
    ko_raw_ab_output.write('ID\t%s\n' % sample)     
    for ko, ab in ko_abundance_dict.items():
        ko_count += 1
        ko_raw_ab_output.write('%s\t%d\n' % (ko, ab))
    ko_raw_ab_output.close()
    
    ko_norm_ab_output = open(ko_norm_ab_file, 'w')
    ko_norm_ab_output.write('ID\t%s\n' % sample)         
    for ko, ab in ko_abundance_dict.items():
        norm_ab = float(ab) / float(read_count)
        ko_norm_ab_output.write('%s\t%e\n' % (ko, norm_ab))
    ko_norm_ab_output.close()

    output = ('Sample: %s\nReads - %d\nKOs - %d\n\n' 
            % (sample, read_count, ko_count))
    print(output)

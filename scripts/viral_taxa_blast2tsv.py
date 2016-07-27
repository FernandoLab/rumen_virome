#! /usr/bin/env python

"""Calculate viral taxonomy abundance from BLAST and Bowtie inputs.
Use '-h' for parameter help.
"""

import os
import re
import argparse
import blast

parser = argparse.ArgumentParser(description = 
        ('Calculate viral taxonomy abundance for each sample from PHAST'
        'viral database BLAST results and bowtie mapping information.'))
parser.add_argument('--blast', 
        help=('BLAST output file from search against the viral db.'
        'Output format should be "-m 8".'),
        type=str,
        required=True)
parser.add_argument('--bowtie_dir', 
        help=('Directory containing output from bowtie in which reads '
        'were mapped to ORFs. Filenames should end in ".psl".'),
        type=str,
        required=True)
parser.add_argument('--viral_db',
        help=('Phast viral database used for assigning taxonomy '
        'downloaded via'
        'wget http://www.phantome.org/Downloads/Viruses/2016_05_01.tgz'),
        type=str,
        required=True)
args = parser.parse_args()

orf_gene_dict = blast.get_orf_gene_dict(args.blast)
blast.setup_viral_taxa_sql_db()
viral_orf_fam_dict, viral_orf_genus_dict = blast.get_viral_orf_taxa_dicts(
        orf_gene_dict, args.viral_db)
bowtie_files = [os.path.join(args.bowtie_dir, f) 
        for f in os.listdir(args.bowtie_dir) if f.endswith('.psl')]
count_mapped_reads = 0
count_fam_reads = 0
count_genus_reads = 0
count_total_reads = 0
for bowtie_file in bowtie_files:
    match = re.search(r'([\w.]+/)([\w.]+)([.]psl)', bowtie_file)
    sample = match.group(2)
    read_count, mapped_reads  = blast.count_sample_reads_bowtie(bowtie_file)
    count_mapped_reads = count_mapped_reads + mapped_reads
    count_total_reads = count_total_reads + read_count

    orf_abundance_dict = blast.get_orf_abundance_dict(bowtie_file)
    
    sample_viral_fam_dict, sample_viral_genus_dict = (
            blast.get_sample_viral_taxa_abundances(
            orf_abundance_dict,viral_orf_fam_dict,viral_orf_genus_dict)
            )

    fam_count = 0
    viral_raw_fam_ab_file = sample + '.fam_raw_abundance.txt'
    viral_raw_fam_ab_output = open(viral_raw_fam_ab_file, 'w')
    viral_raw_fam_ab_output.write('ID\t%s\n' % sample)     
    for fam, ab in sample_viral_fam_dict.items():
        fam_count += 1
        viral_raw_fam_ab_output.write('%s\t%d\n' % (fam, ab))
    viral_raw_fam_ab_output.close()

    #viral_norm_fam_ab_file = sample + '.fam_norm_abundance.txt'
    #viral_norm_fam_ab_output = open(viral_norm_fam_ab_file, 'w')    
    #viral_norm_fam_ab_output.write('ID\t%s\n' % sample)         
    #for fam, ab in sample_viral_fam_dict.items():
    #    norm_ab = float(ab) / float(count_total_reads)
    #    viral_norm_fam_ab_output.write('%s\t%e\n' % (fam, norm_ab))
    #viral_norm_fam_ab_output.close()

    
    genus_count = 0
    viral_raw_genus_ab_file = sample + '.genus_raw_abundance.txt'
    viral_raw_genus_ab_output = open(viral_raw_genus_ab_file, 'w')
    viral_raw_genus_ab_output.write('ID\t%s\n' % sample)     
    for genus, ab in sample_viral_genus_dict.items():
        genus_count += 1
        viral_raw_genus_ab_output.write('%s\t%d\n' % (genus, ab))
    viral_raw_genus_ab_output.close()

    count_fam_reads = (count_fam_reads +
            (sum(sample_viral_fam_dict.values())))
    count_genus_reads = (count_genus_reads + 
            (sum(sample_viral_genus_dict.values())))

    
    output = ('Sample: %s\nReads - %d\nFamilies: %s\nGenera: %s\n\n' 
            % (sample, read_count, fam_count, genus_count))
    print(output)

output = ('Number of total reads: %d\nNumber of mapped reads: %d\n'
        'Number of mapped reads with family annotation: %d\n'
        'Number of mapped reads with genus annotation: %d') % (
        count_total_reads, count_mapped_reads, count_fam_reads,
        count_genus_reads)
print(output)


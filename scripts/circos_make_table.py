#! /usr/bin/env python

"""

"""

import argparse
import sys
import os



parser = argparse.ArgumentParser(description='script to generate table containing information on how reads are shared between samples for making circos plot')
parser.add_argument('-d', '--directory', required=True, help='directory containing files to analyze, looks for files with sub_cat.fna.comb suffix.', type=str)
parser.add_argument('-c', '--config', required=True, help='files containing infomration regarding sample metadata', type=str)
parser.add_argument('-s', '--suffix', required=True, help='suffix to use to match files for analysis', type=str)
parser.add_argument('-i', '--indiv', required=True, help='use this argument to include individuals from config file, T or F input', type=str)
args = parser.parse_args()

khmer_files_list = [filename for filename in os.listdir(args.directory) if filename.endswith(args.suffix)]

config_file = open(args.config, 'rU')
line_count = 0;
for line in config_file:
    stripline = line.rstrip()
    line_count +=1
    if line_count == 2:
        sample_list = stripline.split()
    if line_count == 4:
        diet_list = stripline.split()
    if line_count == 5:
        if args.indiv == 'T':
            indiv_list = stripline.split()

if line_count < 4:
    sys.exit('\n\nError: Not enough lines in config file, likely constructed incorrectly.\nShould contain five lines:\n' \
    '\tlines 1-3: as defined by khmer when generating .comb files\n' \
    '\tline 4: diet corresponding to sample\n' \
    '\tline 5 (optional): id for indivdual corresponding to sample\n\n')

id_list = list()
for i in range(0, len(diet_list)):  
    if args.indiv == 'T':
        cur_id = diet_list[i] + '-' + indiv_list[i]
    else:
        cur_id = diet_list[i]
    id_list.append(cur_id)


pairwise_dict = {}
for filename in khmer_files_list:
    f = open(args.directory + filename, 'rU')
    prefix = filename.split('.comb')

    sample_count = 0
    for sample in sample_list:
        if sample == prefix[0]:
            sample_column = sample_count 
        sample_count += 1
    
    cur_id = id_list[sample_column]
    for line in f:
        stripline = line.rstrip()
        columns = stripline.split()
        column_count = -2
        for column in columns: 
            column_count += 1
            if column_count != -1:
                column_value = int(column)
                if column_value > 1:
                    comp_id = id_list[column_count]
                    key_tup = (cur_id, comp_id)
                    if key_tup in pairwise_dict:
                        pairwise_dict[key_tup] += 1
                    else:
                        pairwise_dict[key_tup] = 1
         

    f.close()


for key_tup in pairwise_dict:
    print key_tup[0], '\t', key_tup[1], '\t', pairwise_dict[key_tup]


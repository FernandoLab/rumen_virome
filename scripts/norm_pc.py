#! /usr/bin/env python


"""Normalize abundance of viral and microbial PCs by sequencing depth.
"""

import pandas as pd

pc_df = pd.read_table('viral_ion_illumina_pc_cluster_raw_abundance.txt', header=0, index_col=0)
pc_df.drop(['V4_illumina','V9_illumina','V10_illumina','V11_illumina','V12_illumina','V13_illumina'], axis=1, inplace=True)
mapping = pd.read_table('viral_mapping.txt', header=0)

pc_df.to_csv('viral_pc_cluster_raw_abundance.txt', sep='\t')


pc_df['V1'] = pc_df['V1'] / int(mapping.loc[mapping['#SampleID'] == 'V1']['reads'])
pc_df['V2'] = pc_df['V2'] / int(mapping.loc[mapping['#SampleID'] == 'V2']['reads'])
pc_df['V3'] = pc_df['V3'] / int(mapping.loc[mapping['#SampleID'] == 'V3']['reads'])
pc_df['V4'] = pc_df['V4'] / int(mapping.loc[mapping['#SampleID'] == 'V4']['reads'])
pc_df['V5'] = pc_df['V5'] / int(mapping.loc[mapping['#SampleID'] == 'V5']['reads'])
pc_df['V6'] = pc_df['V6'] / int(mapping.loc[mapping['#SampleID'] == 'V6']['reads'])
pc_df['V7'] = pc_df['V7'] / int(mapping.loc[mapping['#SampleID'] == 'V7']['reads'])
pc_df['V8'] = pc_df['V8'] / int(mapping.loc[mapping['#SampleID'] == 'V8']['reads'])
pc_df['V9'] = pc_df['V9'] / int(mapping.loc[mapping['#SampleID'] == 'V9']['reads'])
pc_df['V10'] = pc_df['V10'] / int(mapping.loc[mapping['#SampleID'] == 'V10']['reads'])
pc_df['V11'] = pc_df['V11'] / int(mapping.loc[mapping['#SampleID'] == 'V11']['reads'])
pc_df['V12'] = pc_df['V12'] / int(mapping.loc[mapping['#SampleID'] == 'V12']['reads'])
pc_df['V13'] = pc_df['V13'] / int(mapping.loc[mapping['#SampleID'] == 'V13']['reads'])
pc_df['V14'] = pc_df['V14'] / int(mapping.loc[mapping['#SampleID'] == 'V14']['reads'])
pc_df['V15'] = pc_df['V15'] / int(mapping.loc[mapping['#SampleID'] == 'V15']['reads'])


pc_df.to_csv('viral_pc_cluster_norm_abundance.txt', sep='\t')


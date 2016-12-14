"""Find abundance of hits in blast search
"""

import pandas as pd
import os
import sys
import argparse
import re


viral_pc_df = pd.read_table('viral_ion_illumina_pc_cluster_raw_abundance.txt', header=0, index_col=0, na_values='0')
microbial_pc_df = pd.read_table('microbial_pc_cluster_raw_abundance.txt', header=0, index_col=0, na_values='0')

viral_pc_df.drop(['V4_illumina','V9_illumina','V10_illumina','V11_illumina','V12_illumina','V13_illumina'], axis=1, inplace=True)



viral_core_df = viral_pc_df.dropna(how='any')
viral_core_df.to_csv('viral_core_pc_table.txt', sep='\t')
viral_num_core = len(viral_core_df.index)
viral_num_pan = len(viral_pc_df.index) - viral_num_core

output = 'Number of Total Viral PCs: %d' % len(viral_pc_df.index)
print(output)
output = 'Number of Core Viral PCs: %d' % viral_num_core
print(output)
output = 'Number of Pan Viral PCs: %d' % viral_num_pan
print(output)

microbial_core_df = microbial_pc_df.dropna(how='any')
microbial_core_df.to_csv('microbial_core_pc_table.txt', sep='\t')
microbial_num_core = len(microbial_core_df.index)
microbial_num_pan = len(microbial_pc_df.index) - microbial_num_core

output = 'Number of Total Microbial PCs: %d' % len(microbial_pc_df.index)
print(output)
output = 'Number of Core Microbial PCs: %d' % microbial_num_core
print(output)
output = 'Number of Pan Microbial PCs: %d' % microbial_num_pan
print(output)


viral_cs = viral_pc_df[['V5','V9','V11','V12']]
viral_cs_core = viral_cs.dropna(how='any')
viral_cs_core.to_csv('viral_55cs_core_pc_table.txt', sep='\t')
num_viral_cs_core = len(viral_cs_core.index)
output = '\nNumber of Core Viral PCs in 55CS: %d' % num_viral_cs_core
print(output)

viral_cds = viral_pc_df[['V4','V10','V13']]
viral_cds_core = viral_cds.dropna(how='any')
viral_cds_core.to_csv('viral_27cds_core_pc_table.txt', sep='\t')
num_viral_cds_core = len(viral_cds_core.index)
output = 'Number of Core Viral PCs in 27CDS: %d' % num_viral_cds_core
print(output)

viral_corn = viral_pc_df[['V1','V2','V8','V14']]
viral_corn_core = viral_corn.dropna(how='any')
viral_corn_core.to_csv('viral_corn_core_pc_table.txt', sep='\t')
num_viral_corn_core = len(viral_corn_core.index)
output = 'Number of Core Viral PCs in Corn: %d' % num_viral_corn_core
print(output)

viral_mdgs = viral_pc_df[['V3','V6','V7','V15']]
viral_mdgs_core = viral_mdgs.dropna(how='any')
viral_mdgs_core.to_csv('viral_40mdgs_core_pc_table.txt', sep='\t')
num_viral_mdgs_core = len(viral_mdgs_core.index)
output = 'Number of Core Viral PCs in 40MDGS: %d' % num_viral_mdgs_core
print(output)

viral_non_cs = viral_pc_df.drop(viral_cs_core.columns, axis=1)
viral_non_cs_nan = viral_non_cs.isnull().all(axis=1)
viral_cs_core_exc = viral_cs_core[viral_non_cs.isnull().all(axis=1)]
viral_cs_core_exc.to_csv('viral_55cs_core_exc_pc_table.txt', sep='\t')
num_viral_cs_core_exc = len(viral_cs_core_exc.index)
output = 'Number of Exclusive Core Viral PCs in 55CS: %d' % num_viral_cs_core_exc
print(output)

viral_non_cds = viral_pc_df.drop(viral_cds_core.columns, axis=1)
viral_non_cds_nan = viral_non_cds.isnull().all(axis=1)
viral_cds_core_exc = viral_cds_core[viral_non_cds.isnull().all(axis=1)]
viral_cds_core_exc.to_csv('viral_27cds_core_exc_pc_table.txt', sep='\t')
num_viral_cds_core_exc = len(viral_cds_core_exc.index)
output = 'Number of Exclusive Core Viral PCs in 27CDS: %d' % num_viral_cds_core_exc
print(output)

viral_non_corn = viral_pc_df.drop(viral_corn_core.columns, axis=1)
viral_non_corn_nan = viral_non_corn.isnull().all(axis=1)
viral_corn_core_exc = viral_corn_core[viral_non_corn.isnull().all(axis=1)]
viral_corn_core_exc.to_csv('viral_corn_core_exc_pc_table.txt', sep='\t')
num_viral_corn_core_exc = len(viral_corn_core_exc.index)
output = 'Number of Exclusive Core Viral PCs in Corn: %d' % num_viral_corn_core_exc
print(output)

viral_non_mdgs = viral_pc_df.drop(viral_mdgs_core.columns, axis=1)
viral_non_mdgs_nan = viral_non_mdgs.isnull().all(axis=1)
viral_mdgs_core_exc = viral_mdgs_core[viral_non_mdgs.isnull().all(axis=1)]
viral_mdgs_core_exc.to_csv('viral_40mdgs_core_exc_pc_table.txt', sep='\t')
num_viral_mdgs_core_exc = len(viral_mdgs_core_exc.index)
output = 'Number of Exclusive Core Viral PCs in 40MDGS: %d' % num_viral_mdgs_core_exc
print(output)

microbial_cs = microbial_pc_df[['B5','B9','B11','B12', 'B18']]
microbial_cs_core = microbial_cs.dropna(how='any')
microbial_cs_core.to_csv('microbial_55cs_core_pc_table.txt', sep='\t')
num_microbial_cs_core = len(microbial_cs_core.index)
output = '\nNumber of Core Microbial PCs in 55CS: %d' % num_microbial_cs_core
print(output)

microbial_cds = microbial_pc_df[['B4','B10','B13', 'B16', 'B17']]
microbial_cds_core = microbial_cds.dropna(how='any')
microbial_cds_core.to_csv('microbial_27cds_core_pc_table.txt', sep='\t')
num_microbial_cds_core = len(microbial_cds_core.index)
output = 'Number of Core Microbial PCs in 27CDS: %d' % num_microbial_cds_core
print(output)

microbial_corn = microbial_pc_df[['B1','B2','B8','B14', 'B20']]
microbial_corn_core = microbial_corn.dropna(how='any')
microbial_corn_core.to_csv('microbial_corn_core_pc_table.txt', sep='\t')
num_microbial_corn_core = len(microbial_corn_core.index)
output = 'Number of Core Microbial PCs in Corn: %d' % num_microbial_corn_core
print(output)

microbial_mdgs = microbial_pc_df[['B3','B6','B7','B15', 'B19']]
microbial_mdgs_core = microbial_mdgs.dropna(how='any')
microbial_mdgs_core.to_csv('microbial_40mdgs_core_pc_table.txt', sep='\t')
num_microbial_mdgs_core = len(microbial_mdgs_core.index)
output = 'Number of Core Microbial PCs in 40MDGS: %d' % num_microbial_mdgs_core
print(output)

microbial_non_cs = microbial_pc_df.drop(microbial_cs_core.columns, axis=1)
microbial_non_cs_nan = microbial_non_cs.isnull().all(axis=1)
microbial_cs_core_exc = microbial_cs_core[microbial_non_cs.isnull().all(axis=1)]
microbial_cs_core_exc.to_csv('microbial_55cs_core_exc_pc_table.txt', sep='\t')
num_microbial_cs_core_exc = len(microbial_cs_core_exc.index)
output = 'Number of Exclusive Core Microbial PCs in 55CS: %d' % num_microbial_cs_core_exc
print(output)

microbial_non_cds = microbial_pc_df.drop(microbial_cds_core.columns, axis=1)
microbial_non_cds_nan = microbial_non_cds.isnull().all(axis=1)
microbial_cds_core_exc = microbial_cds_core[microbial_non_cds.isnull().all(axis=1)]
microbial_cds_core_exc.to_csv('microbial_27cds_core_exc_pc_table.txt', sep='\t')
num_microbial_cds_core_exc = len(microbial_cds_core_exc.index)
output = 'Number of Exclusive Core Microbial PCs in 27CDS: %d' % num_microbial_cds_core_exc
print(output)

microbial_non_corn = microbial_pc_df.drop(microbial_corn_core.columns, axis=1)
microbial_non_corn_nan = microbial_non_corn.isnull().all(axis=1)
microbial_corn_core_exc = microbial_corn_core[microbial_non_corn.isnull().all(axis=1)]
microbial_corn_core_exc.to_csv('microbial_corn_core_exc_pc_table.txt', sep='\t')
num_microbial_corn_core_exc = len(microbial_corn_core_exc.index)
output = 'Number of Exclusive Core Microbial PCs in Corn: %d' % num_microbial_corn_core_exc
print(output)

microbial_non_mdgs = microbial_pc_df.drop(microbial_mdgs_core.columns, axis=1)
microbial_non_mdgs_nan = microbial_non_mdgs.isnull().all(axis=1)
microbial_mdgs_core_exc = microbial_mdgs_core[microbial_non_mdgs.isnull().all(axis=1)]
microbial_mdgs_core_exc.to_csv('microbial_40mdgs_core_exc_pc_table.txt', sep='\t')
num_microbial_mdgs_core_exc = len(microbial_mdgs_core_exc.index)
output = 'Number of Exclusive Core Microbial PCs in 40MDGS: %d' % num_microbial_mdgs_core_exc
print(output)
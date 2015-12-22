
"""
Get median k-mer counting of each read in multiple samples.

% scripts/get_comb_multi.py <config file> 

Use '-h' for parameter help.

For each reads input file, there is an output file, with read name and the median k-mer
counting of this read in multiple samples.

NOTE: All 'N's in the input sequences are converted to 'G's.

# format of config file:

First line, the ht files for each reads data set.
Second line, the reads data set files
Third line, the memory to use

The order of ht files and reads files do not have to be the same.
Just the count in outputfile is ordered according to the order of ht files in the first line
of config file.

Example:
$ more config.txt 
1AE.ht  1A.ht  1BE.ht  1B.ht  1CE.ht  1C.ht
1A.fna.gz  1B.fna.gz  1C.fna.gz  sample1A-Exact-1X.fna.gz  sample1B-Exact-1X.fna.gz  sample1C-Exact-1X.fna.gz
100000000

Ht files should be generated by load-into-counting.py.
Hint: use ls -1 | awk '{ ORS=" "; print; }' to get the list of files

"""
import sys, screed, os
import khmer
import argparse

###

def main():
    parser = argparse.ArgumentParser(description='Get coverage of sequences in multiple samples')

    parser.add_argument('config_file')

    args = parser.parse_args()

    file_config = args.config_file

    f_config = open(file_config, 'r')
    lines = f_config.readlines()
    f_config.close()
    print lines
    htfiles = lines[0].split()
    filenames = lines[1].split()
    mem_size = int(lines[2])
    
    round_size = 0
    all_ht = []
    round_ht = []
    for htfile in htfiles:
        size = os.path.getsize(htfile)
        print round_size
        if round_size + size < mem_size:
            round_size = round_size + size
            round_ht.append(htfile)
        else:
            all_ht.append(round_ht)
            round_size = size
            round_ht = []
            round_ht.append(htfile)
    all_ht.append(round_ht)
    print all_ht
    
    for round in range(len(all_ht)): # deal with the ht in memory
        htfiles = all_ht[round]
        hts = []
        print hts
        for htfile in htfiles:
            print 'loading counting hash from', htfile
            ht = khmer.load_counting_hash(htfile)
            K = ht.ksize()
            print K
            hts.append(ht)
            
        for n, input_filename in enumerate(filenames):
            if round >0:
                os.rename(input_filename+'.comb', input_filename+'.temp')
                temp_old = open(input_filename+'.temp','r')
            output = open(input_filename+'.comb','w')
            print 'consuming input', input_filename
            for record in screed.open(input_filename):
                seq = record.sequence.upper()
                if 'N' in seq:
                    seq = seq.replace('N', 'G')
                if round == 0:
                    to_print = record.name
                else:
                    to_print = temp_old.readline().rstrip()
                    

                for ht in hts:
                    med, _, _ = ht.get_median_count(seq)
                    to_print = to_print + ' '+str(med)
                        
                output.write(to_print+'\n')



if __name__ == '__main__':
    main()

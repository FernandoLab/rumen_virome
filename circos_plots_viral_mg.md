# Viral Metagenome Circos Plots

# Generate Tables for Circos Plots - Viral Metagenome
The code below calculates how reads are shared between samples and diets to generate tables that can be used to construct the Circos plots from the manuscript with a GUI (website - directions to do so at the bottom of this notebook). Results may differ slighly from what is observed in the manuscript due to differences in subsampling events.

## Sample-Sample Comparisons


```bash
cd prinseq_output/
 
for f in *trimm_finalQC.fna
do
  filename=$(basename "$f")
  filename="${filename%_*}"
  mothur "#sub.sample(fasta=$f, size=100000)"
done

for f in *subsample.fna
do
    filename=$(basename "$f")
    filename="${filename%.*}"
    load-into-counting.py -k 20 -N 4 -x 1e9 --report-total-kmers -s tsv "$filename""_subsample_k20.kh" $f
done


ls *subsample_k20.kh | awk '{ ORS=" "; print; }' > config_sub.txt
printf "\n" >> config_sub.txt
ls *finalQC.subsample.fna | awk '{ ORS=" "; print; }' >> config_sub.txt
printf "\n" >> config_sub.txt
printf "50000000" >> config_sub.txt

python ../scripts/get_comb_multi_old_median_kmer.py config_sub.txt

printf "\n27CDS 55CS 55CS 27CDS Corn 40MDGS Corn Corn 40MDGS 27CDS 55CS 40MDGS 40MDGS Corn 55CS" >> config_sub.txt
printf "\n259 346 3244 222 3257 259 346 3244 222 3257 259 346 3244 222 3257" >> config_sub.txt

cd ..
```

```
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.10_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 629489.
## 
## Output File Names: 
## VMG.10_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.11_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 2712098.
## 
## Output File Names: 
## VMG.11_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.12_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 470915.
## 
## Output File Names: 
## VMG.12_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.13_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 920397.
## 
## Output File Names: 
## VMG.13_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.14_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 593047.
## 
## Output File Names: 
## VMG.14_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.15_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 328660.
## 
## Output File Names: 
## VMG.15_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.1_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 593530.
## 
## Output File Names: 
## VMG.1_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.2_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 100696.
## 
## Output File Names: 
## VMG.2_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.3_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 1350530.
## 
## Output File Names: 
## VMG.3_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.4_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 1881525.
## 
## Output File Names: 
## VMG.4_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.5_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 2138335.
## 
## Output File Names: 
## VMG.5_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.6_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 276069.
## 
## Output File Names: 
## VMG.6_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.7_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 499466.
## 
## Output File Names: 
## VMG.7_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.8_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 526337.
## 
## Output File Names: 
## VMG.8_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.9_trimm_finalQC.fna, size=100000)
## Sampling 100000 from 703720.
## 
## Output File Names: 
## VMG.9_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.10_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.10_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.10_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 9140310
## saving VMG.10_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.10_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.10_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.11_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.11_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.11_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 9951144
## saving VMG.11_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.11_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.11_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.12_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.12_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.12_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 7477223
## saving VMG.12_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.12_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.12_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.13_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.13_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.13_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 8142170
## saving VMG.13_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.13_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.13_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.14_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.14_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.14_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 6159259
## saving VMG.14_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.14_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.14_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.15_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.15_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.15_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 6315264
## saving VMG.15_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.15_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.15_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.1_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.1_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.1_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 6171616
## saving VMG.1_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.1_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.1_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.2_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.2_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.2_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 2596663
## saving VMG.2_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.2_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.2_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.3_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.3_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.3_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 7186462
## saving VMG.3_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.3_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.3_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.4_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.4_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.4_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 9082962
## saving VMG.4_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.4_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.4_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.5_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.5_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.5_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 5881898
## saving VMG.5_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.5_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.5_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.6_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.6_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.6_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 4381995
## saving VMG.6_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.6_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.6_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.7_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.7_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.7_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 5987065
## saving VMG.7_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.7_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.7_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.8_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.8_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.8_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 5752438
## saving VMG.8_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.8_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.8_trimm_finalQC.subsample_subsample_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to VMG.9_trimm_finalQC.subsample_subsample_k20.kh
## Loading kmers from sequences in ['VMG.9_trimm_finalQC.subsample.fna']
## making k-mer counting table
## consuming input VMG.9_trimm_finalQC.subsample.fna
## Total number of unique k-mers: 7277492
## saving VMG.9_trimm_finalQC.subsample_subsample_k20.kh
## Writing summmary info to VMG.9_trimm_finalQC.subsample_subsample_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: VMG.9_trimm_finalQC.subsample_subsample_k20.kh.info
## ['VMG.10_trimm_finalQC.subsample_subsample_k20.kh VMG.11_trimm_finalQC.subsample_subsample_k20.kh VMG.12_trimm_finalQC.subsample_subsample_k20.kh VMG.13_trimm_finalQC.subsample_subsample_k20.kh VMG.14_trimm_finalQC.subsample_subsample_k20.kh VMG.15_trimm_finalQC.subsample_subsample_k20.kh VMG.1_trimm_finalQC.subsample_subsample_k20.kh VMG.2_trimm_finalQC.subsample_subsample_k20.kh VMG.3_trimm_finalQC.subsample_subsample_k20.kh VMG.4_trimm_finalQC.subsample_subsample_k20.kh VMG.5_trimm_finalQC.subsample_subsample_k20.kh VMG.6_trimm_finalQC.subsample_subsample_k20.kh VMG.7_trimm_finalQC.subsample_subsample_k20.kh VMG.8_trimm_finalQC.subsample_subsample_k20.kh VMG.9_trimm_finalQC.subsample_subsample_k20.kh \n', 'VMG.10_trimm_finalQC.subsample.fna VMG.11_trimm_finalQC.subsample.fna VMG.12_trimm_finalQC.subsample.fna VMG.13_trimm_finalQC.subsample.fna VMG.14_trimm_finalQC.subsample.fna VMG.15_trimm_finalQC.subsample.fna VMG.1_trimm_finalQC.subsample.fna VMG.2_trimm_finalQC.subsample.fna VMG.3_trimm_finalQC.subsample.fna VMG.4_trimm_finalQC.subsample.fna VMG.5_trimm_finalQC.subsample.fna VMG.6_trimm_finalQC.subsample.fna VMG.7_trimm_finalQC.subsample.fna VMG.8_trimm_finalQC.subsample.fna VMG.9_trimm_finalQC.subsample.fna \n', '50000000']
## 0
## 4000000328
## 4000000118
## 4000000118
## 4000000118
## 4000000158
## 4000008618
## 4000007808
## 4000069188
## 4000000118
## 4000000118
## 4000000118
## 4000050908
## 4000000118
## 4000000118
## [[], ['VMG.10_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.11_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.12_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.13_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.14_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.15_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.1_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.2_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.3_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.4_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.5_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.6_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.7_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.8_trimm_finalQC.subsample_subsample_k20.kh'], ['VMG.9_trimm_finalQC.subsample_subsample_k20.kh']]
## []
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.10_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.11_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.12_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.13_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.14_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.15_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.1_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.2_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.3_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.4_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.5_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.6_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.7_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.8_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
## []
## loading counting hash from VMG.9_trimm_finalQC.subsample_subsample_k20.kh
## 20
## consuming input VMG.10_trimm_finalQC.subsample.fna
## consuming input VMG.11_trimm_finalQC.subsample.fna
## consuming input VMG.12_trimm_finalQC.subsample.fna
## consuming input VMG.13_trimm_finalQC.subsample.fna
## consuming input VMG.14_trimm_finalQC.subsample.fna
## consuming input VMG.15_trimm_finalQC.subsample.fna
## consuming input VMG.1_trimm_finalQC.subsample.fna
## consuming input VMG.2_trimm_finalQC.subsample.fna
## consuming input VMG.3_trimm_finalQC.subsample.fna
## consuming input VMG.4_trimm_finalQC.subsample.fna
## consuming input VMG.5_trimm_finalQC.subsample.fna
## consuming input VMG.6_trimm_finalQC.subsample.fna
## consuming input VMG.7_trimm_finalQC.subsample.fna
## consuming input VMG.8_trimm_finalQC.subsample.fna
## consuming input VMG.9_trimm_finalQC.subsample.fna
```

Following script calcualtes the median shared k-mer content between each pair of samples. The script only considers reads that are shared with at least one other sample to avoid potential impact of sequencing errors on the interpretation.


```bash
perl scripts/khmer_multi_threshold.pl -khmer_multi_dir=prinseq_output -threshold=1

mv *.comb.keep.txt prinseq_output/

python scripts/circos_make_table.py -d prinseq_output/ -c prinseq_output/config_sub.txt -i T -s subsample.fna.comb.keep.txt > circos_list_sub.txt
```

```
## VMG.6_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 86904 sequences
## VMG.12_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 61983 sequences
## VMG.8_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 79607 sequences
## VMG.7_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 80437 sequences
## VMG.11_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 47348 sequences
## VMG.9_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 69913 sequences
## VMG.10_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 56901 sequences
## VMG.2_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 84456 sequences
## VMG.3_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 35342 sequences
## VMG.4_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 26875 sequences
## VMG.14_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 78813 sequences
## VMG.5_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 47836 sequences
## VMG.15_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 74167 sequences
## VMG.13_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 62734 sequences
## VMG.1_trimm_finalQC.subsample.fna.comb had 100000 sequences and kept 78518 sequences
## 
## Total number of files analyzed: 15
```


```r
options(stringsAsFactors=FALSE)
matrix_out <- data.frame(matrix(NA, nrow = 15, ncol = 15))
pw <- read.table("circos_list_sub.txt", sep="\t", header=FALSE)
r_c_names <- unique(pw$V1)
row.names(matrix_out) <- r_c_names
colnames(matrix_out) <- r_c_names

for(i in 1:nrow(pw)) {
	matrix_out[pw$V1[i], pw$V2[i]] = pw$V3[i]
}

write.table(matrix_out, file="circos_table_sub.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
```


```bash
sed -E $'s/\tCorn-3257/labels\tCorn-3257/g' circos_table_sub.txt > circos_table2_sub.txt

printf "labels\n255,211,0\n0,159,107\n255,211,0\n196,2,51\n196,2,51\n196,2,51\n255,211,0\n0,159,107\n196,2,51\n0,159,107\n0,135,189\n0,159,107\n0,135,189\n255,211,0\n0,135,189" > circos_row_colors.txt

paste circos_row_colors.txt circos_table2_sub.txt > circos_viral_all_samples_upload.txt

printf "labels\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n0,135,189\n211,211,211\n0,135,189\n211,211,211\n0,135,189" > circos_row_colors_27cds.txt
paste circos_row_colors_27cds.txt circos_table2_sub.txt > circos_viral_27cds_upload.txt

printf "labels\n211,211,211\n211,211,211\n211,211,211\n196,2,51\n196,2,51\n196,2,51\n211,211,211\n211,211,211\n196,2,51\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n211,211,211" > circos_row_colors_55cs.txt
paste circos_row_colors_55cs.txt circos_table2_sub.txt > circos_viral_55cs_upload.txt

printf "labels\n211,211,211\n0,159,107\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n0,159,107\n211,211,211\n0,159,107\n211,211,211\n0,159,107\n211,211,211\n211,211,211\n211,211,211" > circos_row_colors_40mdgs.txt
paste circos_row_colors_40mdgs.txt circos_table2_sub.txt > circos_viral_40mdgs_upload.txt

printf "labels\n255,211,0\n211,211,211\n255,211,0\n211,211,211\n211,211,211\n211,211,211\n255,211,0\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n211,211,211\n255,211,0\n211,211,211\n" > circos_row_colors_corn.txt
paste circos_row_colors_corn.txt circos_table2_sub.txt > circos_viral_corn_upload.txt
```

## Diet-Diet Comparisons

```bash
cd prinseq_output/
mkdir sub_sample_comb
mv *.comb.keep.txt sub_sample_comb
mv *subsample.fna.comb sub_sample_comb

mothur "#sub.sample(fasta=VMG.13_trimm_finalQC.fna, size=133333)"
mothur "#sub.sample(fasta=VMG.10_trimm_finalQC.fna, size=133333)"
mothur "#sub.sample(fasta=VMG.4_trimm_finalQC.fna, size=133334)"

cat VMG.4_trimm_finalQC.subsample.fna VMG.10_trimm_finalQC.subsample.fna VMG.13_trimm_finalQC.subsample.fna > 27cds_sub_cat.fna

cat VMG.14_trimm_finalQC.subsample.fna VMG.1_trimm_finalQC.subsample.fna VMG.2_trimm_finalQC.subsample.fna VMG.8_trimm_finalQC.subsample.fna > corn_sub_cat.fna

cat VMG.15_trimm_finalQC.subsample.fna VMG.3_trimm_finalQC.subsample.fna VMG.6_trimm_finalQC.subsample.fna VMG.7_trimm_finalQC.subsample.fna> 40mdgs_sub_cat.fna

cat VMG.5_trimm_finalQC.subsample.fna VMG.11_trimm_finalQC.subsample.fna VMG.12_trimm_finalQC.subsample.fna VMG.9_trimm_finalQC.subsample.fna > 55cs_sub_cat.fna

for f in *sub_cat.fna
do
    filename=$(basename "$f")
    filename="${filename%.*}"
    load-into-counting.py -k 20 -N 4 -x 1e9 --report-total-kmers -s tsv "$filename""_sub_cat_k20.kh" $f
done

ls *sub_cat_k20.kh | awk '{ ORS=" "; print; }' > config_cat_sub.txt
printf "\n" >> config_cat_sub.txt
ls *sub_cat.fna | awk '{ ORS=" "; print; }' >> config_cat_sub.txt
printf "\n" >> config_cat_sub.txt
printf "50000000" >> config_cat_sub.txt

python ../scripts/get_comb_multi_old_median_kmer.py config_cat_sub.txt

printf "\n27CDS 40MDGS 55CS Corn" >> config_cat_sub.txt

cd ..
```

```
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.13_trimm_finalQC.fna, size=133333)
## Sampling 133333 from 920397.
## 
## Output File Names: 
## VMG.13_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.10_trimm_finalQC.fna, size=133333)
## Sampling 133333 from 629489.
## 
## Output File Names: 
## VMG.10_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## [H[2J
## 
## 
## 
## 
## 
## mothur v.1.36.1
## Last updated: 7/27/2015
## 
## by
## Patrick D. Schloss
## 
## Department of Microbiology & Immunology
## University of Michigan
## pschloss@umich.edu
## http://www.mothur.org
## 
## When using, please cite:
## Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.
## 
## Distributed under the GNU General Public License
## 
## Type 'help()' for information on the commands that are available
## 
## Type 'quit()' to exit program
## 
## 
## 
## mothur > sub.sample(fasta=VMG.4_trimm_finalQC.fna, size=133334)
## Sampling 133334 from 1881525.
## 
## Output File Names: 
## VMG.4_trimm_finalQC.subsample.fna
## 
## [WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.
## 
## mothur > quit()
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to 27cds_sub_cat_sub_cat_k20.kh
## Loading kmers from sequences in ['27cds_sub_cat.fna']
## making k-mer counting table
## consuming input 27cds_sub_cat.fna
## Total number of unique k-mers: 30227960
## saving 27cds_sub_cat_sub_cat_k20.kh
## Writing summmary info to 27cds_sub_cat_sub_cat_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: 27cds_sub_cat_sub_cat_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to 40mdgs_sub_cat_sub_cat_k20.kh
## Loading kmers from sequences in ['40mdgs_sub_cat.fna']
## making k-mer counting table
## consuming input 40mdgs_sub_cat.fna
## Total number of unique k-mers: 20802890
## saving 40mdgs_sub_cat_sub_cat_k20.kh
## Writing summmary info to 40mdgs_sub_cat_sub_cat_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: 40mdgs_sub_cat_sub_cat_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to 55cs_sub_cat_sub_cat_k20.kh
## Loading kmers from sequences in ['55cs_sub_cat.fna']
## making k-mer counting table
## consuming input 55cs_sub_cat.fna
## Total number of unique k-mers: 27354463
## saving 55cs_sub_cat_sub_cat_k20.kh
## Writing summmary info to 55cs_sub_cat_sub_cat_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: 55cs_sub_cat_sub_cat_k20.kh.info
## 
## || This is the script 'load-into-counting.py' in khmer.
## || You are running khmer version 1.4.1
## || You are also using screed version 0.8
## ||
## || If you use this script in a publication, please cite EACH of the following:
## ||
## ||   * MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
## ||   * Q Zhang et al., http://dx.doi.org/10.1371/journal.pone.0101271
## ||   * A. Döring et al. http://dx.doi.org:80/10.1186/1471-2105-9-11
## ||
## || Please see http://khmer.readthedocs.org/en/latest/citations.html for details.
## 
## 
## PARAMETERS:
##  - kmer size =    20 		(-k)
##  - n tables =     4 		(-N)
##  - min tablesize = 1e+09 	(-x)
## 
## Estimated memory usage is 4e+09 bytes (n_tables x min_tablesize)
## --------
## Saving k-mer counting table to corn_sub_cat_sub_cat_k20.kh
## Loading kmers from sequences in ['corn_sub_cat.fna']
## making k-mer counting table
## consuming input corn_sub_cat.fna
## Total number of unique k-mers: 18923802
## saving corn_sub_cat_sub_cat_k20.kh
## Writing summmary info to corn_sub_cat_sub_cat_k20.kh.info.tsv
## fp rate estimated to be 0.000
## DONE.
## wrote to: corn_sub_cat_sub_cat_k20.kh.info
## ['27cds_sub_cat_sub_cat_k20.kh 40mdgs_sub_cat_sub_cat_k20.kh 55cs_sub_cat_sub_cat_k20.kh corn_sub_cat_sub_cat_k20.kh \n', '27cds_sub_cat.fna 40mdgs_sub_cat.fna 55cs_sub_cat.fna corn_sub_cat.fna \n', '50000000']
## 0
## 4000000378
## 4000059868
## 4000000248
## [[], ['27cds_sub_cat_sub_cat_k20.kh'], ['40mdgs_sub_cat_sub_cat_k20.kh'], ['55cs_sub_cat_sub_cat_k20.kh'], ['corn_sub_cat_sub_cat_k20.kh']]
## []
## consuming input 27cds_sub_cat.fna
## consuming input 40mdgs_sub_cat.fna
## consuming input 55cs_sub_cat.fna
## consuming input corn_sub_cat.fna
## []
## loading counting hash from 27cds_sub_cat_sub_cat_k20.kh
## 20
## consuming input 27cds_sub_cat.fna
## consuming input 40mdgs_sub_cat.fna
## consuming input 55cs_sub_cat.fna
## consuming input corn_sub_cat.fna
## []
## loading counting hash from 40mdgs_sub_cat_sub_cat_k20.kh
## 20
## consuming input 27cds_sub_cat.fna
## consuming input 40mdgs_sub_cat.fna
## consuming input 55cs_sub_cat.fna
## consuming input corn_sub_cat.fna
## []
## loading counting hash from 55cs_sub_cat_sub_cat_k20.kh
## 20
## consuming input 27cds_sub_cat.fna
## consuming input 40mdgs_sub_cat.fna
## consuming input 55cs_sub_cat.fna
## consuming input corn_sub_cat.fna
## []
## loading counting hash from corn_sub_cat_sub_cat_k20.kh
## 20
## consuming input 27cds_sub_cat.fna
## consuming input 40mdgs_sub_cat.fna
## consuming input 55cs_sub_cat.fna
## consuming input corn_sub_cat.fna
```


```bash
perl scripts/khmer_multi_threshold.pl -khmer_multi_dir=prinseq_output -threshold=1

mv *.comb.keep.txt prinseq_output/

python scripts/circos_make_table.py -d prinseq_output/ -c prinseq_output/config_cat_sub.txt -i F -s sub_cat.fna.comb.keep.txt > circos_list_sub_cat.txt

mkdir prinseq_output/sub_diet_comb
mv prinseq_output/*sub_cat.fna.comb.keep.txt prinseq_output/sub_diet_comb
```

```
## corn_sub_cat.fna.comb had 400000 sequences and kept 333310 sequences
## 55cs_sub_cat.fna.comb had 400000 sequences and kept 236703 sequences
## sub_sample_comb had 0 sequences and kept 0 sequences
## 27cds_sub_cat.fna.comb had 400000 sequences and kept 253562 sequences
## 40mdgs_sub_cat.fna.comb had 400000 sequences and kept 321958 sequences
## 
## Total number of files analyzed: 5
```


```r
options(stringsAsFactors=FALSE)
matrix_out <- data.frame(matrix(NA, nrow = 4, ncol = 4))
pw <- read.table("circos_list_sub_cat.txt", sep="\t", header=FALSE)
r_c_names <- unique(pw$V1)
row.names(matrix_out) <- r_c_names
colnames(matrix_out) <- r_c_names

for(i in 1:nrow(pw)) {
	matrix_out[pw$V1[i], pw$V2[i]] = pw$V3[i]
}

write.table(matrix_out, file="circos_table_sub_cat.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
```


```bash
sed -E $'s/\40MDGS/labels\40MDGS/g' circos_table_sub_cat.txt > circos_table_sub_cat2.txt

printf "labels\n0,159,107\n196,2,51\n0,135,189\n255,211,0" > circos_row_colors_cat.txt

paste circos_row_colors_cat.txt circos_table_sub_cat2.txt > circos_table_sub_cat3.txt
```

## Circos GUI

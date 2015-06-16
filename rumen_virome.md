Impact of dietary change on auxillary viral metabolism and viral-bacterial interactions in the rumen of beef cattle
===============
Author: Christopher L. Anderson (canderson30@unl.edu)


##Introduction
To recreate the analysis from Anderson et al. manuscript.  All commands below were done in a linux enviornment and memory intensive commands were carried out on a high performance cluster at UNL.


##Clone the repository
Clone the repository to get some intermediate files and certain scripts for the analysis:


    git clone https://github.com/chrisLanderson/rumen_virome.git
    cd rumen_virome

##Retrieve Raw Data
After seqeuncing the Torrent Server (version?) software demulitplexed samples, trimed off adaptors and barcodes, and removed reads less than 100 basepairs. More precisley the following commands were used within the Torrent Server:

--barcode-mode 1 --barcode-cutoff 0 --min-read-length 100 --trim-min-read-len 100


To download this raw data in FASTQ format:

    mkdir raw_viral
    mkdir raw_total
    scp canderson3@crane.unl.edu:/work/samodha/canderson3/raw_viral.tgz raw_viral/
    scp canderson3@crane.unl.edu:/work/samodha/canderson3/raw_total.tgz raw_total/
    tar -zxvf raw_viral.tgz
    tar -zxvf raw_total.tgz

##Retrieve Intermidiate Files
In order to skip some long computational steps, you can use the outputs in this directory as you go along. The output names corresponding to the intermediate files in this directory will be provided after the commands if available.

    scp -r canderson3@crane.unl.edu:/work/samodha/canderson3/interm ./

##Viral Metagenome QC:
Due to the different styles of library prep (transposon vs linker amplified), the QC steps for the viral metagenomes and total metagenomes are a bit different. 

To ensure the data was properly trimmed by the Torrent Server and to deal with known biases such as artifical duplications created by the transposon library preps used for the viral metagenomes, we performed additional QC steps.


We use Trimmomatic for adaptor and transposon removal.  To get the software:

    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
    unzip Trimmomatic-0.33.zip 

The file transp_adapt_remove.cat.txt (in intermediate file directory, iterm/transp_adapt_remove.cat.txt) was supplied to direct removal of adaptors, barcodes, and transposon seqeunces. Further, observed GC and k-mer bias was noted in the 5' and 3' ends.  In order to begin to alleviate those biases, the first 20 bp of the read was trimmed and a sample dependent trimming from the 3' end was done (more 3' trimming later).  To trim all samples:

    mkdir trimmomatic_output
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.1.fastq trimmomatic_output/VMG.1_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.2.fastq trimmomatic_output/VMG.2_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.3.fastq trimmomatic_output/VMG.3_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.4.fastq trimmomatic_output/VMG.4_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.5.fastq trimmomatic_output/VMG.5_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:100 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.6.fastq trimmomatic_output/VMG.6_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.7.fastq trimmomatic_output/VMG.7_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.8.fastq trimmomatic_output/VMG.8_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.9.fastq trimmomatic_output/VMG.9_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.10.fastq trimmomatic_output/VMG.10_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.11.fastq trimmomatic_output/VMG.11_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.12.fastq trimmomatic_output/VMG.12_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.13.fastq trimmomatic_output/VMG.13_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.14.fastq trimmomatic_output/VMG.14_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75
	java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.15.fastq trimmomatic_output/VMG.15_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 CROP:175 MINLEN:75

To remove artifical duplicates we use cd-hit-454.  We noticed early on the transposon library preps with the Ion Torrent Museek kits introduces (unfortunately sometimes many) duplicates (associated with viral sample?). While the software was originally designed for 454 sequencing, Ion Torrent shares many of the same error profiles as 454 and cd-hit-454 appears to work fairly well for removing duplicates in our data.

First download cd-hit-454:

    wget https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz
    tar -xvf cd-hit-v4.6.1-2012-08-27.tgz
    cd cd-hit-v4.6.1-2012-08-27
    make openmp=yes
    cd ..

Then remove the duplicates from all viral metagenome samples and write output to cd_hit_454_output:

    export OMP_NUM_THREADS=10

    mkdir cd_hit_454_output
    cd trimmomatic_output
    for f in *.fastq
    do
        filename=$(basename "$f")
        filename="${filename%.*}"
        cd ~
        cd-hit-v4.6.1-2012-08-27/./cd-hit-454 -i trimmomatic_output/$f -o cd_hit_454_output/"$filename""_cd454.fastq" -M 6100 -T 10
    done


Now, we wanted to check for other types of duplicates and some that may have been missed by cd-hit-454. Further, we can use prinseq to remove transposon associated seqeunces that kept showing up in the 3' end.  We determined that trimming the last 5% of a read helped a lot.


Download prinseq-lite:

    wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz
    tar -xvf prinseq-lite-0.20.4.tar.gz

Run prinseq on each viral metagenome:

    mkdir prinseq_output
    cd cd_hit_454_output
    for f in *.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        cd ~
        perl prinseq-lite-0.20.4/./prinseq-lite.pl -custom_params "TGAACTG 1;GAACTGA 1;AACTGAC 1;ACTGACG 1;CTGACGC 1;TGACGCA 1;GACGCAC 1;ACGCACG 1;CGCACGA 1;GCACGAA 1;" -derep 14 -lc_method dust -lc_threshold 7 -fastq cd_hit_454_output/$f -out_format 3 -out_good prinseq_output/"$filename""_prinseq"
    done


##Total Metagenome QC:
For the total metagenome, we observe less artifical duplications and bias due to the more standard linker amplification library prep.  As a result, we had fewer QC steps.  First, to ensure removal of adaptors and barcodes for Torrent server:

    mkdir trimmomatic_output_total
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.1.fastq trimmomatic_output_total/BMG.1_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.2.fastq trimmomatic_output_total/BMG.2_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.3.fastq trimmomatic_output_total/BMG.3_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.4.fastq trimmomatic_output_total/BMG.4_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.5.fastq trimmomatic_output_total/BMG.5_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.6.fastq trimmomatic_output_total/BMG.6_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.7.fastq trimmomatic_output_total/BMG.7_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.8.fastq trimmomatic_output_total/BMG.8_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.9.fastq trimmomatic_output_total/BMG.9_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.10.fastq trimmomatic_output_total/BMG.10_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.11.fastq trimmomatic_output_total/BMG.11_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.12.fastq trimmomatic_output_total/BMG.12_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.13.fastq trimmomatic_output_total/BMG.13_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.14.fastq trimmomatic_output_total/BMG.14_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.15.fastq trimmomatic_output_total/BMG.15_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.16.fastq trimmomatic_output_total/BMG.16_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.17.fastq trimmomatic_output_total/BMG.17_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.18.fastq trimmomatic_output_total/BMG.18_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.19.fastq trimmomatic_output_total/BMG.19_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_total/BMG.20.fastq trimmomatic_output_total/BMG.20_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/total_adapt_remove.cat.txt:2:20:10 MINLEN:85

Now remove all artificial duplicates:

    export OMP_NUM_THREADS=10

    mkdir cd_hit_454_output_total
    cd trimmomatic_output_total
    for f in *.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        cd ~
        cd-hit-v4.6.1-2012-08-27/./cd-hit-454 -i trimmomatic_output_total/$f -o cd_hit_454_output_total/"$filename""_cd454.fastq" -M 6100 -T 10
    done


##rRNA Contamination
Convert all the viral finalQC FASTQ files to fasta to use for rRNA predictor.
	
	cd prinseq_output
	for f in *_.fastq
	do
    	convert_fastaqual_fastq.py -f $f -c fastq_to_fastaqual
	done

Download the rRNA predictor:
	
	wget http://weizhong-lab.ucsd.edu/meta_rna/rRNA_prediction.tar.bz2
	bzip2 -d rRNA_prediction.tar.bz2
	tar -xvf rRNA_prediction.tar
	chmod 777 -R rRNA_prediction/
	
For some reason I can only get this to run from within the examples folder.  So, move all the fasta files generated to the examples folder within rRNA predictor.
	
	mkdir rRNA_prediction/examples/e1/viral_input
	mkdir rRNA_prediction/examples/e1/viral_output
	cp prinseq_output/*_prinseq.fna rRNA_prediction/examples/e1/input_viral

Run the predictor for all samples.  Will likely need to put the full path to get it running:	

	cd rRNA_prediction/examples/e1
	export PATH=//work/samodha/canderson3/rRNA_prediction/rRNA_hmm_fs_wst_v0:$PATH
	/work/samodha/canderson3/rRNA_prediction/scripts/rRNA_hmm_run_wst_v0.pl  input_viral output_viral

Check the number of rRNA hits for each sample using custom script to parse outputs:

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.1_finalQC.fna.coord 

Total number of rRNA: 149
Total number of prokaryptic 16S and 23S rRNA: 74

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.2_finalQC.fna.coord 

Total number of rRNA: 6
Total number of prokaryptic 16S and 23S rRNA: 5

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.3_finalQC.fna.coord 

Total number of rRNA: 125
Total number of prokaryptic 16S and 23S rRNA: 112

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.4_finalQC.fna.coord 

Total number of rRNA: 828
Total number of prokaryptic 16S and 23S rRNA: 760

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.5_finalQC.fna.coord 

Total number of rRNA: 1030
Total number of prokaryptic 16S and 23S rRNA: 1010

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.6_finalQC.fna.coord 

Total number of rRNA: 38
Total number of prokaryptic 16S and 23S rRNA: 29

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.7_finalQC.fna.coord 

Total number of rRNA: 27
Total number of prokaryptic 16S and 23S rRNA: 23

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.8_finalQC.fna.coord 

Total number of rRNA: 139
Total number of prokaryptic 16S and 23S rRNA: 117

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.9_finalQC.fna.coord 

Total number of rRNA: 249
Total number of prokaryptic 16S and 23S rRNA: 244

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.10_finalQC.fna.coord 

Total number of rRNA: 579
Total number of prokaryptic 16S and 23S rRNA: 513

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.11_finalQC.fna.coord 

Total number of rRNA: 3110
Total number of prokaryptic 16S and 23S rRNA: 2981

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.12_finalQC.fna.coord 

Total number of rRNA: 428
Total number of prokaryptic 16S and 23S rRNA: 393

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.13_finalQC.fna.coord 

Total number of rRNA: 642
Total number of prokaryptic 16S and 23S rRNA: 625

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.14_finalQC.fna.coord 

Total number of rRNA: 96
Total number of prokaryptic 16S and 23S rRNA: 82

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.15_finalQC.fna.coord 

Total number of rRNA: 358
Total number of prokaryptic 16S and 23S rRNA: 353


Cat all the .coord and fasta files together to get a total count for each:

	cat rRNA_prediction/examples/e1/input_viral/*.fna > VMG.cat.fasta
	grep -c ">" VMG.cat.fasta
13232896
	
	cat rRNA_prediction/examples/e1/output_viral/*.coord > rRNA_prediction/examples/e1/output_viral/VMG.cat.coord
	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/output_viral/VMG.cat.coord 

Total number of rRNA: 7818

Total number of prokaryptic 16S and 23S rRNA: 7321
	
7818 total rRNA detected / 13232896 = 0.059%

7321 prokaryotic rRNA detected / 13232896 = 0.055%

##K-mer Profiles
Get khmer software to compare k-mer profiles of different metagenomes:

    cd virtualenv-*; python virtualenv.py ../khmerEnv; cd ..
    source khmerEnv/bin/activate
    pip install khmer

Now, load all viral metagenomes into counting:

    source khmerEnv/bin/activate
    cd /work/samodha/canderson3/prinseq_output
    for f in *_finalQC.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        load-into-counting.py -k 20 -N 4 -x 8e9 --report-total-kmers -s tsv "$filename""_k20.ct" $f
    done


Now, load all total metagenomes into counting:

    source khmerEnv/bin/activate
    cd cd_hit_454_output_total
    for f in *.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        load-into-counting.py -k 20 -N 4 -x 8e9 --report-total-kmers -s tsv "$filename""_k20.ct" $f
    done
        

Now, need to look at how reads are shared across samples in a pairwise fashion within each dataset. Currently, we need to fetch the script from authors github and do some altering to it:

	wget https://raw.githubusercontent.com/qingpeng/igs-diversity/master/scripts/get_comb_multi.py

To get it to work in our current khmer environment, use vi and alter the first line to be as follows, minus the single quotes (will need to change the path to match yours):

'#!/lustre/work/samodha/canderson3/khmerEnv/bin/python'
'#remove if K <= line and move everything over??

	mv get_comb_multi.py khmerEnv/bin/
	chmod 775 khmerEnv/bin/get_comb_multi.py 
	
Viral:

Setup the config file used by the get_comb_multi.py script. To do so, we need to know the average read length:
	
	cat prinseq_output/*_finalQC.fastq > VMG.cat.fastq
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fastq VMG.cat.fastq -stats_len
	
stats_len	max	150
stats_len	mean	131.98
stats_len	median	150
stats_len	min	75
stats_len	mode	150
stats_len	modeval	8748392
stats_len	range	76
stats_len	stddev	29.02
	
	cd prinseq_output
	ls *_k20.ct | awk '{ ORS=" "; print; }' > config.txt
	printf "\n" >>config.txt
	ls *_finalQC.fastq | awk '{ ORS=" "; print; }' >> config.txt
	printf "\n" >>config.txt
	printf "50000000" >> config.txt
	printf ""
	printf ""

	source khmerEnv/bin/activate
	get_comb_multi.py config.txt


Total:

Setup the config file used by the get_comb_multi.py script. To do so, we need to know the average read length:

	cat cd_hit_454_output_total/*.fastq > BMG.cat.fastq
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fastq BMG.cat.fastq -stats_len

stats_len	max	400
stats_len	mean	274.39
stats_len	median	293
stats_len	min	85
stats_len	mode	350
stats_len	modeval	2356329
stats_len	range	316
stats_len	stddev	86.91




##IGS Pipeline
Do threshold of 1 (present in one other sample), whatever that was before
THen do the pipeline with notes from test run
community Stats in R. Then look at pairwise statistics

##Assembly
Viral:

	wget http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0-Linux.tar.gz
	tar -zxvf SPAdes-3.5.0-Linux.tar.gz
	
	python SPAdes-3.5.0-Linux/bin/spades.py -k 21,33,55,77 --only-assembler --sc -s VMG.cat.fastq -m 60000 -t 2 -o vmg_cat_sc

Total:

	python SPAdes-3.5.0-Linux/bin/spades.py -k 21,33,55,77,99,127 --only-assembler --sc -s BMG.cat.fastq -m 200000 -t 2 -o bmg_cat_sc


##Predict ORFs
Download prodigal:

	wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.2/prodigal.linux
	
	wget https://github.com/hyattpd/Prodigal/archive/v2.6.2.zip
	unzip v2.6.2.zip
	cd Prodigal-2.6.2/
	#will need to use full path to install where needed
	make install INSTALLDIR=/work/samodha/canderson3/
	cd ..
	
Viral:
	
	mv vmg_cat_sc/contigs.fasta /work/samodha/canderson3/
	mv contigs.fasta vmg.contigs.fasta
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fasta vmg.contigs.fasta -min_len 250 -out_format 1 -out_good vmg.contigs.filter
	./prodigal -i vmg.contigs.filter.fasta -d vmg.orfs_nt.fasta -a vmg.orfs_aa.fasta -p meta 

Remove later:
perl_scripts/./length_distribution.pl vmg.orfs_nt.fasta 
Total reads:	268554
Total nt:	89488467
Mean length:	333.223362899082
Median length:	261
Mode length:	252
Max length:	11589
Min length:	60
Length	# Seqs
10	0
35	0
50	0
75	7219
99	12564
100	0
150	30751
200	35456
250	37120
300	43523
350	24599
400	17958
450	12368
500	8447
550	6684
600	5228
650	3646
700	3173
750	2515
800	2071
850	1824
900	1581
950	1333
1000	1197

Filter the ORFs based on length of nt, then use those sequences as a guide to remove the short seqeunces from aa file:
	
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fasta vmg.orfs_nt.fasta -min_len 200 -out_format 1 -out_good vmg.orfs_nt.filter
	
	
Download and install basic version of QIIME, used more lately

	wget https://pypi.python.org/packages/source/v/virtualenv/virtualenv-12.0.7.tar.gz#md5=e08796f79d112f3bfa6653cc10840114
    tar -xzf virtualenv-12.0.7.tar.gz
    cd virtualenv-*; python virtualenv.py ../qiimeEnv; cd ..
    source qiimeEnv/bin/activate
    pip install numpy
    pip install qiime

	source qiimeEnv/bin/activate
	filter_fasta.py -f vmg.orfs_aa.fasta -a vmg.orfs_nt.filter.fasta -o vmg.orfs_aa.filter.fasta
	



Total:

	mv vmg_cat_sc/contigs.fasta /work/samodha/canderson3/
	mv contigs.fasta bmg.contigs.fasta
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fasta vmg.contigs.fasta -min_len 250 -out_format 1 -out_good vmg.contigs.filter
	./prodigal -i bmg.contigs.filter.fasta -d bmg.orfs_nt.fasta -a bmg.orfs_aa.fasta -p meta 


perl_scripts/length_distribution.pl bmg.orfs_nt.fasta 
Total reads:	1637008
Total nt:	653989896
Mean length:	399.503176526932
Median length:	315
Mode length:	204
Max length:	11283
Min length:	60
Length	# Seqs
10	0
35	0
50	0
75	26390
99	52880
100	0
150	159350
200	175694
250	191802
300	171692
350	134710
400	121632
450	104511
500	82397
550	74355
600	61459
650	45533
700	38119
750	30394
800	23881
850	20938
900	17838
950	14061
1000	12274

	
Filter the ORFs based on length of nt, then use those sequences as a guide to remove the short seqeunces from aa file:
	
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fasta bmg.orfs_nt.fasta -min_len 200 -out_format 1 -out_good bmg.orfs_nt.filter
	
	source qiimeEnv/bin/activate
	filter_fasta.py -f bmg.orfs_aa.fasta -a bmg.orfs_nt.filter.fasta -o bmg.orfs_aa.filter.fasta
	
##Search against KEGG and PHAST VIRAL DBs

Unfortunately, we do use a licensed version of KEGG and usearch for large searches and I can not share those, but the outputs from the search are available in the interm directory. If you wish to use the free version of usearch (up to 4GB RAM), it can be downloaded from http://www.drive5.com/usearch/download.html and is straightforward to install and get running.

We prefer ublast from userach to do our large scale blasting.  Here, we will be doing protein-protein searches.  First, we need to configure the usearch datbase from the KEGG sequeneces:

	./usearch7.0.10 -makeudb_ublast keggV69.genes.pep.txt -output keggV69.genes.pep.udb

I included the exact version I used (April 15, 2015) into the interm directory. You can download the latest version though using:

	wget http://phast.wishartlab.com/phage_finder/DB/prophage_virus.db
	
I like to replace all spaces in PHAST database with underscores so usearch outputs the whole entry in the output file.  Did so manually in a text editor.

Create the userach database for PHAST as well:

	./usearch7.0.10 -makeudb_ublast prophage_virus.db -output prophage_virus.udb


Search Viral against KEGG:

	./usearch7.0.10 -ublast vmg.orfs_aa.filter.fasta -db keggV69.genes.pep.udb -evalue 1e-5 -blast6out vmg.ublast.kegg.txt -strand both -top_hits_only -threads 15


Search Viral against PHAST:

	./usearch7.0.10 -ublast vmg.orfs_aa.filter.fasta -db prophage_virus.udb -evalue 1e-5 -blast6out vmg.ublast.phast.txt -strand both -top_hits_only -threads 15


Search Total against KEGG:
	
	./usearch7.0.10 -ublast bmg.orfs_aa.filter.fasta -db keggV69.genes.pep.udb -evalue 1e-5 -blast6out vmg.ublast.kegg.txt -strand both -top_hits_only -threads 15


Search Total against PHAST:
	
	./usearch7.0.10 -ublast bmg.orfs_aa.filter.fasta -db prophage_virus.udb -evalue 1e-5 -blast6out vmg.ublast.phast.txt -strand both -top_hits_only -threads 15

##Map Reads to ORFs
Download bowtie2:

	wget http://sourceforge.net/projects/bowtie-bio/files/latest/download?source=files
	mv download?source=files bowtie-2.2.5.source.zip
	unzip bowtie-2.2.5.source.zip
	cd bowtie2-2.2.5/
	make
	cd ..

Make bowtie database:
	
	mkdir vmg_orfs_bowtie
	bowtie2-2.2.5/bowtie2-build vmg.orfs_nt.filter.fasta vmg_orfs_bowtie/vmg_orfs_bowtie_db
	mkdir bmg_orfs_bowtie
	bowtie2-2.2.5/bowtie2-build bmg.orfs_nt.filter.fasta bmg_orfs_bowtie/bmg_orfs_bowtie_db


Align all viral reads to nt ORFs:

    for f in prinseq_output/*_finalQC.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        bowtie2-2.2.5/bowtie2 -U $f --end-to-end --sensitive -x vmg_orfs_bowtie/vmg_orfs_bowtie_db -S vmg_orfs_bowtie/$filename.psl --un vmg_orfs_bowtie/$filename.unaligned.txt --al vmg_orfs_bowtie/$filename.aligned.txt 
    done
 
All all total metagenome reads to nt ORFs:
	
	for f in cd_hit_454_output_total/*_cd454.fastq
	do
		filename=$(basename "$f")
    	filename="${filename%_*}"
    	bowtie2-2.2.5/bowtie2 -U $f --end-to-end --sensitive -x bmg_orfs_bowtie/bmg_orfs_bowtie_db -S bmg_orfs_bowtie/$filename.psl --un bmg_orfs_bowtie/$filename.unaligned.txt --al bmg_orfs_bowtie/$filename.aligned.txt 
	done


##Create abundance tables
Wrote a custom perl script that calcualtes the abundance of each ORF and KEGG KO based on mapping reads to ORFs and blasting results against KEGG.  Once again, you need to use KEGG and since we are using a licensed version, I cannot provide that file, but the outputs of the script below are available.

The perl script uses the first hit (top hit) because multiple entries could be a top hit in the case that they have the same e-value. So, you could also input a blast file with multiple hits for each entry. The script will output a lot of different infomration, a lot of which we don't end up using downstream.

Viral:

	./blast2tsv.pl -blast_file=vmg.ublast.kegg.txt -bowtie_dir=vmg_orfs_bowtie -gene_ko=kegg_2/genes/ko/ko_genes.list

Number of top hits from BLAST file: 72840

Alignment count for VMG.6: 169706
Total Read count VMG.6: 253755
Redundant gene count with aligned reads for VMG.6: 3515
Non-Redundant gene count with aligned reads for VMG.6: 2363
Total abundance in gene hits from aligned reads for VMG.6: 51913
Non-Redundant ko count with aligned reads for VMG.6: 206
Number of nr genes with ko hits: 437
Total abundance in ko hits from aligned reads for VMG.6: 4736
Number of nr genes without ko hits: 1926
Total abundance of reads from genes without ko hits: 47177

Alignment count for VMG.13: 517178
Total Read count VMG.13: 894911
Redundant gene count with aligned reads for VMG.13: 16101
Non-Redundant gene count with aligned reads for VMG.13: 11115
Total abundance in gene hits from aligned reads for VMG.13: 172042
Non-Redundant ko count with aligned reads for VMG.13: 1481
Number of nr genes with ko hits: 3777
Total abundance in ko hits from aligned reads for VMG.13: 26669
Number of nr genes without ko hits: 7338
Total abundance of reads from genes without ko hits: 145373

Alignment count for VMG.8: 318319
Total Read count VMG.8: 499499
Redundant gene count with aligned reads for VMG.8: 4748
Non-Redundant gene count with aligned reads for VMG.8: 3360
Total abundance in gene hits from aligned reads for VMG.8: 108708
Non-Redundant ko count with aligned reads for VMG.8: 317
Number of nr genes with ko hits: 715
Total abundance in ko hits from aligned reads for VMG.8: 11198
Number of nr genes without ko hits: 2645
Total abundance of reads from genes without ko hits: 97510

Alignment count for VMG.4: 1061166
Total Read count VMG.4: 1803796
Redundant gene count with aligned reads for VMG.4: 21345
Non-Redundant gene count with aligned reads for VMG.4: 12504
Total abundance in gene hits from aligned reads for VMG.4: 431386
Non-Redundant ko count with aligned reads for VMG.4: 1243
Number of nr genes with ko hits: 4123
Total abundance in ko hits from aligned reads for VMG.4: 80943
Number of nr genes without ko hits: 8381
Total abundance of reads from genes without ko hits: 350443

Alignment count for VMG.2: 52531
Total Read count VMG.2: 83104
Redundant gene count with aligned reads for VMG.2: 1109
Non-Redundant gene count with aligned reads for VMG.2: 818
Total abundance in gene hits from aligned reads for VMG.2: 15401
Non-Redundant ko count with aligned reads for VMG.2: 73
Number of nr genes with ko hits: 122
Total abundance in ko hits from aligned reads for VMG.2: 557
Number of nr genes without ko hits: 696
Total abundance of reads from genes without ko hits: 14844

Alignment count for VMG.14: 340451
Total Read count VMG.14: 558324
Redundant gene count with aligned reads for VMG.14: 6012
Non-Redundant gene count with aligned reads for VMG.14: 3782
Total abundance in gene hits from aligned reads for VMG.14: 151391
Non-Redundant ko count with aligned reads for VMG.14: 308
Number of nr genes with ko hits: 765
Total abundance in ko hits from aligned reads for VMG.14: 17462
Number of nr genes without ko hits: 3017
Total abundance of reads from genes without ko hits: 133929

Alignment count for VMG.1: 355977
Total Read count VMG.1: 555195
Redundant gene count with aligned reads for VMG.1: 7609
Non-Redundant gene count with aligned reads for VMG.1: 4622
Total abundance in gene hits from aligned reads for VMG.1: 144310
Non-Redundant ko count with aligned reads for VMG.1: 432
Number of nr genes with ko hits: 1007
Total abundance in ko hits from aligned reads for VMG.1: 28557
Number of nr genes without ko hits: 3615
Total abundance of reads from genes without ko hits: 115753

Alignment count for VMG.10: 317750
Total Read count VMG.10: 622454
Redundant gene count with aligned reads for VMG.10: 13480
Non-Redundant gene count with aligned reads for VMG.10: 9812
Total abundance in gene hits from aligned reads for VMG.10: 115941
Non-Redundant ko count with aligned reads for VMG.10: 1491
Number of nr genes with ko hits: 3473
Total abundance in ko hits from aligned reads for VMG.10: 17950
Number of nr genes without ko hits: 6339
Total abundance of reads from genes without ko hits: 97991

Alignment count for VMG.12: 214299
Total Read count VMG.12: 448057
Redundant gene count with aligned reads for VMG.12: 9853
Non-Redundant gene count with aligned reads for VMG.12: 7295
Total abundance in gene hits from aligned reads for VMG.12: 58612
Non-Redundant ko count with aligned reads for VMG.12: 992
Number of nr genes with ko hits: 2625
Total abundance in ko hits from aligned reads for VMG.12: 9671
Number of nr genes without ko hits: 4670
Total abundance of reads from genes without ko hits: 48941

Alignment count for VMG.3: 894624
Total Read count VMG.3: 1280825
Redundant gene count with aligned reads for VMG.3: 12409
Non-Redundant gene count with aligned reads for VMG.3: 6657
Total abundance in gene hits from aligned reads for VMG.3: 332528
Non-Redundant ko count with aligned reads for VMG.3: 517
Number of nr genes with ko hits: 1480
Total abundance in ko hits from aligned reads for VMG.3: 40368
Number of nr genes without ko hits: 5177
Total abundance of reads from genes without ko hits: 292160

Alignment count for VMG.9: 439638
Total Read count VMG.9: 673709
Redundant gene count with aligned reads for VMG.9: 12439
Non-Redundant gene count with aligned reads for VMG.9: 7948
Total abundance in gene hits from aligned reads for VMG.9: 177056
Non-Redundant ko count with aligned reads for VMG.9: 892
Number of nr genes with ko hits: 2526
Total abundance in ko hits from aligned reads for VMG.9: 23309
Number of nr genes without ko hits: 5422
Total abundance of reads from genes without ko hits: 153747

Alignment count for VMG.11: 1566769
Total Read count VMG.11: 2667035
Redundant gene count with aligned reads for VMG.11: 39741
Non-Redundant gene count with aligned reads for VMG.11: 20721
Total abundance in gene hits from aligned reads for VMG.11: 621550
Non-Redundant ko count with aligned reads for VMG.11: 1490
Number of nr genes with ko hits: 8225
Total abundance in ko hits from aligned reads for VMG.11: 138390
Number of nr genes without ko hits: 12496
Total abundance of reads from genes without ko hits: 483160

Alignment count for VMG.15: 185755
Total Read count VMG.15: 312234
Redundant gene count with aligned reads for VMG.15: 5956
Non-Redundant gene count with aligned reads for VMG.15: 4398
Total abundance in gene hits from aligned reads for VMG.15: 78894
Non-Redundant ko count with aligned reads for VMG.15: 776
Number of nr genes with ko hits: 1439
Total abundance in ko hits from aligned reads for VMG.15: 12678
Number of nr genes without ko hits: 2959
Total abundance of reads from genes without ko hits: 66216

Alignment count for VMG.7: 306421
Total Read count VMG.7: 468317
Redundant gene count with aligned reads for VMG.7: 5508
Non-Redundant gene count with aligned reads for VMG.7: 3403
Total abundance in gene hits from aligned reads for VMG.7: 109234
Non-Redundant ko count with aligned reads for VMG.7: 233
Number of nr genes with ko hits: 628
Total abundance in ko hits from aligned reads for VMG.7: 12777
Number of nr genes without ko hits: 2775
Total abundance of reads from genes without ko hits: 96457

Alignment count for VMG.5: 1181469
Total Read count VMG.5: 2111681
Redundant gene count with aligned reads for VMG.5: 20957
Non-Redundant gene count with aligned reads for VMG.5: 12000
Total abundance in gene hits from aligned reads for VMG.5: 534394
Non-Redundant ko count with aligned reads for VMG.5: 1170
Number of nr genes with ko hits: 3997
Total abundance in ko hits from aligned reads for VMG.5: 112021
Number of nr genes without ko hits: 8003
Total abundance of reads from genes without ko hits: 422373


Total nucleotides for all samples combined: 1746518841
Average nucleotides per samples: 116434589.4
Number of ORFs with assignments printed to orf_assignment.txt: 24232

	mkdir viral_blast2tsv
	mv *no_ko_hit.txt viral_blast2tsv/
	mv *_abundances.txt viral_blast2tsv/
	mv *_abundance.txt viral_blast2tsv/
	mv blast_top_hit.txt viral_blast2tsv/
	mv orf_ko_assignment.txt viral_blast2tsv/


Merge all the normalized abundance tables together:
	
	source qiimeEnv/bin/activate
	for f in viral_blast2tsv/*ko_corrected_abundances.txt
	do
		filename=$(basename "$f")
		filename="${filename%.txt}"
		biom convert --table-type="OTU table" --to-json -i $f -o viral_blast2tsv/$filename.biom
	done

	merge_otu_tables.py -i viral_blast2tsv/VMG.1.ko_corrected_abundances.biom,viral_blast2tsv/VMG.2.ko_corrected_abundances.biom,viral_blast2tsv/VMG.3.ko_corrected_abundances.biom,viral_blast2tsv/VMG.4.ko_corrected_abundances.biom,viral_blast2tsv/VMG.5.ko_corrected_abundances.biom,viral_blast2tsv/VMG.6.ko_corrected_abundances.biom,viral_blast2tsv/VMG.7.ko_corrected_abundances.biom,viral_blast2tsv/VMG.8.ko_corrected_abundances.biom,viral_blast2tsv/VMG.9.ko_corrected_abundances.biom,viral_blast2tsv/VMG.10.ko_corrected_abundances.biom,viral_blast2tsv/VMG.11.ko_corrected_abundances.biom,viral_blast2tsv/VMG.12.ko_corrected_abundances.biom,viral_blast2tsv/VMG.13.ko_corrected_abundances.biom,viral_blast2tsv/VMG.14.ko_corrected_abundances.biom,viral_blast2tsv/VMG.15.ko_corrected_abundances.biom -o viral_blast2tsv/VMG.ko_corrected_abundances.biom

	biom summarize-table -i viral_blast2tsv/VMG.ko_corrected_abundances.biom -o viral_blast2tsv/VMG.ko_corrected_abundances_summarize.txt
	
	biom convert -i viral_blast2tsv/VMG.ko_corrected_abundances.biom -o viral_blast2tsv/VMG.ko_corrected_abundances.txt --table-type="OTU table" --to-tsv


Total:

	./blast2tsv.pl -blast_file=bmg.ublast.kegg.txt -bowtie_dir=bmg_orfs_bowtie -gene_ko=kegg_2/genes/ko/ko_genes.list
	
Number of top hits from BLAST file: 830683

Alignment count for BMG.1: 468360
Total Read count BMG.1: 931099
Redundant gene count with aligned reads for BMG.1: 108316
Non-Redundant gene count with aligned reads for BMG.1: 52631
Total abundance in gene hits from aligned reads for BMG.1: 401605
Non-Redundant ko count with aligned reads for BMG.1: 2707
Number of nr genes with ko hits: 25688
Total abundance in ko hits from aligned reads for BMG.1: 173612
Number of nr genes without ko hits: 26943
Total abundance of reads from genes without ko hits: 227993

Alignment count for BMG.10: 309109
Total Read count BMG.10: 631225
Redundant gene count with aligned reads for BMG.10: 98022
Non-Redundant gene count with aligned reads for BMG.10: 55013
Total abundance in gene hits from aligned reads for BMG.10: 265551
Non-Redundant ko count with aligned reads for BMG.10: 2685
Number of nr genes with ko hits: 27803
Total abundance in ko hits from aligned reads for BMG.10: 131065
Number of nr genes without ko hits: 27210
Total abundance of reads from genes without ko hits: 134486

Alignment count for BMG.11: 1119959
Total Read count BMG.11: 2261004
Redundant gene count with aligned reads for BMG.11: 165603
Non-Redundant gene count with aligned reads for BMG.11: 68806
Total abundance in gene hits from aligned reads for BMG.11: 931976
Non-Redundant ko count with aligned reads for BMG.11: 2871
Number of nr genes with ko hits: 33467
Total abundance in ko hits from aligned reads for BMG.11: 401052
Number of nr genes without ko hits: 35339
Total abundance of reads from genes without ko hits: 530924

Alignment count for BMG.12: 770771
Total Read count BMG.12: 1591514
Redundant gene count with aligned reads for BMG.12: 133470
Non-Redundant gene count with aligned reads for BMG.12: 61010
Total abundance in gene hits from aligned reads for BMG.12: 632660
Non-Redundant ko count with aligned reads for BMG.12: 2682
Number of nr genes with ko hits: 29825
Total abundance in ko hits from aligned reads for BMG.12: 258092
Number of nr genes without ko hits: 31185
Total abundance of reads from genes without ko hits: 374568

Alignment count for BMG.13: 355193
Total Read count BMG.13: 972805
Redundant gene count with aligned reads for BMG.13: 107361
Non-Redundant gene count with aligned reads for BMG.13: 49668
Total abundance in gene hits from aligned reads for BMG.13: 289533
Non-Redundant ko count with aligned reads for BMG.13: 2681
Number of nr genes with ko hits: 23799
Total abundance in ko hits from aligned reads for BMG.13: 123909
Number of nr genes without ko hits: 25869
Total abundance of reads from genes without ko hits: 165624

Alignment count for BMG.14: 1103838
Total Read count BMG.14: 2475541
Redundant gene count with aligned reads for BMG.14: 182944
Non-Redundant gene count with aligned reads for BMG.14: 86298
Total abundance in gene hits from aligned reads for BMG.14: 912395
Non-Redundant ko count with aligned reads for BMG.14: 3229
Number of nr genes with ko hits: 43523
Total abundance in ko hits from aligned reads for BMG.14: 416138
Number of nr genes without ko hits: 42775
Total abundance of reads from genes without ko hits: 496257

Alignment count for BMG.15: 338712
Total Read count BMG.15: 734576
Redundant gene count with aligned reads for BMG.15: 89556
Non-Redundant gene count with aligned reads for BMG.15: 45650
Total abundance in gene hits from aligned reads for BMG.15: 275246
Non-Redundant ko count with aligned reads for BMG.15: 2525
Number of nr genes with ko hits: 22108
Total abundance in ko hits from aligned reads for BMG.15: 117376
Number of nr genes without ko hits: 23542
Total abundance of reads from genes without ko hits: 157870

Alignment count for BMG.16: 812199
Total Read count BMG.16: 2011286
Redundant gene count with aligned reads for BMG.16: 147022
Non-Redundant gene count with aligned reads for BMG.16: 79450
Total abundance in gene hits from aligned reads for BMG.16: 604995
Non-Redundant ko count with aligned reads for BMG.16: 3166
Number of nr genes with ko hits: 40386
Total abundance in ko hits from aligned reads for BMG.16: 296067
Number of nr genes without ko hits: 39064
Total abundance of reads from genes without ko hits: 308928

Alignment count for BMG.17: 211026
Total Read count BMG.17: 509528
Redundant gene count with aligned reads for BMG.17: 70402
Non-Redundant gene count with aligned reads for BMG.17: 44428
Total abundance in gene hits from aligned reads for BMG.17: 170837
Non-Redundant ko count with aligned reads for BMG.17: 2605
Number of nr genes with ko hits: 21933
Total abundance in ko hits from aligned reads for BMG.17: 78803
Number of nr genes without ko hits: 22495
Total abundance of reads from genes without ko hits: 92034

Alignment count for BMG.18: 419898
Total Read count BMG.18: 1240812
Redundant gene count with aligned reads for BMG.18: 113752
Non-Redundant gene count with aligned reads for BMG.18: 51085
Total abundance in gene hits from aligned reads for BMG.18: 340162
Non-Redundant ko count with aligned reads for BMG.18: 2694
Number of nr genes with ko hits: 24375
Total abundance in ko hits from aligned reads for BMG.18: 142549
Number of nr genes without ko hits: 26710
Total abundance of reads from genes without ko hits: 197613

Alignment count for BMG.19: 485087
Total Read count BMG.19: 995310
Redundant gene count with aligned reads for BMG.19: 116604
Non-Redundant gene count with aligned reads for BMG.19: 59400
Total abundance in gene hits from aligned reads for BMG.19: 408014
Non-Redundant ko count with aligned reads for BMG.19: 2790
Number of nr genes with ko hits: 29668
Total abundance in ko hits from aligned reads for BMG.19: 180907
Number of nr genes without ko hits: 29732
Total abundance of reads from genes without ko hits: 227107

Alignment count for BMG.2: 749692
Total Read count BMG.2: 1331228
Redundant gene count with aligned reads for BMG.2: 145268
Non-Redundant gene count with aligned reads for BMG.2: 73768
Total abundance in gene hits from aligned reads for BMG.2: 660700
Non-Redundant ko count with aligned reads for BMG.2: 3166
Number of nr genes with ko hits: 37811
Total abundance in ko hits from aligned reads for BMG.2: 338181
Number of nr genes without ko hits: 35957
Total abundance of reads from genes without ko hits: 322519

Alignment count for BMG.20: 408000
Total Read count BMG.20: 888555
Redundant gene count with aligned reads for BMG.20: 106270
Non-Redundant gene count with aligned reads for BMG.20: 55476
Total abundance in gene hits from aligned reads for BMG.20: 319884
Non-Redundant ko count with aligned reads for BMG.20: 2702
Number of nr genes with ko hits: 26648
Total abundance in ko hits from aligned reads for BMG.20: 130613
Number of nr genes without ko hits: 28828
Total abundance of reads from genes without ko hits: 189271

Alignment count for BMG.3: 1086151
Total Read count BMG.3: 2070496
Redundant gene count with aligned reads for BMG.3: 166442
Non-Redundant gene count with aligned reads for BMG.3: 68879
Total abundance in gene hits from aligned reads for BMG.3: 930384
Non-Redundant ko count with aligned reads for BMG.3: 2767
Number of nr genes with ko hits: 33171
Total abundance in ko hits from aligned reads for BMG.3: 387823
Number of nr genes without ko hits: 35708
Total abundance of reads from genes without ko hits: 542561

Alignment count for BMG.4: 1559650
Total Read count BMG.4: 3165966
Redundant gene count with aligned reads for BMG.4: 212337
Non-Redundant gene count with aligned reads for BMG.4: 98306
Total abundance in gene hits from aligned reads for BMG.4: 1299794
Non-Redundant ko count with aligned reads for BMG.4: 3190
Number of nr genes with ko hits: 49641
Total abundance in ko hits from aligned reads for BMG.4: 630815
Number of nr genes without ko hits: 48665
Total abundance of reads from genes without ko hits: 668979

Alignment count for BMG.5: 259823
Total Read count BMG.5: 666970
Redundant gene count with aligned reads for BMG.5: 79992
Non-Redundant gene count with aligned reads for BMG.5: 37560
Total abundance in gene hits from aligned reads for BMG.5: 212833
Non-Redundant ko count with aligned reads for BMG.5: 2280
Number of nr genes with ko hits: 17377
Total abundance in ko hits from aligned reads for BMG.5: 87361
Number of nr genes without ko hits: 20183
Total abundance of reads from genes without ko hits: 125472

Alignment count for BMG.6: 1173763
Total Read count BMG.6: 2344185
Redundant gene count with aligned reads for BMG.6: 230737
Non-Redundant gene count with aligned reads for BMG.6: 98932
Total abundance in gene hits from aligned reads for BMG.6: 1017524
Non-Redundant ko count with aligned reads for BMG.6: 3232
Number of nr genes with ko hits: 49970
Total abundance in ko hits from aligned reads for BMG.6: 490093
Number of nr genes without ko hits: 48962
Total abundance of reads from genes without ko hits: 527431

Alignment count for BMG.7: 230721
Total Read count BMG.7: 484852
Redundant gene count with aligned reads for BMG.7: 89489
Non-Redundant gene count with aligned reads for BMG.7: 52921
Total abundance in gene hits from aligned reads for BMG.7: 194331
Non-Redundant ko count with aligned reads for BMG.7: 2746
Number of nr genes with ko hits: 26741
Total abundance in ko hits from aligned reads for BMG.7: 94508
Number of nr genes without ko hits: 26180
Total abundance of reads from genes without ko hits: 99823

Alignment count for BMG.8: 469011
Total Read count BMG.8: 884957
Redundant gene count with aligned reads for BMG.8: 80567
Non-Redundant gene count with aligned reads for BMG.8: 42988
Total abundance in gene hits from aligned reads for BMG.8: 382214
Non-Redundant ko count with aligned reads for BMG.8: 2336
Number of nr genes with ko hits: 20776
Total abundance in ko hits from aligned reads for BMG.8: 149320
Number of nr genes without ko hits: 22212
Total abundance of reads from genes without ko hits: 232894

Alignment count for BMG.9: 251814
Total Read count BMG.9: 667046
Redundant gene count with aligned reads for BMG.9: 87964
Non-Redundant gene count with aligned reads for BMG.9: 49582
Total abundance in gene hits from aligned reads for BMG.9: 188350
Non-Redundant ko count with aligned reads for BMG.9: 2503
Number of nr genes with ko hits: 24125
Total abundance in ko hits from aligned reads for BMG.9: 85108
Number of nr genes without ko hits: 25457
Total abundance of reads from genes without ko hits: 103242


Total nucleotides for all samples combined: 7369829541
Average nucleotides per samples: 491321969.4
Number of ORFs with assignments printed to orf_assignment.txt: 355187	
	
	mkdir bmg_blast2tsv
	mv *no_ko_hit.txt bmg_blast2tsv/
	mv *_abundances.txt bmg_blast2tsv/
	mv *_abundance.txt bmg_blast2tsv/
	mv blast_top_hit.txt bmg_blast2tsv/
	mv orf_ko_assignment.txt bmg_blast2tsv/
	
Merge all the normalized abundance tables together:
	
	source qiimeEnv/bin/activate
	for f in bmg_blast2tsv/*ko_corrected_abundances.txt
	do
		filename=$(basename "$f")
		filename="${filename%.txt}"
		biom convert --table-type="OTU table" --to-json -i $f -o bmg_blast2tsv/$filename.biom
	done

	merge_otu_tables.py -i bmg_blast2tsv/BMG.1.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.2.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.3.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.4.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.5.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.6.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.7.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.8.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.9.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.10.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.11.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.12.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.13.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.14.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.15.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.16.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.17.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.18.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.19.ko_corrected_abundances.biom,bmg_blast2tsv/BMG.20.ko_corrected_abundances.biom -o bmg_blast2tsv/BMG.ko_corrected_abundances.biom

	biom summarize-table -i bmg_blast2tsv/BMG.ko_corrected_abundances.biom -o bmg_blast2tsv/BMG.ko_corrected_abundances_summarize.txt
	
	biom convert -i bmg_blast2tsv/BMG.ko_corrected_abundances.biom -o bmg_blast2tsv/BMG.ko_corrected_abundances.txt --table-type="OTU table" --to-tsv



##Filter Viral KOs for viral homology
We used the results from searching against the PHAST viral database to filter KOs, saying that for a KO to be considered downstream, an ORF that the KO came from had to have a singificant PHAST database hit somewhere on the contig that the ORF came from.  This helps to ensure the metabolism genes we are looking at are actually of viral origin.

NEED TO DOUBLE CHECK THIS SCRIPT!!!!???!!!!

First, need to create a list of viral KOs from abundance table just created:

	awk 'NR>2{print $1}' viral_blast2tsv/VMG.ko_corrected_abundances.txt > viral_all_kos.txt

Now, use custom perl script to filter based on criteria described above:

	./get_knum_phast_hits.pl -phast_blast=vmg.ublast.phast.txt -orf_ko_assign=viral_blast2tsv/orf_ko_assignment.txt -knum=viral_all_kos.txt

Number of unique differential k numbers found with ORFs: 2473
Total number of PHAST top hits: 30429
Contigs with phast and KO hit: 5481 
ORFs on a contig with phast hit and ORF has a KO: 6698
K numbers printed that were on a contig with PHAST hit and KO: 700

	mv orf_ko_phast_hits.txt viral_blast2tsv/
	mv ko_w_phast_hits_list.txt viral_blast2tsv

Use the ko_w_phast_hits_list.txt file to filter the KO abundance table generated previosuly:

	source qiimeEnv/bin/activate
	filter_otus_from_otu_table.py -i viral_blast2tsv/VMG.ko_corrected_abundances.biom -o viral_blast2tsv/VMG.ko_corrected_abundances_phast_filter.biom -e viral_blast2tsv/ko_w_phast_hits_list.txt --negate_ids_to_exclude
	biom convert -i viral_blast2tsv/VMG.ko_corrected_abundances_phast_filter.biom -o viral_blast2tsv/VMG.ko_corrected_abundances_phast_filter.txt --table-type="OTU table" --to-tsv



##Find Differential Features KOs:
We used LEfSe to identify differential KOs based on diet within the viral metagenome and total metagenome.  At the time, I had troubles getting LEfSe to operate via the command line so I used the version hosted on Galaxy.  If you try to recreate this step it may differ if there was changes in the LEfSe version hosted on galaxy, but the outputs are all available.

Viral:

To format the KO table for upload to run with LEfSe, I manually edited the phast filtered KO abundance table by adding Diet and Animal ID information.  The file I uploaded is found at viral_blast2tsv/VMG.ko_corrected_abundance_phast_filter_lefse.txt

Go to http://huttenhower.sph.harvard.edu/galaxy/ (hyperlink this?) to upload the file using the Get Data Module on the left and upload that file as tabular. 

Then hit LEfSe module on left side and A) Format Data for LEfSe.   Use DIET as class, no sublass, and ANIMAL as subject. Per-sample normalization should be set to 'Yes' as well. Execute the command.

Now, choose B) LDA Effect Size (LEfSe), from the left side.  Use default values and execute the command.

Click on the result on the right and download them.  Save the result to the viral_blast2tsv directory as VMG.ko_phast_lefse_output.txt



Total:

To format the KO table for upload to run with LEfSe, I manually edited the KO abundance table by adding Diet and Animal ID information.  The file I uploaded is found at bmg_blast2tsv/BMG.ko_corrected_abundance_lefse.txt

Output can be found as BMG.ko_corrected_abundance_lefse_output.txt.


##Create Metabolic Networks
We used mmnet package in R to create two metabolic netwroks: one with all the KOs found in viral and total metagenome and another that represents all of global metabolism as defined by KEGG. Because igraph does not seem to be able to calculate clustering coeffecient of a directed graph and network analyzer in cytoscape provides a few other features, we are going to use cytoscape to calcualte topological properties of the nodes.  However, I was not sure how to run cytoscape on linux, so I installed R on both linux and mac and cytoscape on mac.  I only used the mac versions to push the mmnet graphs from R to cytoscape on mac and calulcate the topological properties.  I try my best to describe that process below and as always, the inputs and outputs are available.

To further complicate things, RCytoscape is only compatible with  older versions of cytoscape, so I had to install two versions of cytoscape (2.8.2 and 3.1) and initally push the R graphs to cytoscape 2.8.2, save the graph, then open it in cytoscape 3.1 to use the most up to date ntework analyzer algorithms.

Install R-3.2.0 on mac and linux systems using same installation

Linux:

	wget http://rweb.quant.ku.edu/cran/src/base/R-3/R-3.2.0.tar.gz	
	tar -zxvf R-3.2.0.tar.gz 
	cd R-3.2.0
	./configure
	make
	cd ..

Mac:
	
	wget http://cran.r-project.org/bin/macosx/R-3.2.0.pkg

Double click on the installer and follow install directions
	

Install some libraries we are going to need on mac version:
	
	R
	source("http://bioconductor.org/biocLite.R")
	biocLite("mmnet")
	biocLite ("RCytoscape")
	library(mmnet)
	library(RCytoscape)

	
Create metabolic network for global metabolism and get topology features of each enzyme:

Cytoscape 2.8.2, install RCytoscape v.1.8
Plugins --> Cytoscape RPC --> Activate Cytoscape RPC
Accept default values, hit OK

Now, in R on the mac:
	
	data(RefDbcache)
	if (require(RCytoscape)) {
	refnet <- RefDbcache$network
	net <- igraph.to.graphNEL(refnet)
	node.attr = list.vertex.attributes(refnet)
	if (length(node.attr)) {
	node.attr.class = sapply(node.attr, class)
	node.attr.class[node.attr.class == "character"] = "char"
	for (i in 1:length(node.attr)) net <- initNodeAttribute(net,
	attribute.name = node.attr[i], attribute.type = node.attr.class[i],
	default.value = "0")
	}
	net <- initEdgeAttribute(net, attribute.name = "weight", attribute.type = "numeric", default.value = "1")
	cw <- new.CytoscapeWindow("net", graph = net, overwriteWindow = TRUE)
	displayGraph(cw)
	}


This should push the network to cytoscape 2.8.2Apply force-directed layoutSaved file as refnet.cysOpen network in cytoscape 3.1.0 (File --> open)Select the largest connected portion of the network (Tools --> network analyzer --> subnetwork creation --> extract connected components, select first component --> extract)
This should generate a network with 3616 nodes.
Analyzed subset network as directed (Tools --> network analyzer --> network analysis --> Analyze network --> treat the network as directed --> ok )
This will take awhile.File -- Export --> Table --> node(1) default node
Save as:refnet_attributes.csv
Open in excel and convert to tab delimited text file. Save as refnet_attributes.txt



Need to repeat the workflow above that was for the reference net for a network containing only the KOs found in the viral and total metagenome.

	source qiimeEnv/bin/activate
	merge_otu_tables.py -i bmg_blast2tsv/BMG.ko_corrected_abundances.biom,viral_blast2tsv/VMG.ko_corrected_abundances_phast_filter.biom -o VMG_BMG.ko_corrected_abundances.biom
	
For some reason, wont read these biom files in R, convert to a dense biom. Troubles with new biom format and quickly being able to convert from sparse to dense (needed in mmnet R package). But very easy in earlier biom format versions:

	source qiimeEnv/bin/activate
	biom convert -i VMG_BMG.ko_corrected_abundances.biom -o VMG_BMG.ko_corrected_abundances.txt --to-tsv --table-type "OTU table"
	
Only need the KO identifiers downstream and can only have one column of data entering mmnet package.
	
	awk -v OFS="\t" 'BEGIN { print "OTU" "\t" "sample" "\t" "taxonomy"} NR>2 {print $1, 100, $1}' VMG_BMG.ko_corrected_abundances.txt > VMG_BMG.ko_corrected_abundances.collapse.txt
	
	
Convert to dense biom file (have to quit qiime or new window):
	
	biom-format-1.3.1/scripts/./biom convert -i VMG_BMG.ko_corrected_abundances.collapse.txt -o VMG_BMG.ko_corrected_abundances.collapse.dense.biom --table-type "OTU table" -t dense


using the R version on mac (will need to transfer VMG_BMG.ko_corrected_abundances.dense.biom to the mac):
	
	R
	library(mmnet)
	library(RCytoscape)
	vmg_bmg <- read_biom("VMG_BMG.ko_corrected_abundances.collapse.dense.biom")
	ssn <- constructSSN(vmg_bmg)
	
	if (require(RCytoscape)) {
	net <- igraph.to.graphNEL(ssn)
	node.attr = list.vertex.attributes(ssn)
	if (length(node.attr)) {
	node.attr.class = sapply(node.attr, class)
	node.attr.class[node.attr.class == "character"] = "char"
	for (i in 1:length(node.attr)) net <- initNodeAttribute(net,
	attribute.name = node.attr[i], attribute.type = node.attr.class[i],
	default.value = "0")
	}
	net <- initEdgeAttribute(net, attribute.name = "weight", attribute.type = "numeric", default.value = "1")
	cw <- new.CytoscapeWindow("net", graph = net, overwriteWindow = TRUE)
	displayGraph(cw)
	}

Apply force-directed layoutSaved file as ssn.cysOpen network in cytoscape 3.1.0 (File --> open)Select the largest connected portion of the network (Tools --> network analyzer --> subnetwork creation --> extract connected components, select first component --> extract)
This should generate a network with 1266 nodes.
Analyzed subset network as directed (Tools --> network analyzer --> network analysis --> Analyze network --> treat the network as directed --> ok )
This will take awhile.File -- Export --> Table --> node(1) default node
Save as:ssn_attributes.csv
Open in excel and convert to tab delimited text file. Save as ssn_attributes.txt

##Topology of Differential KOs in Metabolic Network

Want to investigate the topology of differential features in both the total metagenome and viral metagenome in the context of a metabolic network.  In our metabolic network KOs (enzyme families) are nodes and their connections are established by KEGG and implemented into a network previously with the mmnet package in R.

Need to take the network attribute tables and identify the topological features for the differential KOs and compare the values of differential to non-differntial and viral differential to total differential.

	mkdir network_plots
	R-3.2.0/bin/./R
	vmg_lefse <- read.table("viral_blast2tsv/VMG.ko_phast_lefse_output.txt", sep="\t", header=FALSE)
	vmg_lefse$split <- (vmg_lefse$V4 > 1)
	vmg_diff <- split(vmg_lefse, vmg_lefse$split)
	vmg_diff <- as.data.frame(vmg_diff)
	dim(vmg_diff)
	names(vmg_diff)[1] <- "canonicalName"
	names(vmg_diff)[2] <- "KruskallWallis"
	names(vmg_diff)[3] <- "Diet"
	names(vmg_diff)[4] <- "LDA"
	names(vmg_diff)[5] <- "PValue"
	names(vmg_diff)[6] <- "Dataset"
	
	bmg_lefse <- read.table("bmg_blast2tsv/BMG.ko_corrected_abundance_lefse_output.txt", sep="\t", header=FALSE)
	bmg_lefse$split <- (bmg_lefse$V4 > 1)
	bmg_diff <- split(bmg_lefse, bmg_lefse$split)
	bmg_diff <- as.data.frame(bmg_diff)
	dim(bmg_diff)
	names(bmg_diff)[1] <- "canonicalName"
	names(bmg_diff)[2] <- "KruskallWallis"
	names(bmg_diff)[3] <- "Diet"
	names(bmg_diff)[4] <- "LDA"
	names(bmg_diff)[5] <- "PValue"
	names(bmg_diff)[6] <- "Dataset"
	
	
	ref_topo <- read.table("refnet_attributes.txt", header=TRUE, sep="\t")

Not everything that is differential is within the the meatbolic network.  So, we will lose about half of the KOs when we search the metabolic network:

	vmg_diff_topo <- merge(x=vmg_diff, y=ref_topo, by="canonicalName")
	bmg_diff_topo <- merge(x=bmg_diff, y=ref_topo, by="canonicalName")

Get KOs that are not differential in either dataset:

	diff <- rbind(vmg_diff, bmg_diff)
	total_topo <- merge(x=ref_topo, y=diff, all.x=TRUE)
	total_topo[is.na(total_topo)] <- 0
	non_diff_topo <- subset(total_topo, total_topo$LDA == 0)
	
If you look at unique names of each dataset, should equal the total number of nodes (3616).

Statistically compare topological properties of interest, but first make groups to compare:

	vmg_diff_topo[vmg_diff_topo=="TRUE"] <- "VMG"
	bmg_diff_topo[bmg_diff_topo=="TRUE"] <- "BMG"
	non_diff_topo$Dataset[non_diff_topo$Dataset == 0] <- "Non"
	
	vmg_bmg <- rbind(vmg_diff_topo, bmg_diff_topo)
	vmg_non <- rbind(vmg_diff_topo, non_diff_topo)
	bmg_non <- rbind(bmg_diff_topo, non_diff_topo)	
	vmg_non_edge <- wilcox.test(EdgeCount ~ Dataset, data=vmg_non, paired=FALSE)	
	vmg_non_in <- wilcox.test(Indegree ~ Dataset, data=vmg_non, paired=FALSE)	
	vmg_non_out <- wilcox.test(Outdegree ~ Dataset, data=vmg_non, paired=FALSE)
	vmg_non_path <- wilcox.test(AverageShortestPathLength ~ Dataset, data=vmg_non, paired=FALSE)
	vmg_non_bc <- wilcox.test(BetweennessCentrality ~ Dataset, data=vmg_non, paired=FALSE)
	vmg_non_cc <- wilcox.test(ClosenessCentrality ~ Dataset, data=vmg_non, paired=FALSE)
	vmg_non_clust <- wilcox.test(ClusteringCoefficient ~ Dataset, data=vmg_non, paired=FALSE)
	vmg_non_ecc <- wilcox.test(Eccentricity ~ Dataset, data=vmg_non, paired=FALSE)
	vmg_non_nc <- wilcox.test(NeighborhoodConnectivity ~ Dataset, data=vmg_non, paired=FALSE)
	
	bmg_non_edge <- wilcox.test(EdgeCount ~ Dataset, data=bmg_non, paired=FALSE)
	bmg_non_in <- wilcox.test(Indegree ~ Dataset, data=bmg_non, paired=FALSE)
	bmg_non_out <- wilcox.test(Outdegree ~ Dataset, data=bmg_non, paired=FALSE)
	bmg_non_path <- wilcox.test(AverageShortestPathLength ~ Dataset, data=bmg_non, paired=FALSE)
	bmg_non_bc <- wilcox.test(BetweennessCentrality ~ Dataset, data=bmg_non, paired=FALSE)
	bmg_non_cc <- wilcox.test(ClosenessCentrality ~ Dataset, data=bmg_non, paired=FALSE)
	bmg_non_clust <- wilcox.test(ClusteringCoefficient ~ Dataset, data=bmg_non, paired=FALSE)
	bmg_non_ecc <- wilcox.test(Eccentricity ~ Dataset, data=bmg_non, paired=FALSE)
	bmg_non_nc <- wilcox.test(NeighborhoodConnectivity ~ Dataset, data=bmg_non, paired=FALSE)
	
	vmg_bmg_edge <- wilcox.test(EdgeCount ~ Dataset, data=vmg_bmg, paired=FALSE)
	vmg_bmg_in <- wilcox.test(Indegree ~ Dataset, data=vmg_bmg, paired=FALSE)
	vmg_bmg_out <- wilcox.test(Outdegree ~ Dataset, data=vmg_bmg, paired=FALSE)
	vmg_bmg_path <- wilcox.test(AverageShortestPathLength ~ Dataset, data=vmg_bmg, paired=FALSE)
	vmg_bmg_bc <- wilcox.test(BetweennessCentrality ~ Dataset, data=vmg_bmg, paired=FALSE)
	vmg_bmg_cc <- wilcox.test(ClosenessCentrality ~ Dataset, data=vmg_bmg, paired=FALSE)
	vmg_bmg_clust <- wilcox.test(ClusteringCoefficient ~ Dataset, data=vmg_bmg, paired=FALSE)
	vmg_bmg_ecc <- wilcox.test(Eccentricity ~ Dataset, data=vmg_bmg, paired=FALSE)
	vmg_bmg_nc <- wilcox.test(NeighborhoodConnectivity ~ Dataset, data=vmg_bmg, paired=FALSE)

Plot these results:

	install.packages("ggplot2")
	library(ggplot2)
	
	vmg_bmg_non <- rbind(vmg_bmg, non_diff_topo)
	vmg_bmg_non$Dataset <- factor(vmg_bmg_non$Dataset,levels=c('Non','VMG','BMG'))

Total Degree:

	VMG <- boxplot.stats(vmg_diff_topo$EdgeCount)$stats
	BMG <- boxplot.stats(bmg_diff_topo$EdgeCount)$stats
	Non <- boxplot.stats(non_diff_topo$EdgeCount)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$EdgeCount)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))	#}
	
	plot_degree <- ggplot(vmg_bmg_non, aes(x=Dataset, y=EdgeCount)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Total Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=89)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_degree, file="network_plots/degree_grey.pdf", w=5, h=5)

In-Degree:

	VMG <- boxplot.stats(vmg_diff_topo$Indegree)$stats
	BMG <- boxplot.stats(bmg_diff_topo$Indegree)$stats
	Non <- boxplot.stats(non_diff_topo$Indegree)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$Indegree)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))	#}
	
	plot_indegree <- ggplot(vmg_bmg_non, aes(x=Dataset, y=Indegree)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "In-Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=89)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_indegree, file="network_plots/indegree_grey.pdf", w=5, h=5)

Out-Degree:

	VMG <- boxplot.stats(vmg_diff_topo$Outdegree)$stats
	BMG <- boxplot.stats(bmg_diff_topo$Outdegree)$stats
	Non <- boxplot.stats(non_diff_topo$Outdegree)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("?","?","?")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$Outdegree)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))	#}
	
	plot_outdegree <- ggplot(vmg_bmg_non, aes(x=Dataset, y=Outdegree)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Out-Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=89)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_outdegree, file="network_plots/outdegree_grey.pdf", w=5, h=5)

Average Path:

	VMG <- boxplot.stats(vmg_diff_topo$AverageShortestPathLength)$stats
	BMG <- boxplot.stats(bmg_diff_topo$AverageShortestPathLength)$stats
	Non <- boxplot.stats(non_diff_topo$AverageShortestPathLength)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$AverageShortestPathLength)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))	#}
	
	plot_avgpath <- ggplot(vmg_bmg_non, aes(x=Dataset, y=AverageShortestPathLength)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Average Shortest Path Length",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=89)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_avgpath, file="network_plots/avgpath_grey.pdf", w=5, h=5)

Betweenness Centrality:

	VMG <- boxplot.stats(vmg_diff_topo$BetweennessCentrality)$stats
	BMG <- boxplot.stats(bmg_diff_topo$BetweennessCentrality)$stats
	Non <- boxplot.stats(non_diff_topo$BetweennessCentrality)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$BetweennessCentrality)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))	#}
	
	plot_bc <- ggplot(vmg_bmg_non, aes(x=Dataset, y=BetweennessCentrality)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Betweenness Centrality",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=89)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_bc, file="network_plots/bc_grey.pdf", w=5, h=5)





##Differential KOs mapped to PWs

##Community Stats on all KEGG functions, ORF abundance, metabolic genes


##Viral-Bacterial metagenome read-sharing


##16S alpha and beta-diverstiy

##Illumina Viral Data QC
Different QC steps (minor, mostly in how we trim off adaptors and dealing with ambiguous bases/homopolymers that were dealt with by the torrent server previously) were used for the illumina viral metagenome data. Once again, there apperas to be duplication issues associated with the transposon preps presumambly, so we must be careful in dealing with them to ensure their removal.

Download the raw illumina metagenome reads:

	scp .... or wget

Instead of trimmomatic, for Illumina data I have had better luck with cutadapt (version 1.8.1). Install cutadapt:

	cd virtualenv-*; python virtualenv.py ../cutadapt; cd ..
	source cutadapt/bin/activate
	pip install cutadapt

Trim the adaptors and quality trim:

	source cutadapt/bin/activate
	mkdir cutadapt_qc
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V4_S1_L001_R1_001.trim.fastq.gz raw_illumina_viral_run1/V4_S1_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V9_S2_L001_R1_001.trim.fastq.gz raw_illumina_viral_run1/V9_S2_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V10_S3_L001_R1_001.trim.fastq.gz raw_illumina_viral_run1/V10_S3_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V11_S4_L001_R1_001.trim.fastq.gz raw_illumina_viral_run1/V11_S4_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V12_S5_L001_R1_001.trim.fastq.gz raw_illumina_viral_run1/V12_S5_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V13_S6_L001_R1_001.trim.fastq.gz raw_illumina_viral_run1/V13_S6_L001_R1_001.fastq.gz

	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V4_S1_L001_R2_001.trim.fastq.gz raw_illumina_viral_run1/V4_S1_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V9_S2_L001_R2_001.trim.fastq.gz raw_illumina_viral_run1/V9_S2_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V10_S3_L001_R2_001.trim.fastq.gz raw_illumina_viral_run1/V10_S3_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V11_S4_L001_R2_001.trim.fastq.gz raw_illumina_viral_run1/V11_S4_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V12_S5_L001_R2_001.trim.fastq.gz raw_illumina_viral_run1/V12_S5_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc/V13_S6_L001_R2_001.trim.fastq.gz raw_illumina_viral_run1/V13_S6_L001_R2_001.fastq.gz

	
	mkdir cutadapt_qc2
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V4_S1_L001_R1_001.trim.fastq.gz raw_illumina_viral_run2/V4_S1_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V9_S2_L001_R1_001.trim.fastq.gz raw_illumina_viral_run2/V9_S2_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V10_S3_L001_R1_001.trim.fastq.gz raw_illumina_viral_run2/V10_S3_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V11_S4_L001_R1_001.trim.fastq.gz raw_illumina_viral_run2/V11_S4_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V12_S5_L001_R1_001.trim.fastq.gz raw_illumina_viral_run2/V12_S5_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V13_S6_L001_R1_001.trim.fastq.gz raw_illumina_viral_run2/V13_S6_L001_R1_001.fastq.gz

	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V4_S1_L001_R2_001.trim.fastq.gz raw_illumina_viral_run2/V4_S1_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V9_S2_L001_R2_001.trim.fastq.gz raw_illumina_viral_run2/V9_S2_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V10_S3_L001_R2_001.trim.fastq.gz raw_illumina_viral_run2/V10_S3_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V11_S4_L001_R2_001.trim.fastq.gz raw_illumina_viral_run2/V11_S4_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V12_S5_L001_R2_001.trim.fastq.gz raw_illumina_viral_run2/V12_S5_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc2/V13_S6_L001_R2_001.trim.fastq.gz raw_illumina_viral_run2/V13_S6_L001_R2_001.fastq.gz

	
Prefer userach quality filter based on expected errors per read for filtering illumina data.  Need a free liscence to download it.  Go to http://www.drive5.com/usearch/download.html and download the linux version of USEARCH v8.0.1623 and enter your email address.  Copy the download link in the email address and use below:

	wget <download_link>
	mv <name> usearch8.0.1623
	chmod 775 usearch8.0.1623 

Remove seqeunces that have an estimated error rate >1%:
	
	mkdir usearch_qc
	gunzip -d cutadapt_qc/*
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V4_S1_L001_R1_001.trim.fastq -fastqout usearch_qc/V4_S1_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V4_S1_L001_R2_001.trim.fastq -fastqout usearch_qc/V4_S1_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V9_S2_L001_R1_001.trim.fastq -fastqout usearch_qc/V9_S2_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V9_S2_L001_R2_001.trim.fastq -fastqout usearch_qc/V9_S2_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V10_S3_L001_R1_001.trim.fastq -fastqout usearch_qc/V10_S3_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V10_S3_L001_R2_001.trim.fastq -fastqout usearch_qc/V10_S3_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V11_S4_L001_R1_001.trim.fastq -fastqout usearch_qc/V11_S4_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V11_S4_L001_R2_001.trim.fastq -fastqout usearch_qc/V11_S4_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V12_S5_L001_R1_001.trim.fastq -fastqout usearch_qc/V12_S5_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V12_S5_L001_R2_001.trim.fastq -fastqout usearch_qc/V12_S5_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V13_S6_L001_R1_001.trim.fastq -fastqout usearch_qc/V13_S6_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V13_S6_L001_R2_001.trim.fastq -fastqout usearch_qc/V13_S6_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01

	mkdir usearch_qc2
	gunzip -d cutadapt_qc2/*
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V4_S1_L001_R1_001.trim.fastq -fastqout usearch_qc2/V4_S1_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V4_S1_L001_R2_001.trim.fastq -fastqout usearch_qc2/V4_S1_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V9_S2_L001_R1_001.trim.fastq -fastqout usearch_qc2/V9_S2_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V9_S2_L001_R2_001.trim.fastq -fastqout usearch_qc2/V9_S2_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V10_S3_L001_R1_001.trim.fastq -fastqout usearch_qc2/V10_S3_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V10_S3_L001_R2_001.trim.fastq -fastqout usearch_qc2/V10_S3_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V11_S4_L001_R1_001.trim.fastq -fastqout usearch_qc2/V11_S4_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V11_S4_L001_R2_001.trim.fastq -fastqout usearch_qc2/V11_S4_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V12_S5_L001_R1_001.trim.fastq -fastqout usearch_qc2/V12_S5_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V12_S5_L001_R2_001.trim.fastq -fastqout usearch_qc2/V12_S5_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V13_S6_L001_R1_001.trim.fastq -fastqout usearch_qc2/V13_S6_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V13_S6_L001_R2_001.trim.fastq -fastqout usearch_qc2/V13_S6_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01

	
Removed duplicates:

	./usearch8.0.1623 -derep_prefix usearch_qc/V4_S1_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc/V4.R1.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc/V9_S2_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc/V9.R1.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc/V10_S3_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc/V10.R1.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc/V11_S4_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc/V11.R1.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc/V12_S5_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc/V12.R1.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc/V13_S6_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc/V13.R1.unique.fasta -minseqlength 80

	./usearch8.0.1623 -derep_prefix usearch_qc/V4_S1_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc/V4.R2.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc/V9_S2_L001_R2_001_usearch_error.fastq  -fastaout usearch_qc/V9.R2.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc/V10_S3_L001_R2_001_usearch_error.fastq  -fastaout usearch_qc/V10.R2.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc/V11_S4_L001_R2_001_usearch_error.fastq  -fastaout usearch_qc/V11.R2.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc/V12_S5_L001_R2_001_usearch_error.fastq  -fastaout usearch_qc/V12.R2.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc/V13_S6_L001_R2_001_usearch_error.fastq  -fastaout usearch_qc/V13.R2.unique.fasta -minseqlength 80


	./usearch8.0.1623 -derep_prefix usearch_qc2/V4_S1_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc2/V4.R1.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc2/V9_S2_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc2/V9.R1.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc2/V10_S3_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc2/V10.R1.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc2/V11_S4_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc2/V11.R1.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc2/V12_S5_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc2/V12.R1.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc2/V13_S6_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc2/V13.R1.unique.fasta -minseqlength 80

	./usearch8.0.1623 -derep_prefix usearch_qc2/V4_S1_L001_R1_001_usearch_error.fastq  -fastaout usearch_qc2/V4.R2.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc2/V9_S2_L001_R2_001_usearch_error.fastq  -fastaout usearch_qc2/V9.R2.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc2/V10_S3_L001_R2_001_usearch_error.fastq  -fastaout usearch_qc2/V10.R2.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc2/V11_S4_L001_R2_001_usearch_error.fastq  -fastaout usearch_qc2/V11.R2.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc2/V12_S5_L001_R2_001_usearch_error.fastq  -fastaout usearch_qc2/V12.R2.unique.fasta -minseqlength 80
	./usearch8.0.1623 -derep_prefix usearch_qc2/V13_S6_L001_R2_001_usearch_error.fastq  -fastaout usearch_qc2/V13.R2.unique.fasta -minseqlength 80
	


Now, we can use QIIME to filter the original FASTQ files with the seqeunce identifiers from the truncated/dereplicated FASTQ files:

    source qiimeEnv/bin/activate
    source qiimeEnv/bin/activate
	filter_fasta.py -f usearch_qc/V4_S1_L001_R1_001_usearch_error.fastq -s usearch_qc/V4_S1_L001_R1_001_trunc75_ids.txt -o usearch_qc/V4_S1_L001_R1_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc/V9_S2_L001_R1_001_usearch_error.fastq -s usearch_qc/V9_S2_L001_R1_001_trunc75_ids.txt -o usearch_qc/V9_S2_L001_R1_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc/V10_S3_L001_R1_001_usearch_error.fastq -s usearch_qc/V10_S3_L001_R1_001_trunc75_ids.txt -o usearch_qc/V10_S3_L001_R1_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc/V11_S4_L001_R1_001_usearch_error.fastq -s usearch_qc/V11_S4_L001_R1_001_trunc75_ids.txt -o usearch_qc/V11_S4_L001_R1_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc/V12_S5_L001_R1_001_usearch_error.fastq -s usearch_qc/V12_S5_L001_R1_001_trunc75_ids.txt -o usearch_qc/V12_S5_L001_R1_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc/V13_S6_L001_R1_001_usearch_error.fastq -s usearch_qc/V13_S6_L001_R1_001_trunc75_ids.txt -o usearch_qc/V13_S6_L001_R1_001_usearch_finalqc.fastq

	filter_fasta.py -f usearch_qc/V4_S1_L001_R2_001_usearch_error.fastq -s usearch_qc/V4_S1_L001_R2_001_trunc75_ids.txt -o usearch_qc/V4_S1_L001_R2_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc/V9_S2_L001_R2_001_usearch_error.fastq -s usearch_qc/V9_S2_L001_R2_001_trunc75_ids.txt -o usearch_qc/V9_S2_L001_R2_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc/V10_S3_L001_R2_001_usearch_error.fastq -s usearch_qc/V10_S3_L001_R2_001_trunc75_ids.txt -o usearch_qc/V10_S3_L001_R2_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc/V11_S4_L001_R2_001_usearch_error.fastq -s usearch_qc/V11_S4_L001_R2_001_trunc75_ids.txt -o usearch_qc/V11_S4_L001_R2_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc/V12_S5_L001_R2_001_usearch_error.fastq -s usearch_qc/V12_S5_L001_R2_001_trunc75_ids.txt -o usearch_qc/V12_S5_L001_R2_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc/V13_S6_L001_R2_001_usearch_error.fastq -s usearch_qc/V13_S6_L001_R2_001_trunc75_ids.txt -o usearch_qc/V13_S6_L001_R2_001_usearch_finalqc.fastq

	filter_fasta.py -f usearch_qc2/V4_S1_L001_R1_001_usearch_error.fastq -s usearch_qc2/V4_S1_L001_R1_001_trunc75_ids.txt -o usearch_qc2/V4_S1_L001_R1_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc2/V9_S2_L001_R1_001_usearch_error.fastq -s usearch_qc2/V9_S2_L001_R1_001_trunc75_ids.txt -o usearch_qc2/V9_S2_L001_R1_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc2/V10_S3_L001_R1_001_usearch_error.fastq -s usearch_qc2/V10_S3_L001_R1_001_trunc75_ids.txt -o usearch_qc2/V10_S3_L001_R1_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc2/V11_S4_L001_R1_001_usearch_error.fastq -s usearch_qc2/V11_S4_L001_R1_001_trunc75_ids.txt -o usearch_qc2/V11_S4_L001_R1_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc2/V12_S5_L001_R1_001_usearch_error.fastq -s usearch_qc2/V12_S5_L001_R1_001_trunc75_ids.txt -o usearch_qc2/V12_S5_L001_R1_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc2/V13_S6_L001_R1_001_usearch_error.fastq -s usearch_qc2/V13_S6_L001_R1_001_trunc75_ids.txt -o usearch_qc2/V13_S6_L001_R1_001_usearch_finalqc.fastq

	filter_fasta.py -f usearch_qc2/V4_S1_L001_R2_001_usearch_error.fastq -s usearch_qc2/V4_S1_L001_R2_001_trunc75_ids.txt -o usearch_qc2/V4_S1_L001_R2_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc2/V9_S2_L001_R2_001_usearch_error.fastq -s usearch_qc2/V9_S2_L001_R2_001_trunc75_ids.txt -o usearch_qc2/V9_S2_L001_R2_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc2/V10_S3_L001_R2_001_usearch_error.fastq -s usearch_qc2/V10_S3_L001_R2_001_trunc75_ids.txt -o usearch_qc2/V10_S3_L001_R2_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc2/V11_S4_L001_R2_001_usearch_error.fastq -s usearch_qc2/V11_S4_L001_R2_001_trunc75_ids.txt -o usearch_qc2/V11_S4_L001_R2_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc2/V12_S5_L001_R2_001_usearch_error.fastq -s usearch_qc2/V12_S5_L001_R2_001_trunc75_ids.txt -o usearch_qc2/V12_S5_L001_R2_001_usearch_finalqc.fastq
	filter_fasta.py -f usearch_qc2/V13_S6_L001_R2_001_usearch_error.fastq -s usearch_qc2/V13_S6_L001_R2_001_trunc75_ids.txt -o usearch_qc2/V13_S6_L001_R2_001_usearch_finalqc.fastq

Added one more step to ensure removal of duplicates.  Also, combined read1 and read2 together because we had many orphans as it was due to quality falling quickly on v3 MiSeq chemsitry on read2:

	mkdir prinseq_illumina
    cd usearch_qc
    for f in *_R2_finalqc.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        cd ~
        perl prinseq-lite-0.20.4/./prinseq-lite.pl -derep 14 -lc_method dust -lc_threshold 7 -fastq usearch_qc/$f -out_format 3 -out_good prinseq_illumina/"$filename""_prinseq"
    done
    
    mkdir prinseq_illumina2
    cd usearch_qc2
    for f in *_R2_finalqc.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        cd ~
        perl prinseq-lite-0.20.4/./prinseq-lite.pl -derep 14 -lc_method dust -lc_threshold 7 -fastq usearch_qc2/$f -out_format 3 -out_good prinseq_illumina2/"$filename""_prinseq"
    done


Combine reads from all runs together:

   
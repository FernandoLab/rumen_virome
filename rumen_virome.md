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
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.1.fastq trimmomatic_output/VMG.1_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.2.fastq trimmomatic_output/VMG.2_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.3.fastq trimmomatic_output/VMG.3_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.4.fastq trimmomatic_output/VMG.4_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.5.fastq trimmomatic_output/VMG.5_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:75 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.6.fastq trimmomatic_output/VMG.6_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.7.fastq trimmomatic_output/VMG.7_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.8.fastq trimmomatic_output/VMG.8_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.9.fastq trimmomatic_output/VMG.9_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.10.fastq trimmomatic_output/VMG.10_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.11.fastq trimmomatic_output/VMG.11_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.12.fastq trimmomatic_output/VMG.12_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.13.fastq trimmomatic_output/VMG.13_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.14.fastq trimmomatic_output/VMG.14_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.15.fastq trimmomatic_output/VMG.15_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:15 CROP:150 MINLEN:75

Alternatively, the outputs from the above command are available at interm/trimmomatic_output.

To remove artifical duplicates we use cd-hit-454.  We noticed early on the transposon library preps with the Ion Torrent Museek kits introduces (unfortunately sometimes many) duplicates (associated with viral sample?). While the software was originally designed for 454 sequencing, Ion Torrent shares many of the same error profiles as 454 and cd-hit-454 appears to work fairly well for removing duplicates in our data.

First download cd-hit-454:

    wget https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz
    tar xvf cd-hit-v4.6.1-2012-08-27.tgz
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


Alternatively, the outputs from the above command are available at interm/cd_hit_454_output.


Now, we wanted to double check once again for other types of duplicates and some that may have been missed by cd-hit-454. Further, we can use prinseq to remove transposon associated seqeunces that kept showing up in the 3' end.  We dtermined that trimming the last 5% of a read did the trick. 


Download prinseq-lite:

    wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz
    tar xvf prinseq-lite-0.20.4.tar.gz

Run prinseq on each viral metagenome:

    mkdir prinseq_output
    cd cd_hit_454_output
    for f in *.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        cd ~
        perl prinseq-lite-0.20.4/./prinseq-lite.pl -custom_params "TGAACTG 1;GAACTGA 1;AACTGAC 1;ACTGACG 1;CTGACGC 1;TGACGCA 1;GACGCAC 1;ACGCACG 1;CGCACGA 1;GCACGAA 1;" -derep 12345 -lc_method dust -lc_threshold 7 -fastq cd_hit_454_output/$f -out_format 3 -out_good prinseq_output/"$filename""_prinseq"
    done


Alternative, the outputs for these commands are in interm/prinseq_outputs with extensions prinseq.fastq.

To apply one additional filter to remove artifcats/duplicates, we used 60 bp prefix duplicate filtering.  This can be done again using prinseq and then removing the resulting sequences from the previous file for downstream analysis. We used prinseq again and truncated seqeunces to 60 bp and used exact and exact reverse complement dereplicating:

    cd prinseq_output
    for f in *.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        perl ~/prinseq-lite-0.20.4/prinseq-lite.pl -trim_to_len 60 -derep 14 -fastq $f -out_format 3 -out_good "$filename""_truncatederep"
    done

Alternatively, the outputs are located at:?

Now, we can use the generated truncate/dereplicated FASTQ files to pull teh seqeunce records out of the original FASTQ file to keep for downstream analysis.  QIIME has a command to do this quickly, but downloading it can be a pain.  However, we will use it downstream again later, so here is how I downloaded it into a linux system. If you already have access to QIIME I would just run the commands below, but in the effort to recreaste the analysis, this is how I would install QIIME on a linux (and likely works on mac?) system by creating a virtual enviornment:

    wget https://pypi.python.org/packages/source/v/virtualenv/virtualenv-12.0.7.tar.gz#md5=e08796f79d112f3bfa6653cc10840114
    tar xzf virtualenv-12.0.7.tar.gz
    cd virtualenv-*; python virtualenv.py ../qiimeEnv; cd ..
    source qiimeEnv/bin/activate
    pip install numpy
    pip install qiime


Now, we can use QIIME to filter the original FASTQ files with the seqeunce identifiers from the truncated/dereplicated FASTQ files:

    source qiimeEnv/bin/activate
    cd prinseq_output
    for f in *_truncatederep.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        convert_fastaqual_fastq.py -f $f -c fastq_to_fastaqual
        grep ">" "$filename"_truncatederep.fna | cut -c 2- > "$filename"_keep_ids.txt
        filter_fasta.py -f "$filename"_prinseq.fastq.fastq -s "$filename"_keep_ids.txt -o "$filename"_finalQC.fastq
    done

Alternatively, the final QC files generated are available at interm/prinseq_output with the extensions finalQC.fastq

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
	for f in *_finalQC.fastq
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
	cp prinseq_output/*_finalQC.fna rRNA_prediction/examples/e1/input_viral

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

Take out later:
perl_scripts/./length_distribution.pl vmg.contigs.fasta 
Total reads:	199788
Total nt:	115856508
Mean length:	579.897231064929
Median length:	327
Mode length:	228
Max length:	56361
Min length:	31
Length	# Seqs
10	0
35	49
50	236
75	331
99	246
100	9
150	301
200	335
250	34413
300	48936
350	25638
400	17492
450	11863
500	8825
550	6801
600	5401
650	4179
700	3524
750	2943
800	2481
850	2215
900	1854
950	1580
1000	1421



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
	
	source qiimeEnv/bin/activate
	filter_fasta.py -f vmg.orfs_aa.fasta -a vmg.orfs_nt.filter.fasta -o vmg.orfs_aa.filter.fasta
	



Total:

	
	




	


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


Search Total against PHAST:



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


Align all viral reads to nt ORFs:

    for f in prinseq_output/*_finalQC.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        bowtie2-2.2.5/bowtie2 -U $f --end-to-end --sensitive -x vmg_orfs_bowtie/vmg_orfs_bowtie_db -S vmg_orfs_bowtie/$filename.psl --un vmg_orfs_bowtie/$filename.unaligned.txt --al vmg_orfs_bowtie/$filename.aligned.txt 
    done
 
All all total metagenome reads to nt ORFs:
	



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




##Create Metabolic Networks
We used mmnet package in R to create two metabolic netwroks: one with all the KOs found in viral and total metagenome and another that represents all of global metabolism as defined by KEGG.




##Viral-Bacterial metagenome read-sharing


##16S alpha and beta-diverstiy


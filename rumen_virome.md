Impact of dietary change on auxillary viral metabolism and viral-bacterial interactions in the rumen of beef cattle
===============
Author: Christopher L. Anderson (canderson30@unl.edu)


##Introduction
To recreate the analysis from Anderson et al. manuscript.  All commands below were done in a linux enviornment and memory intensive commands were carried out on HCC resources at UNL.


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
    for f in trimmomatic_output/*.fastq
    do
        filename=$(basename "$f")
        filename="${filename%.*}"
        cd ~
        cd-hit-v4.6.1-2012-08-27/./cd-hit-454 -i $f -o cd_hit_454_output/"$filename""_cd454.fastq" -M 6100 -T 10
    done


Now, we wanted to check for other types of duplicates and some that may have been missed by cd-hit-454. Further, we can use prinseq to remove transposon associated seqeunces that kept showing up in the 3' end.

Download prinseq-lite:

    wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz
    tar -xvf prinseq-lite-0.20.4.tar.gz

Run prinseq on each viral metagenome:

    mkdir prinseq_output
    for f in cd_hit_454_output/*.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        cd ~
        perl prinseq-lite-0.20.4/./prinseq-lite.pl -custom_params "TGAACTG 1;GAACTGA 1;AACTGAC 1;ACTGACG 1;CTGACGC 1;TGACGCA 1;GACGCAC 1;ACGCACG 1;CGCACGA 1;GCACGAA 1;" -derep 14 -lc_method dust -lc_threshold 7 -fastq $f -out_format 3 -out_good prinseq_output/"$filename""_prinseq"
    done

Use prefix length to remove duplicates with exact matches over first 75 basepairs:
	
    for f in prinseq_output/*.fastq
	do
		filename=$(basename "$f")
    	filename="${filename%_*}"
    	perl prinseq-lite-0.20.4/prinseq-lite.pl -trim_to_len 75 -derep 1 -fastq $f -out_format 2 -out_good prinseq_output/"$filename""_truncatederep"
	done


Use qiime to remove those identified as prefix duplicates. Easiest way to download QIIME I have found is with anaconda package manager.
Here, we download QIIME along with some other programs we need later:
	
	wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda-2.3.0-Linux-x86_64.sh
	bash Anaconda-2.3.0-Linux-x86_64.sh
	anaconda/bin/conda create -n rumenVirome python pip numpy matplotlib scipy pandas cython mock nose
	source anaconda/bin/activate rumenVirome
	anaconda/bin/conda install -c https://conda.binstar.org/jorge qiime
	anaconda/bin/conda install -c r r
	pip install Cython --upgrade
	pip install cutadapt
	pip install khmer
	
	source anaconda/bin/activate rumenVirome
	cd prinseq_output
	for f in *truncatederep.fasta
	do
   		filename=$(basename "$f")
    	filename="${filename%_*}"
    	grep ">" $f | cut -c 2- > "$filename"_keep_ids.txt
    	filter_fasta.py -f "$filename"_prinseq.fastq -s "$filename"_keep_ids.txt -o "$filename"_finalQC.fastq
		convert_fastaqual_fastq.py -f "$filename"_finalQC.fastq -c fastq_to_fastaqual
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
    for f in trimmomatic_output_total/*.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        cd-hit-v4.6.1-2012-08-27/./cd-hit-454 -i $f -o cd_hit_454_output_total/"$filename""_cd454.fastq" -M 6100 -T 10
    done

##Illumina Viral Data QC
Different QC steps (minor, mostly in how we trim off adaptors and dealing with ambiguous bases/homopolymers that were dealt with by the torrent server previously) were used for the illumina viral metagenome data. Once again, there apperas to be duplication issues associated with the transposon preps presumambly, so we must be careful in dealing with them to ensure their removal.

Download the raw illumina metagenome reads, put on our public server for now:

	scp .... or wget
	
Instead of trimmomatic, for Illumina data I have had better luck with cutadapt (version 1.8.1).
Trim the adaptors and quality trim:

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
	
	mkdir cutadapt_qc3
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc3/v9_S2_L001_R1_001.trim.fastq.gz raw_illumina_viral_run3/v9_S2_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc3/v10_S3_L001_R1_001.trim.fastq.gz raw_illumina_viral_run3/v10_S3_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc3/v12_S4_L001_R1_001.trim.fastq.gz raw_illumina_viral_run3/v12_S4_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc3/v13_S5_L001_R1_001.trim.fastq.gz raw_illumina_viral_run3/v13_S5_L001_R1_001.fastq.gz

	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc3/v9_S2_L001_R2_001.trim.fastq.gz raw_illumina_viral_run3/v9_S2_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc3/v10_S3_L001_R2_001.trim.fastq.gz raw_illumina_viral_run3/v10_S3_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc3/v12_S4_L001_R2_001.trim.fastq.gz raw_illumina_viral_run3/v12_S4_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc3/v13_S5_L001_R2_001.trim.fastq.gz raw_illumina_viral_run3/v13_S5_L001_R2_001.fastq.gz

	mkdir cutadapt_qc4
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc4/v4_S1_L001_R1_001.trim.fastq.gz raw_illumina_viral_run4/v4_S1_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc4/v10_S2_L001_R1_001.trim.fastq.gz raw_illumina_viral_run4/v10_S2_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc4/v11_S3_L001_R1_001.trim.fastq.gz raw_illumina_viral_run4/v11_S3_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc4/v12_S4_L001_R1_001.trim.fastq.gz raw_illumina_viral_run4/v12_S4_L001_R1_001.fastq.gz
	cutadapt -n 2 -u -25 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc4/v13_S5_L001_R1_001.trim.fastq.gz raw_illumina_viral_run4/v13_S5_L001_R1_001.fastq.gz

	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc4/v4_S1_L001_R2_001.trim.fastq.gz raw_illumina_viral_run4/v4_S1_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc4/v10_S2_L001_R2_001.trim.fastq.gz raw_illumina_viral_run4/v10_S2_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc4/v11_S3_L001_R2_001.trim.fastq.gz raw_illumina_viral_run4/v11_S3_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc4/v12_S4_L001_R2_001.trim.fastq.gz raw_illumina_viral_run4/v12_S4_L001_R2_001.fastq.gz
	cutadapt -n 2 -u -100 -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTGACGCTGCCGACGA -b GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -b AGATGTGTATAAGAGACAG -b CTGTCTCTTATACACATCT -o cutadapt_qc4/v13_S5_L001_R2_001.trim.fastq.gz raw_illumina_viral_run4/v13_S5_L001_R2_001.fastq.gz


Prefer userach quality filter based on expected errors per read for filtering illumina data.  Need a free liscence to download it.  Go to http://www.drive5.com/usearch/download.html and download the linux version of USEARCH v8.0.1623 and enter your email address.  Copy the download link in the email address and use below:

	wget <download_link> -o ./usearch8.0.1623
	chmod 775 ./usearch8.0.1623 
	mv ./usearch8.0.1623 anaconda/envs/rumenVirome/bin/
	
Remove seqeunces that have an estimated error rate >1%:
	
	mkdir usearch_qc
	gunzip -d cutadapt_qc/*
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V4_S1_L001_R1_001.trim.fastq -fastqout usearch_qc/V4_S1_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V4_S1_L001_R2_001.trim.fastq -fastqout usearch_qc/V4_S1_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V9_S2_L001_R1_001.trim.fastq -fastqout usearch_qc/V9_S2_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V9_S2_L001_R2_001.trim.fastq -fastqout usearch_qc/V9_S2_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V10_S3_L001_R1_001.trim.fastq -fastqout usearch_qc/V10_S3_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V10_S3_L001_R2_001.trim.fastq -fastqout usearch_qc/V10_S3_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V11_S4_L001_R1_001.trim.fastq -fastqout usearch_qc/V11_S4_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V11_S4_L001_R2_001.trim.fastq -fastqout usearch_qc/V11_S4_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V12_S5_L001_R1_001.trim.fastq -fastqout usearch_qc/V12_S5_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V12_S5_L001_R2_001.trim.fastq -fastqout usearch_qc/V12_S5_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V13_S6_L001_R1_001.trim.fastq -fastqout usearch_qc/V13_S6_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc/V13_S6_L001_R2_001.trim.fastq -fastqout usearch_qc/V13_S6_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80

	mkdir usearch_qc2
	gunzip -d cutadapt_qc2/*
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V4_S1_L001_R1_001.trim.fastq -fastqout usearch_qc2/V4_S1_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V4_S1_L001_R2_001.trim.fastq -fastqout usearch_qc2/V4_S1_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V9_S2_L001_R1_001.trim.fastq -fastqout usearch_qc2/V9_S2_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V9_S2_L001_R2_001.trim.fastq -fastqout usearch_qc2/V9_S2_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V10_S3_L001_R1_001.trim.fastq -fastqout usearch_qc2/V10_S3_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V10_S3_L001_R2_001.trim.fastq -fastqout usearch_qc2/V10_S3_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V11_S4_L001_R1_001.trim.fastq -fastqout usearch_qc2/V11_S4_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V11_S4_L001_R2_001.trim.fastq -fastqout usearch_qc2/V11_S4_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V12_S5_L001_R1_001.trim.fastq -fastqout usearch_qc2/V12_S5_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V12_S5_L001_R2_001.trim.fastq -fastqout usearch_qc2/V12_S5_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V13_S6_L001_R1_001.trim.fastq -fastqout usearch_qc2/V13_S6_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc2/V13_S6_L001_R2_001.trim.fastq -fastqout usearch_qc2/V13_S6_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	
	mkdir usearch_qc3
	gunzip -d cutadapt_qc3/*
	./usearch8.0.1623 -fastq_filter cutadapt_qc3/v9_S2_L001_R1_001.trim.fastq -fastqout usearch_qc3/v9_S2_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc3/v9_S2_L001_R2_001.trim.fastq -fastqout usearch_qc3/v9_S2_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc3/v10_S3_L001_R1_001.trim.fastq -fastqout usearch_qc3/v10_S3_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc3/v10_S3_L001_R2_001.trim.fastq -fastqout usearch_qc3/v10_S3_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc3/v12_S4_L001_R1_001.trim.fastq -fastqout usearch_qc3/v12_S4_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc3/v12_S4_L001_R2_001.trim.fastq -fastqout usearch_qc3/v12_S4_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc3/v13_S5_L001_R1_001.trim.fastq -fastqout usearch_qc3/v13_S5_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc3/v13_S5_L001_R2_001.trim.fastq -fastqout usearch_qc3/v13_S5_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80

	mkdir usearch_qc4
	gunzip -d cutadapt_qc4/*
	./usearch8.0.1623 -fastq_filter cutadapt_qc4/v4_S1_L001_R1_001.trim.fastq -fastqout usearch_qc4/v4_S1_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc4/v4_S1_L001_R2_001.trim.fastq -fastqout usearch_qc4/v4_S1_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc4/v10_S2_L001_R1_001.trim.fastq -fastqout usearch_qc4/v10_S2_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc4/v10_S2_L001_R2_001.trim.fastq -fastqout usearch_qc4/v10_S2_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc4/v11_S3_L001_R1_001.trim.fastq -fastqout usearch_qc4/v11_S3_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc4/v11_S3_L001_R2_001.trim.fastq -fastqout usearch_qc4/v11_S3_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc4/v12_S4_L001_R1_001.trim.fastq -fastqout usearch_qc4/v12_S4_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc4/v12_S4_L001_R2_001.trim.fastq -fastqout usearch_qc4/v12_S4_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc4/v13_S5_L001_R1_001.trim.fastq -fastqout usearch_qc4/v13_S5_L001_R1_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80
	./usearch8.0.1623 -fastq_filter cutadapt_qc4/v13_S5_L001_R2_001.trim.fastq -fastqout usearch_qc4/v13_S5_L001_R2_001_usearch_error.fastq -fastq_maxns 1 -fastq_maxee_rate 0.01 -fastq_minlen 80


merge read1 and read2 together, many of reverse reads lost due to quality falling off early in V3 chemistry, not worth getting pairs back together:

	cat usearch_qc/V4_S1_L001_R1_001_usearch_error.fastq usearch_qc/V4_S1_L001_R2_001_usearch_error.fastq > usearch_qc/V4_S1_L001_R1_R2.fastq
	cat usearch_qc/V9_S2_L001_R1_001_usearch_error.fastq usearch_qc/V9_S2_L001_R2_001_usearch_error.fastq > usearch_qc/V9_S2_L001_R1_R2.fastq
	cat usearch_qc/V10_S3_L001_R1_001_usearch_error.fastq usearch_qc/V10_S3_L001_R2_001_usearch_error.fastq > usearch_qc/V10_S3_L001_R1_R2.fastq
	cat usearch_qc/V11_S4_L001_R1_001_usearch_error.fastq usearch_qc/V11_S4_L001_R2_001_usearch_error.fastq > usearch_qc/V11_S4_L001_R1_R2.fastq
	cat usearch_qc/V12_S5_L001_R1_001_usearch_error.fastq usearch_qc/V12_S5_L001_R2_001_usearch_error.fastq > usearch_qc/V12_S5_L001_R1_R2.fastq
	cat usearch_qc/V13_S6_L001_R1_001_usearch_error.fastq usearch_qc/V13_S6_L001_R2_001_usearch_error.fastq > usearch_qc/V13_S6_L001_R1_R2.fastq

	cat usearch_qc2/V4_S1_L001_R1_001_usearch_error.fastq usearch_qc2/V4_S1_L001_R2_001_usearch_error.fastq > usearch_qc2/V4_S1_L001_R1_R2.fastq
	cat usearch_qc2/V9_S2_L001_R1_001_usearch_error.fastq usearch_qc2/V9_S2_L001_R2_001_usearch_error.fastq > usearch_qc2/V9_S2_L001_R1_R2.fastq
	cat usearch_qc2/V10_S3_L001_R1_001_usearch_error.fastq usearch_qc2/V10_S3_L001_R2_001_usearch_error.fastq > usearch_qc2/V10_S3_L001_R1_R2.fastq
	cat usearch_qc2/V11_S4_L001_R1_001_usearch_error.fastq usearch_qc2/V11_S4_L001_R2_001_usearch_error.fastq > usearch_qc2/V11_S4_L001_R1_R2.fastq
	cat usearch_qc2/V12_S5_L001_R1_001_usearch_error.fastq usearch_qc2/V12_S5_L001_R2_001_usearch_error.fastq > usearch_qc2/V12_S5_L001_R1_R2.fastq
	cat usearch_qc2/V13_S6_L001_R1_001_usearch_error.fastq usearch_qc2/V13_S6_L001_R2_001_usearch_error.fastq > usearch_qc2/V13_S6_L001_R1_R2.fastq

	cat usearch_qc3/v9_S2_L001_R1_001_usearch_error.fastq usearch_qc3/v9_S2_L001_R2_001_usearch_error.fastq > usearch_qc3/v9_S2_L001_R1_R2.fastq
	cat usearch_qc3/v10_S3_L001_R1_001_usearch_error.fastq usearch_qc3/v10_S3_L001_R2_001_usearch_error.fastq > usearch_qc3/v10_S3_L001_R1_R2.fastq
	cat usearch_qc3/v12_S4_L001_R1_001_usearch_error.fastq usearch_qc3/v12_S4_L001_R2_001_usearch_error.fastq > usearch_qc3/v12_S4_L001_R1_R2.fastq
	cat usearch_qc3/v13_S5_L001_R1_001_usearch_error.fastq usearch_qc3/v13_S5_L001_R2_001_usearch_error.fastq > usearch_qc3/v13_S5_L001_R1_R2.fastq

	cat usearch_qc4/v4_S1_L001_R1_001_usearch_error.fastq usearch_qc4/v4_S1_L001_R2_001_usearch_error.fastq > usearch_qc4/v4_S1_L001_R1_R2.fastq
	cat usearch_qc4/v10_S2_L001_R1_001_usearch_error.fastq usearch_qc4/v10_S2_L001_R2_001_usearch_error.fastq > usearch_qc4/v10_S2_L001_R1_R2.fastq
	cat usearch_qc4/v11_S3_L001_R1_001_usearch_error.fastq usearch_qc4/v11_S3_L001_R2_001_usearch_error.fastq > usearch_qc4/v11_S3_L001_R1_R2.fastq
	cat usearch_qc4/v12_S4_L001_R1_001_usearch_error.fastq usearch_qc4/v12_S4_L001_R2_001_usearch_error.fastq > usearch_qc4/v12_S4_L001_R1_R2.fastq
	cat usearch_qc4/v13_S5_L001_R1_001_usearch_error.fastq usearch_qc4/v13_S5_L001_R2_001_usearch_error.fastq > usearch_qc4/v13_S5_L001_R1_R2.fastq

Use prinseq to remove duplicates with prefix filter of 80 bp:
	
	for f in usearch_qc/*_R2.fastq
	do
		filename=$(basename "$f")
		filename="${filename%_*}"
		perl prinseq-lite-0.20.4/prinseq-lite.pl -trim_to_len 80 -derep 1 -fastq $f -out_format 2 -out_good usearch_qc/"$filename""_R2_truncatederep"
		grep ">" usearch_qc/"$filename""_R2_truncatederep.fasta" | cut -c 2- > usearch_qc/"$filename""_R2_keep_ids.txt"
		filter_fasta.py -f $f -s usearch_qc/"$filename""_R2_keep_ids.txt" -o usearch_qc/"$filename""_R2_finalQC.fastq"
	done
	
	for f in usearch_qc2/*_R2.fastq
	do
		filename=$(basename "$f")
		filename="${filename%_*}"
		perl prinseq-lite-0.20.4/prinseq-lite.pl -trim_to_len 80 -derep 1 -fastq $f -out_format 2 -out_good usearch_qc2/"$filename""_R2_truncatederep"
		grep ">" usearch_qc2/"$filename""_R2_truncatederep.fasta" | cut -c 2- > usearch_qc2/"$filename""_R2_keep_ids.txt"
		filter_fasta.py -f $f -s usearch_qc2/"$filename""_R2_keep_ids.txt" -o usearch_qc2/"$filename""_R2_finalQC.fastq"
	done
	
	for f in usearch_qc3/*_R2.fastq
	do
		filename=$(basename "$f")
		filename="${filename%_*}"
		perl prinseq-lite-0.20.4/prinseq-lite.pl -trim_to_len 80 -derep 1 -fastq $f -out_format 2 -out_good usearch_qc3/"$filename""_R2_truncatederep"
		grep ">" usearch_qc3/"$filename""_R2_truncatederep.fasta" | cut -c 2- > usearch_qc3/"$filename""_R2_keep_ids.txt"
		filter_fasta.py -f $f -s usearch_qc3/"$filename""_R2_keep_ids.txt" -o usearch_qc3/"$filename""_R2_finalQC.fastq"
	done
	
	for f in usearch_qc4/*_R2.fastq
	do
		filename=$(basename "$f")
		filename="${filename%_*}"
		perl prinseq-lite-0.20.4/prinseq-lite.pl -trim_to_len 80 -derep 1 -fastq $f -out_format 2 -out_good usearch_qc4/"$filename""_R2_truncatederep"
		grep ">" usearch_qc4/"$filename""_R2_truncatederep.fasta" | cut -c 2- > usearch_qc4/"$filename""_R2_keep_ids.txt"
		filter_fasta.py -f $f -s usearch_qc4/"$filename""_R2_keep_ids.txt" -o usearch_qc4/"$filename""_R2_finalQC.fastq"
	done
	

Combine reads from all runs together for each sample:

	cat usearch_qc/V4_S1_L001_R1_R2_finalQC.fastq usearch_qc2/V4_S1_L001_R1_R2_finalQC.fastq usearch_qc4/v4_S1_L001_R1_R2_finalQC.fastq > V4.illumina.cat.fastq
	cat usearch_qc/V9_S2_L001_R1_R2_finalQC.fastq usearch_qc2/V9_S2_L001_R1_R2_finalQC.fastq usearch_qc3/v9_S2_L001_R1_R2_finalQC.fastq > V9.illumina.cat.fastq
	cat usearch_qc/V10_S3_L001_R1_R2_finalQC.fastq usearch_qc2/V10_S3_L001_R1_R2_finalQC.fastq usearch_qc3/v10_S3_L001_R1_R2_finalQC.fastq usearch_qc4/v10_S2_L001_R1_R2_finalQC.fastq > V10.illumina.cat.fastq
	cat usearch_qc/V11_S4_L001_R1_R2_finalQC.fastq usearch_qc2/V11_S4_L001_R1_R2_finalQC.fastq usearch_qc4/v11_S3_L001_R1_R2_finalQC.fastq > V11.illumina.cat.fastq
	cat usearch_qc/V12_S5_L001_R1_R2_finalQC.fastq usearch_qc2/V12_S5_L001_R1_R2_finalQC.fastq usearch_qc3/v12_S4_L001_R1_R2_finalQC.fastq usearch_qc4/v12_S4_L001_R1_R2_finalQC.fastq > V12.illumina.cat.fastq
	cat usearch_qc/V13_S6_L001_R1_R2_finalQC.fastq usearch_qc2/V13_S6_L001_R1_R2_finalQC.fastq usearch_qc3/v13_S5_L001_R1_R2_finalQC.fastq usearch_qc4/v13_S5_L001_R1_R2_finalQC.fastq > V13.illumina.cat.fastq
	
	expr $(cat V4.illumina.cat.fastq | wc -l) / 4

5576608

	expr $(cat V9.illumina.cat.fastq | wc -l) / 4

7415628

	expr $(cat V10.illumina.cat.fastq | wc -l) / 4

7358902

	expr $(cat V11.illumina.cat.fastq | wc -l) / 4

8385479

	expr $(cat V12.illumina.cat.fastq | wc -l) / 4

7455959

	expr $(cat V13.illumina.cat.fastq | wc -l) / 4

5527947



##rRNA Contamination
Convert all the viral finalQC FASTQ files to fasta to use for rRNA predictor. Download the rRNA predictor:
	
	wget http://weizhong-lab.ucsd.edu/meta_rna/rRNA_prediction.tar.bz2
	bzip2 -d rRNA_prediction.tar.bz2
	tar -xvf rRNA_prediction.tar
	chmod 777 -R rRNA_prediction/
	
For some reason I can only get this to run from within the examples folder.  So, move all the fasta files generated to the examples folder within rRNA predictor.
	
	mkdir rRNA_prediction/examples/e1/viral_input
	mkdir rRNA_prediction/examples/e1/viral_output
	cp prinseq_output/*.fna rRNA_prediction/examples/e1/viral_input/

Run the predictor for all samples.  Will likely need to put the full path to get it running:	

	cd rRNA_prediction/examples/e1
	export PATH=//work/samodha/canderson3/rRNA_prediction/rRNA_hmm_fs_wst_v0:$PATH
	/work/samodha/canderson3/rRNA_prediction/scripts/rRNA_hmm_run_wst_v0.pl  viral_input viral_output
	cd ..
	cd ..
	cd ..
	
Check the number of rRNA hits for each sample using custom script to parse outputs:
	
	wget https://raw.githubusercontent.com/chrisLanderson/rumen_virome/master/parse_rRNA_output.pl
	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.1_trimm_finalQC.fna.coord 

Total number of rRNA: 156
Total number of prokaryptic 16S and 23S rRNA: 74

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.2_trimm_finalQC.fna.coord 

Total number of rRNA: 5
Total number of prokaryptic 16S and 23S rRNA: 4

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.3_trimm_finalQC.fna.coord 

Total number of rRNA: 124
Total number of prokaryptic 16S and 23S rRNA: 111

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.4_trimm_finalQC.fna.coord 

Total number of rRNA: 837
Total number of prokaryptic 16S and 23S rRNA: 764

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.5_trimm_finalQC.fna.coord 

Total number of rRNA: 1084
Total number of prokaryptic 16S and 23S rRNA: 1063

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.6_trimm_finalQC.fna.coord 

Total number of rRNA: 43
Total number of prokaryptic 16S and 23S rRNA: 33

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.7_trimm_finalQC.fna.coord 

Total number of rRNA: 27
Total number of prokaryptic 16S and 23S rRNA: 22

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.8_trimm_finalQC.fna.coord 

Total number of rRNA: 142
Total number of prokaryptic 16S and 23S rRNA: 119

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.9_trimm_finalQC.fna.coord 

Total number of rRNA: 238
Total number of prokaryptic 16S and 23S rRNA: 232

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.10_trimm_finalQC.fna.coord 

Total number of rRNA: 582
Total number of prokaryptic 16S and 23S rRNA: 515

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.11_trimm_finalQC.fna.coord 

Total number of rRNA: 3068
Total number of prokaryptic 16S and 23S rRNA: 2929

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.12_trimm_finalQC.fna.coord 

Total number of rRNA: 419
Total number of prokaryptic 16S and 23S rRNA: 384

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.13_trimm_finalQC.fna.coord 

Total number of rRNA: 650
Total number of prokaryptic 16S and 23S rRNA: 630

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.14_trimm_finalQC.fna.coord 

Total number of rRNA: 93
Total number of prokaryptic 16S and 23S rRNA: 79

	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.15_trimm_finalQC.fna.coord 

Total number of rRNA: 368
Total number of prokaryptic 16S and 23S rRNA: 363
	

Cat all the .coord and fasta files together to get a total count for each:

	cat rRNA_prediction/examples/e1/viral_input/*.fna > VMG.cat.fasta
	grep -c ">" VMG.cat.fasta

13724814
	
	cat rRNA_prediction/examples/e1/viral_output/*.coord > rRNA_prediction/examples/e1/viral_output/VMG.cat.coord
	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.cat.coord 

Total number of rRNA: 7850
Total number of prokaryptic 16S and 23S rRNA: 7322


7850 total rRNA detected / 13724814 = 0.057%

7322 prokaryotic rRNA detected / 13724814 = 0.053%


In illumina dataset:

For some reason I can only get this to run from within the examples folder.  So, move all the fasta files generated to the examples folder within rRNA predictor.
	
	mkdir rRNA_prediction/examples/e1/viral_illumina_input
	mkdir rRNA_prediction/examples/e1/viral_illumina_output
	source anaconda/bin/activate rumenVirome
	convert_fastaqual_fastq.py -c fastq_to_fastaqual -f V4.illumina.cat.fastq 
	convert_fastaqual_fastq.py -c fastq_to_fastaqual -f V9.illumina.cat.fastq 
	convert_fastaqual_fastq.py -c fastq_to_fastaqual -f V10.illumina.cat.fastq 
	convert_fastaqual_fastq.py -c fastq_to_fastaqual -f V11.illumina.cat.fastq 
	convert_fastaqual_fastq.py -c fastq_to_fastaqual -f V12.illumina.cat.fastq 
	convert_fastaqual_fastq.py -c fastq_to_fastaqual -f V13.illumina.cat.fastq 	
	
	cp *.illumina.cat.fna rRNA_prediction/examples/e1/viral_illumina_input/
	

Run the predictor for all samples.  Will likely need to put the full path to get it running:	

	cd rRNA_prediction/examples/e1
	export PATH=//work/samodha/canderson3/rRNA_prediction/rRNA_hmm_fs_wst_v0:$PATH
	/work/samodha/canderson3/rRNA_prediction/scripts/rRNA_hmm_run_wst_v0.pl  viral_illumina_input viral_illumina_output
	cd ..
	cd ..
	cd ..
	
Check the number of rRNA hits for each sample using custom script to parse outputs:
	
	wget https://raw.githubusercontent.com/chrisLanderson/rumen_virome/master/parse_rRNA_output.pl
	./parse_rRNA_output.pl -rrna rRNA_prediction/examples/e1/viral_output/VMG.1_trimm_finalQC.fna.coord 




##K-mer Profiles
Use khmer software to compare k-mer profiles of different metagenomes:

Now, load all viral metagenomes into counting:

    cd /work/samodha/canderson3/prinseq_output
    for f in *_finalQC.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        load-into-counting.py -k 20 -N 4 -x 8e9 --report-total-kmers -s tsv "$filename""_k20.kh" $f
    done


Now, load all total metagenomes into counting:

    
    cd cd_hit_454_output_total
    for f in *.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        load-into-counting.py -k 20 -N 4 -x 8e9 --report-total-kmers -s tsv "$filename""_k20.kh" $f
    done

Load all viral illumina datasets into counting:

	##kmer profiles for illumina data

	load-into-counting.py -k 20 -N 4 -x 8e9 --report-total-kmers -s tsv V4.illumina.kh V4.illumina.cat.fastq
	load-into-counting.py -k 20 -N 4 -x 8e9 --report-total-kmers -s tsv V9.illumina.kh V9.illumina.cat.fastq
	load-into-counting.py -k 20 -N 4 -x 8e9 --report-total-kmers -s tsv V10.illumina.kh V10.illumina.cat.fastq
	load-into-counting.py -k 20 -N 4 -x 8e9 --report-total-kmers -s tsv V11.illumina.kh V11.illumina.cat.fastq
	load-into-counting.py -k 20 -N 4 -x 8e9 --report-total-kmers -s tsv V12.illumina.kh V12.illumina.cat.fastq
	load-into-counting.py -k 20 -N 4 -x 8e9 --report-total-kmers -s tsv V13.illumina.kh V13.illumina.cat.fastq

Now, need to look at how reads are shared across samples in a pairwise fashion within each dataset. Currently, we need to fetch the script from authors github and do some altering to it:

	wget https://raw.githubusercontent.com/qingpeng/igs-diversity/master/scripts/get_comb_multi_old_median_kmer.py
	chmod 775 get_comb_multi_old_median_kmer.py 

In case, the script has been modified or moved, the script I used in the analysis is available at:


To get it to work in our current khmer environment, use vi and alter the first line to be as follows, minus the single quotes (will need to change the path to match yours):

'#delete first line
'#remove if K <= line and move everything over??
	
	
Setup the config file used by the get_comb_multi.py script. Ion viral first:
	
	cd prinseq_output/
	ls *k20.kh | awk '{ ORS=" "; print; }' > config.txt
	printf "\n" >> config.txt
	ls *finalQC.fastq | awk '{ ORS=" "; print; }' >> config.txt
	printf "\n" >> config.txt
	printf "55000000" >> config.txt
	
	source ../anaconda/bin/activate rumenVirome
	python ../get_comb_multi_old_median_kmer.py config.txt

	cd ..
	
Illumina Viral:

	ls *illumina.kh | awk '{ ORS=" "; print; }' > config.txt
	printf "\n" >> config.txt
	printf "V10.illumina.cat.fastq V11.illumina.cat.fastq V12.illumina.cat.fastq V13.illumina.cat.fastq V4.illumina.cat.fastq V9.illumina.cat.fastq" >> config.txt
	printf "\n" >> config.txt
	printf "200000000" >> config.txt
	
	source anaconda/bin/activate rumenVirome
	python get_comb_multi_old_median_kmer.py config.txt

Total:

	cd cd_hit_454_output_total/
	ls *k20.kh | awk '{ ORS=" "; print; }' > config.txt
	printf "\n" >> config.txt
	ls *.fastq | awk '{ ORS=" "; print; }' >> config.txt
	printf "\n" >> config.txt
	printf "200000000" >> config.txt
	
	source ../anaconda/bin/activate rumenVirome
	python ../get_comb_multi_old_median_kmer.py config.txt

	cd ..
	


To try and eliminate the imapct of some sequencing errors on the results, say that a read must be shared with at least one other sample in the dataset:

	wget https://raw.githubusercontent.com/chrisLanderson/rumen_virome/master/khmer_multi_threshold.pl
	chmod 775 khmer_multi_threshold.pl
	
	./khmer_multi_threshold.pl -khmer_multi_dir= -threshold=1
	
	./khmer_multi_threshold.pl -khmer_multi_dir= -threshold=1
	
	./khmer_multi_threshold.pl -khmer_multi_dir= -threshold=1


Look at how reads are shared within a diet and within an individual:
Total:

	./khmer_multi_grouping.pl -group=diet -khmer_multi_dir= -threshold=3 -config=config.txt 
	./khmer_multi_grouping.pl -group=animal -khmer_multi_dir= -threshold=1 -config=config.txt 

Ion Viral:


Illumina Viral:


Compare diets pairwise to observe the shared number of reads between them:
Total:

	./khmer_multi_diet_overlap.pl -khmer_multi_dir= -threshold=2 -config=config.txt -diet1=40MDGS -diet2=Control


Ion Viral:


Illumina viral:






Overlap of Viral and total:



##Assembly
Viral:

	wget http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0-Linux.tar.gz
	tar -zxvf SPAdes-3.5.0-Linux.tar.gz
	cat prinseq_output/*_finalQC.fastq > VMG.cat.fastq
	python SPAdes-3.5.0-Linux/bin/spades.py -k 21,33,55,77,99,127 --only-assembler --sc -s VMG.cat.fastq -m 60000 -t 2 -o vmg_cat_sc

Total:

	cat cd_hit_454_output_total/*.fastq > BMG.cat.fastq
	python SPAdes-3.5.0-Linux/bin/spades.py -k 21,33,55,77,99,127 --only-assembler --sc -s BMG.cat.fastq -m 200000 -t 2 -o bmg_cat_sc


Illumina:

	cat *.illumina.cat.fastq > VMG.illumina.cat.fastq
	python SPAdes-3.5.0-Linux/bin/spades.py -k 21,33,55,77,99,127 --sc -s VMG.illumina.cat.fastq -m 250000 -t 2 -o vmg_illumina_cat_sc



##Predict ORFs
Download prodigal:

	wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.2/prodigal.linux
	
	wget https://github.com/hyattpd/Prodigal/archive/v2.6.2.zip
	unzip v2.6.2.zip
	cur=$(pwd)
	cd Prodigal-2.6.2/
	make install INSTALLDIR=$cur
	cd ..
	
Viral:
	
	mv vmg_cat_sc/contigs.fasta ./
	mv contigs.fasta vmg.contigs.fasta
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fasta vmg.contigs.fasta -min_len 300 -out_format 1 -out_good vmg.contigs.filter
	./prodigal -i vmg.contigs.filter.fasta -d vmg.orfs_nt.fasta -a vmg.orfs_aa.fasta -p meta 

Filter the ORFs based on length of nt, then use those sequences as a guide to remove the short seqeunces from aa file:
	
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fasta vmg.orfs_nt.fasta -min_len 300 -out_format 1 -out_good vmg.orfs_nt.filter
	filter_fasta.py -f vmg.orfs_aa.fasta -a vmg.orfs_nt.filter.fasta -o vmg.orfs_aa.filter.fasta
	
Total:

	mv bmg_cat_sc/contigs.fasta ./
	mv contigs.fasta bmg.contigs.fasta
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fasta bmg.contigs.fasta -min_len 300 -out_format 1 -out_good bmg.contigs.filter
	./prodigal -i bmg.contigs.filter.fasta -d bmg.orfs_nt.fasta -a bmg.orfs_aa.fasta -p meta 

Filter the ORFs based on length of nt, then use those sequences as a guide to remove the short seqeunces from aa file:
	
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fasta bmg.orfs_nt.fasta -min_len 300 -out_format 1 -out_good bmg.orfs_nt.filter
	filter_fasta.py -f bmg.orfs_aa.fasta -a bmg.orfs_nt.filter.fasta -o bmg.orfs_aa.filter.fasta
	
Illumina:

	mv vmg_illumina_cat_sc/contigs.fasta ./
	mv contigs.fasta vmg.illumina.contigs.fasta
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fasta vmg.illumina.contigs.fasta -min_len 300 -out_format 1 -out_good vmg.illumina.contigs.filter
	./prodigal -i vmg.illumina.contigs.filter.fasta -d vmg.illumina.orfs_nt.fasta -a vmg.illumina.orfs_aa.fasta -p meta 

Filter the ORFs based on length of nt, then use those sequences as a guide to remove the short seqeunces from aa file:
	
	perl prinseq-lite-0.20.4/./prinseq-lite.pl -fasta vmg.illumina.orfs_nt.fasta -min_len 300 -out_format 1 -out_good vmg.illumina.orfs_nt.filter
	filter_fasta.py -f vmg.illumina.orfs_aa.fasta -a vmg.illumina.orfs_nt.filter.fasta -o vmg.illumina.orfs_aa.filter.fasta
	

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
	
	./usearch7.0.10 -ublast bmg.orfs_aa.filter.fasta -db keggV69.genes.pep.udb -evalue 1e-5 -blast6out bmg.ublast.kegg.txt -strand both -top_hits_only -threads 15


Search Total against PHAST:
	
	./usearch7.0.10 -ublast bmg.orfs_aa.filter.fasta -db prophage_virus.udb -evalue 1e-5 -blast6out bmg.ublast.phast.txt -strand both -top_hits_only -threads 15


Search Illumina Viral against KEGG:

	./usearch7.0.10 -ublast vmg.illumina.orfs_aa.filter.fasta -db keggV69.genes.pep.udb -evalue 1e-5 -blast6out vmg.illumina.ublast.kegg.txt -strand both -top_hits_only -threads 15

Search Illumina Viral against PHAST:
	
	./usearch7.0.10 -ublast vmg.illumina.orfs_aa.filter.fasta -db prophage_virus.udb -evalue 1e-5 -blast6out vmg.illumina.ublast.phast.txt -strand both -top_hits_only -threads 15



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

	mkdir vmg_illumina_orfs_bowtie
	bowtie2-2.2.5/bowtie2-build vmg.illumina.orfs_nt.filter.fasta vmg_illumina_orfs_bowtie/vmg_illumina_orfs_bowtie_db
	
Align all viral reads to nt ORFs:

    for f in prinseq_output/*.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        bowtie2-2.2.5/bowtie2 -U $f --end-to-end --sensitive -x vmg_orfs_bowtie/vmg_orfs_bowtie_db -S vmg_orfs_bowtie/$filename.psl --un vmg_orfs_bowtie/$filename.unaligned.txt --al vmg_orfs_bowtie/$filename.aligned.txt 
    done
 
All total metagenome reads to nt ORFs:
	
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


Number of top hits from BLAST file: 95153

Alignment count for VMG.10_trimm: 335426
Total Read count VMG.10_trimm: 659514
Redundant gene count with aligned reads for VMG.10_trimm: 17361
Non-Redundant gene count with aligned reads for VMG.10_trimm: 13233
Total abundance in gene hits from aligned reads for VMG.10_trimm: 126334
Non-Redundant ko count with aligned reads for VMG.10_trimm: 1880
Number of nr genes with ko hits: 5133
Total abundance in ko hits from aligned reads for VMG.10_trimm: 22533
Number of nr genes without ko hits: 8100
Total abundance of reads from genes without ko hits: 103801

Alignment count for VMG.11_trimm: 1719015
Total Read count VMG.11_trimm: 2890819
Redundant gene count with aligned reads for VMG.11_trimm: 48065
Non-Redundant gene count with aligned reads for VMG.11_trimm: 24214
Total abundance in gene hits from aligned reads for VMG.11_trimm: 675510
Non-Redundant ko count with aligned reads for VMG.11_trimm: 1648
Number of nr genes with ko hits: 9786
Total abundance in ko hits from aligned reads for VMG.11_trimm: 152395
Number of nr genes without ko hits: 14428
Total abundance of reads from genes without ko hits: 523115

Alignment count for VMG.12_trimm: 230192
Total Read count VMG.12_trimm: 497309
Redundant gene count with aligned reads for VMG.12_trimm: 11896
Non-Redundant gene count with aligned reads for VMG.12_trimm: 8813
Total abundance in gene hits from aligned reads for VMG.12_trimm: 64893
Non-Redundant ko count with aligned reads for VMG.12_trimm: 1138
Number of nr genes with ko hits: 3314
Total abundance in ko hits from aligned reads for VMG.12_trimm: 11478
Number of nr genes without ko hits: 5499
Total abundance of reads from genes without ko hits: 53415

Alignment count for VMG.13_trimm: 550866
Total Read count VMG.13_trimm: 958298
Redundant gene count with aligned reads for VMG.13_trimm: 19542
Non-Redundant gene count with aligned reads for VMG.13_trimm: 13660
Total abundance in gene hits from aligned reads for VMG.13_trimm: 185322
Non-Redundant ko count with aligned reads for VMG.13_trimm: 1750
Number of nr genes with ko hits: 4885
Total abundance in ko hits from aligned reads for VMG.13_trimm: 27712
Number of nr genes without ko hits: 8775
Total abundance of reads from genes without ko hits: 157610

Alignment count for VMG.14_trimm: 382173
Total Read count VMG.14_trimm: 652788
Redundant gene count with aligned reads for VMG.14_trimm: 6537
Non-Redundant gene count with aligned reads for VMG.14_trimm: 4105
Total abundance in gene hits from aligned reads for VMG.14_trimm: 170389
Non-Redundant ko count with aligned reads for VMG.14_trimm: 392
Number of nr genes with ko hits: 924
Total abundance in ko hits from aligned reads for VMG.14_trimm: 22306
Number of nr genes without ko hits: 3181
Total abundance of reads from genes without ko hits: 148083

Alignment count for VMG.15_trimm: 204269
Total Read count VMG.15_trimm: 356082
Redundant gene count with aligned reads for VMG.15_trimm: 7150
Non-Redundant gene count with aligned reads for VMG.15_trimm: 5441
Total abundance in gene hits from aligned reads for VMG.15_trimm: 86831
Non-Redundant ko count with aligned reads for VMG.15_trimm: 1021
Number of nr genes with ko hits: 2009
Total abundance in ko hits from aligned reads for VMG.15_trimm: 13781
Number of nr genes without ko hits: 3432
Total abundance of reads from genes without ko hits: 73050

Alignment count for VMG.1_trimm: 388672
Total Read count VMG.1_trimm: 661076
Redundant gene count with aligned reads for VMG.1_trimm: 8354
Non-Redundant gene count with aligned reads for VMG.1_trimm: 5101
Total abundance in gene hits from aligned reads for VMG.1_trimm: 152237
Non-Redundant ko count with aligned reads for VMG.1_trimm: 514
Number of nr genes with ko hits: 1176
Total abundance in ko hits from aligned reads for VMG.1_trimm: 33573
Number of nr genes without ko hits: 3925
Total abundance of reads from genes without ko hits: 118664

Alignment count for VMG.2_trimm: 63883
Total Read count VMG.2_trimm: 141850
Redundant gene count with aligned reads for VMG.2_trimm: 1126
Non-Redundant gene count with aligned reads for VMG.2_trimm: 852
Total abundance in gene hits from aligned reads for VMG.2_trimm: 6453
Non-Redundant ko count with aligned reads for VMG.2_trimm: 84
Number of nr genes with ko hits: 135
Total abundance in ko hits from aligned reads for VMG.2_trimm: 667
Number of nr genes without ko hits: 717
Total abundance of reads from genes without ko hits: 5786

Alignment count for VMG.3_trimm: 1005238
Total Read count VMG.3_trimm: 1472309
Redundant gene count with aligned reads for VMG.3_trimm: 13671
Non-Redundant gene count with aligned reads for VMG.3_trimm: 7270
Total abundance in gene hits from aligned reads for VMG.3_trimm: 367298
Non-Redundant ko count with aligned reads for VMG.3_trimm: 608
Number of nr genes with ko hits: 1717
Total abundance in ko hits from aligned reads for VMG.3_trimm: 45852
Number of nr genes without ko hits: 5553
Total abundance of reads from genes without ko hits: 321446

Alignment count for VMG.4_trimm: 1159901
Total Read count VMG.4_trimm: 2026959
Redundant gene count with aligned reads for VMG.4_trimm: 24337
Non-Redundant gene count with aligned reads for VMG.4_trimm: 14186
Total abundance in gene hits from aligned reads for VMG.4_trimm: 450491
Non-Redundant ko count with aligned reads for VMG.4_trimm: 1366
Number of nr genes with ko hits: 4845
Total abundance in ko hits from aligned reads for VMG.4_trimm: 81635
Number of nr genes without ko hits: 9341
Total abundance of reads from genes without ko hits: 368856

Alignment count for VMG.5_trimm: 1345061
Total Read count VMG.5_trimm: 2166417
Redundant gene count with aligned reads for VMG.5_trimm: 25115
Non-Redundant gene count with aligned reads for VMG.5_trimm: 13928
Total abundance in gene hits from aligned reads for VMG.5_trimm: 605441
Non-Redundant ko count with aligned reads for VMG.5_trimm: 1287
Number of nr genes with ko hits: 4675
Total abundance in ko hits from aligned reads for VMG.5_trimm: 122338
Number of nr genes without ko hits: 9253
Total abundance of reads from genes without ko hits: 483103

Alignment count for VMG.6_trimm: 185311
Total Read count VMG.6_trimm: 305262
Redundant gene count with aligned reads for VMG.6_trimm: 3711
Non-Redundant gene count with aligned reads for VMG.6_trimm: 2524
Total abundance in gene hits from aligned reads for VMG.6_trimm: 54237
Non-Redundant ko count with aligned reads for VMG.6_trimm: 245
Number of nr genes with ko hits: 494
Total abundance in ko hits from aligned reads for VMG.6_trimm: 5160
Number of nr genes without ko hits: 2030
Total abundance of reads from genes without ko hits: 49077

Alignment count for VMG.7_trimm: 325582
Total Read count VMG.7_trimm: 530173
Redundant gene count with aligned reads for VMG.7_trimm: 5940
Non-Redundant gene count with aligned reads for VMG.7_trimm: 3668
Total abundance in gene hits from aligned reads for VMG.7_trimm: 111078
Non-Redundant ko count with aligned reads for VMG.7_trimm: 273
Number of nr genes with ko hits: 707
Total abundance in ko hits from aligned reads for VMG.7_trimm: 13728
Number of nr genes without ko hits: 2961
Total abundance of reads from genes without ko hits: 97350

Alignment count for VMG.8_trimm: 349140
Total Read count VMG.8_trimm: 563401
Redundant gene count with aligned reads for VMG.8_trimm: 5475
Non-Redundant gene count with aligned reads for VMG.8_trimm: 3933
Total abundance in gene hits from aligned reads for VMG.8_trimm: 113900
Non-Redundant ko count with aligned reads for VMG.8_trimm: 468
Number of nr genes with ko hits: 980
Total abundance in ko hits from aligned reads for VMG.8_trimm: 12425
Number of nr genes without ko hits: 2953
Total abundance of reads from genes without ko hits: 101475

Alignment count for VMG.9_trimm: 480977
Total Read count VMG.9_trimm: 748778
Redundant gene count with aligned reads for VMG.9_trimm: 14351
Non-Redundant gene count with aligned reads for VMG.9_trimm: 9052
Total abundance in gene hits from aligned reads for VMG.9_trimm: 198460
Non-Redundant ko count with aligned reads for VMG.9_trimm: 989
Number of nr genes with ko hits: 2944
Total abundance in ko hits from aligned reads for VMG.9_trimm: 26403
Number of nr genes without ko hits: 6108
Total abundance of reads from genes without ko hits: 172057


Total nucleotides for all samples combined: 2284667869
Average nucleotides per samples: 152311191.266667
Number of ORFs with assignments printed to orf_assignment.txt: 32928

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

	merge_otu_tables.py -i viral_blast2tsv/VMG.1_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.2_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.3_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.4_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.5_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.6_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.7_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.8_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.9_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.10_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.11_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.12_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.13_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.14_trimm.ko_corrected_abundances.biom,viral_blast2tsv/VMG.15_trimm.ko_corrected_abundances.biom -o viral_blast2tsv/VMG.ko_corrected_abundances.biom

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

Number of unique differential k numbers found with ORFs: 2877
Total number of PHAST top hits: 36850
Contigs with phast and KO hit: 6810 
ORFs on a contig with phast hit and ORF has a KO: 8520
K numbers printed that were on a contig with PHAST hit and KO: 895

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


This should push the network to cytoscape 2.8.2
Apply force-directed layout
Saved file as refnet.cys
Open network in cytoscape 3.1.0 (File --> open)

Select the largest connected portion of the network (Tools --> network analyzer --> subnetwork creation --> extract connected components, select first component --> extract)


This should generate a network with 3616 nodes.


Analyzed subset network as directed (Tools --> network analyzer --> network analysis --> Analyze network --> treat the network as directed --> ok )


This will take awhile.

File -- Export --> Table --> node(1) default node

Save as:
refnet_attributes.csv

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

Apply force-directed layout
Saved file as ssn.cys
Open network in cytoscape 3.1.0 (File --> open)

Select the largest connected portion of the network (Tools --> network analyzer --> subnetwork creation --> extract connected components, select first component --> extract)


This should generate a network with 1304 nodes.


Analyzed subset network as directed (Tools --> network analyzer --> network analysis --> Analyze network --> treat the network as directed --> ok )


This will take awhile.

File -- Export --> Table --> node(1) default node

Save as:
ssn_attributes.csv

Open in excel and convert to tab delimited text file. Save as ssn_attributes.txt

##Topology of Differential KOs in Metabolic Network

Want to investigate the topology of differential features in both the total metagenome and viral metagenome in the context of a metabolic network.  In our metabolic network KOs (enzyme families) are nodes and their connections are established by KEGG and implemented into a network previously with the mmnet package in R.

Need to take the network attribute tables and identify the topological features for the differential KOs and compare the values of differential to non-differntial and viral differential to total differential.

	mkdir ref_network_plots
	mkdir ssn_network_plots
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
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_degree <- ggplot(vmg_bmg_non, aes(x=Dataset, y=EdgeCount)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Total Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_degree, file="ref_network_plots/degree_grey.pdf", w=5, h=5)

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
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_indegree <- ggplot(vmg_bmg_non, aes(x=Dataset, y=Indegree)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "In-Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_indegree, file="ref_network_plots/indegree_grey.pdf", w=5, h=5)

Out-Degree:

	VMG <- boxplot.stats(vmg_diff_topo$Outdegree)$stats
	BMG <- boxplot.stats(bmg_diff_topo$Outdegree)$stats
	Non <- boxplot.stats(non_diff_topo$Outdegree)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$Outdegree)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_outdegree <- ggplot(vmg_bmg_non, aes(x=Dataset, y=Outdegree)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Out-Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_outdegree, file="ref_network_plots/outdegree_grey.pdf", w=5, h=5)

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
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_avgpath <- ggplot(vmg_bmg_non, aes(x=Dataset, y=AverageShortestPathLength)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Average Shortest Path Length",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_avgpath, file="ref_network_plots/averagepathlength_grey.pdf", w=5, h=5)

Betweenness Centrality:

	VMG <- boxplot.stats(vmg_diff_topo$BetweennessCentrality)$stats
	BMG <- boxplot.stats(bmg_diff_topo$BetweennessCentrality)$stats
	Non <- boxplot.stats(non_diff_topo$BetweennessCentrality)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","a,b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$BetweennessCentrality)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_bc <- ggplot(vmg_bmg_non, aes(x=Dataset, y=BetweennessCentrality)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Betweenness Centrality",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_bc, file="ref_network_plots/betweennesscentrality_grey.pdf", w=5, h=5)

Closeness Centrality:

	VMG <- boxplot.stats(vmg_diff_topo$ClosenessCentrality)$stats
	BMG <- boxplot.stats(bmg_diff_topo$ClosenessCentrality)$stats
	Non <- boxplot.stats(non_diff_topo$ClosenessCentrality)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$ClosenessCentrality)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_cc <- ggplot(vmg_bmg_non, aes(x=Dataset, y=ClosenessCentrality)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Closeness Centrality",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_cc, file="ref_network_plots/closenesscentrality_grey.pdf", w=5, h=5)


Neighborhood Connectivity:

	VMG <- boxplot.stats(vmg_diff_topo$NeighborhoodConnectivity)$stats
	BMG <- boxplot.stats(bmg_diff_topo$NeighborhoodConnectivity)$stats
	Non <- boxplot.stats(non_diff_topo$NeighborhoodConnectivity)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$NeighborhoodConnectivity)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_nc <- ggplot(vmg_bmg_non, aes(x=Dataset, y=NeighborhoodConnectivity)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Neighborhood Connectivity",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_nc, file="ref_network_plots/neighborhoodconnectivity_grey.pdf", w=5, h=5)

Eccentricity:

	VMG <- boxplot.stats(vmg_diff_topo$Eccentricity)$stats
	BMG <- boxplot.stats(bmg_diff_topo$Eccentricity)$stats
	Non <- boxplot.stats(non_diff_topo$Eccentricity)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$Eccentricity)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_ecc <- ggplot(vmg_bmg_non, aes(x=Dataset, y=Eccentricity)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Eccentricity",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_ecc, file="ref_network_plots/eccentricity_grey.pdf", w=5, h=5)

Clustering Coefficient:

	VMG <- boxplot.stats(vmg_diff_topo$ClusteringCoefficient)$stats
	BMG <- boxplot.stats(bmg_diff_topo$ClusteringCoefficient)$stats
	Non <- boxplot.stats(non_diff_topo$ClusteringCoefficient)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("a","a","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$ClusteringCoefficient)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_clust <- ggplot(vmg_bmg_non, aes(x=Dataset, y=ClusteringCoefficient)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Clustering Coefficient",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_clust, file="ref_network_plots/clustering_coeffecient_grey.pdf", w=5, h=5)
	
All on one 9 panel plot:

	install.packages("gridExtra")
	install.packages("grid")
	library(gridExtra)
	library(grid)

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
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_degree <- ggplot(vmg_bmg_non, aes(x=Dataset, y=EdgeCount)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Total Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.ticks = element_blank())


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
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_indegree <- ggplot(vmg_bmg_non, aes(x=Dataset, y=Indegree)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "In-Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.ticks = element_blank())


	VMG <- boxplot.stats(vmg_diff_topo$Outdegree)$stats
	BMG <- boxplot.stats(bmg_diff_topo$Outdegree)$stats
	Non <- boxplot.stats(non_diff_topo$Outdegree)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$Outdegree)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_outdegree <- ggplot(vmg_bmg_non, aes(x=Dataset, y=Outdegree)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Out-Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.ticks = element_blank())


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
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_avgpath <- ggplot(vmg_bmg_non, aes(x=Dataset, y=AverageShortestPathLength)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Average Shortest Path Length",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.ticks = element_blank())


	VMG <- boxplot.stats(vmg_diff_topo$BetweennessCentrality)$stats
	BMG <- boxplot.stats(bmg_diff_topo$BetweennessCentrality)$stats
	Non <- boxplot.stats(non_diff_topo$BetweennessCentrality)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","a,b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$BetweennessCentrality)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_bc <- ggplot(vmg_bmg_non, aes(x=Dataset, y=BetweennessCentrality)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Betweenness Centrality",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.ticks = element_blank())


	VMG <- boxplot.stats(vmg_diff_topo$ClosenessCentrality)$stats
	BMG <- boxplot.stats(bmg_diff_topo$ClosenessCentrality)$stats
	Non <- boxplot.stats(non_diff_topo$ClosenessCentrality)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$ClosenessCentrality)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_cc <- ggplot(vmg_bmg_non, aes(x=Dataset, y=ClosenessCentrality)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Closeness Centrality",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.ticks = element_blank())



	VMG <- boxplot.stats(vmg_diff_topo$NeighborhoodConnectivity)$stats
	BMG <- boxplot.stats(bmg_diff_topo$NeighborhoodConnectivity)$stats
	Non <- boxplot.stats(non_diff_topo$NeighborhoodConnectivity)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$NeighborhoodConnectivity)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_nc <- ggplot(vmg_bmg_non, aes(x=Dataset, y=NeighborhoodConnectivity)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Neighborhood Connectivity",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.ticks = element_blank())


	VMG <- boxplot.stats(vmg_diff_topo$Eccentricity)$stats
	BMG <- boxplot.stats(bmg_diff_topo$Eccentricity)$stats
	Non <- boxplot.stats(non_diff_topo$Eccentricity)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$Eccentricity)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_ecc <- ggplot(vmg_bmg_non, aes(x=Dataset, y=Eccentricity)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Eccentricity",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.ticks = element_blank())


	VMG <- boxplot.stats(vmg_diff_topo$ClusteringCoefficient)$stats
	BMG <- boxplot.stats(bmg_diff_topo$ClusteringCoefficient)$stats
	Non <- boxplot.stats(non_diff_topo$ClusteringCoefficient)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("a","a","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non$ClusteringCoefficient)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	
	plot_clust <- ggplot(vmg_bmg_non, aes(x=Dataset, y=ClusteringCoefficient)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Clustering Coefficient",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=3485)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=47)")) +
	theme(axis.ticks = element_blank())

	pdf(file="ref_network_plots/network_topology_panel_grey.pdf",w=12,h=12)
	grid.arrange(plot_degree,plot_indegree,plot_outdegree,plot_bc,plot_cc,plot_clust,plot_avgpath,plot_nc,plot_ecc,legend,ncol=3,nrow = 3)
	dev.off()



Repeat all of above for the ssn network:

	ssn_topo <- read.table("ssn_attributes.txt", header=TRUE, sep="\t")
	vmg_diff_topo_ssn <- merge(x=vmg_diff, y=ssn_topo, by="canonicalName")
	bmg_diff_topo_ssn <- merge(x=bmg_diff, y=ssn_topo, by="canonicalName")

There is one less KO in the BMG dataset because it was removed when extracting the largest component of the network.

	total_topo_ssn <- merge(x=ssn_topo, y=diff, all.x=TRUE)
	total_topo_ssn[is.na(total_topo_ssn)] <- 0
	non_diff_topo_ssn <- subset(total_topo_ssn, total_topo_ssn$LDA == 0)
	
	vmg_diff_topo_ssn[vmg_diff_topo_ssn=="TRUE"] <- "VMG"
	bmg_diff_topo_ssn[bmg_diff_topo_ssn=="TRUE"] <- "BMG"
	non_diff_topo_ssn$Dataset[non_diff_topo_ssn$Dataset == 0] <- "Non"
	
	vmg_bmg_ssn <- rbind(vmg_diff_topo_ssn, bmg_diff_topo_ssn)
	vmg_non_ssn <- rbind(vmg_diff_topo_ssn, non_diff_topo_ssn)
	bmg_non_ssn <- rbind(bmg_diff_topo_ssn, non_diff_topo_ssn) 
	
	vmg_non_ssn_edge <- wilcox.test(EdgeCount ~ Dataset, data=vmg_non_ssn, paired=FALSE)    
	vmg_non_ssn_in <- wilcox.test(Indegree ~ Dataset, data=vmg_non_ssn, paired=FALSE)   
	vmg_non_ssn_out <- wilcox.test(Outdegree ~ Dataset, data=vmg_non_ssn, paired=FALSE)
	vmg_non_ssn_path <- wilcox.test(AverageShortestPathLength ~ Dataset, data=vmg_non_ssn, paired=FALSE)
	vmg_non_ssn_bc <- wilcox.test(BetweennessCentrality ~ Dataset, data=vmg_non_ssn, paired=FALSE)
	vmg_non_ssn_cc <- wilcox.test(ClosenessCentrality ~ Dataset, data=vmg_non_ssn, paired=FALSE)
	vmg_non_ssn_clust <- wilcox.test(ClusteringCoefficient ~ Dataset, data=vmg_non_ssn, paired=FALSE)
	vmg_non_ssn_ecc <- wilcox.test(Eccentricity ~ Dataset, data=vmg_non_ssn, paired=FALSE)
	vmg_non_ssn_nc <- wilcox.test(NeighborhoodConnectivity ~ Dataset, data=vmg_non_ssn, paired=FALSE)

	bmg_non_ssn_edge <- wilcox.test(EdgeCount ~ Dataset, data=bmg_non_ssn, paired=FALSE)
	bmg_non_ssn_in <- wilcox.test(Indegree ~ Dataset, data=bmg_non_ssn, paired=FALSE)
	bmg_non_ssn_out <- wilcox.test(Outdegree ~ Dataset, data=bmg_non_ssn, paired=FALSE)
	bmg_non_ssn_path <- wilcox.test(AverageShortestPathLength ~ Dataset, data=bmg_non_ssn, paired=FALSE)
	bmg_non_ssn_bc <- wilcox.test(BetweennessCentrality ~ Dataset, data=bmg_non_ssn, paired=FALSE)
	bmg_non_ssn_cc <- wilcox.test(ClosenessCentrality ~ Dataset, data=bmg_non_ssn, paired=FALSE)
	bmg_non_ssn_clust <- wilcox.test(ClusteringCoefficient ~ Dataset, data=bmg_non_ssn, paired=FALSE)
	bmg_non_ssn_ecc <- wilcox.test(Eccentricity ~ Dataset, data=bmg_non_ssn, paired=FALSE)
	bmg_non_ssn_nc <- wilcox.test(NeighborhoodConnectivity ~ Dataset, data=bmg_non_ssn, paired=FALSE)

	vmg_bmg_ssn_edge <- wilcox.test(EdgeCount ~ Dataset, data=vmg_bmg_ssn, paired=FALSE)
	vmg_bmg_ssn_in <- wilcox.test(Indegree ~ Dataset, data=vmg_bmg_ssn, paired=FALSE)
	vmg_bmg_ssn_out <- wilcox.test(Outdegree ~ Dataset, data=vmg_bmg_ssn, paired=FALSE)
	vmg_bmg_ssn_path <- wilcox.test(AverageShortestPathLength ~ Dataset, data=vmg_bmg_ssn, paired=FALSE)
	vmg_bmg_ssn_bc <- wilcox.test(BetweennessCentrality ~ Dataset, data=vmg_bmg_ssn, paired=FALSE)
	vmg_bmg_ssn_cc <- wilcox.test(ClosenessCentrality ~ Dataset, data=vmg_bmg_ssn, paired=FALSE)
	vmg_bmg_ssn_clust <- wilcox.test(ClusteringCoefficient ~ Dataset, data=vmg_bmg_ssn, paired=FALSE)
	vmg_bmg_ssn_ecc <- wilcox.test(Eccentricity ~ Dataset, data=vmg_bmg_ssn, paired=FALSE)
	vmg_bmg_ssn_nc <- wilcox.test(NeighborhoodConnectivity ~ Dataset, data=vmg_bmg_ssn, paired=FALSE)

Plot these ssn results:

Total Degree:

	VMG <- boxplot.stats(vmg_diff_topo_ssn$EdgeCount)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$EdgeCount)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$EdgeCount)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("a","a","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$EdgeCount)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_degree <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=EdgeCount)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Total Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_degree, file="ssn_network_plots/degree_grey.pdf", w=5, h=5)
	
In-Degree:

	VMG <- boxplot.stats(vmg_diff_topo_ssn$Indegree)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$Indegree)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$Indegree)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","a,b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$Indegree)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_indegree <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=Indegree)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "In-Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_indegree, file="ssn_network_plots/indegree_grey.pdf", w=5, h=5)
	
Out-Degree:

	VMG <- boxplot.stats(vmg_diff_topo_ssn$Outdegree)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$Outdegree)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$Outdegree)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("a","a","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$Outdegree)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_outdegree <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=Outdegree)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Out-Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_outdegree, file="ssn_network_plots/outdegree_grey.pdf", w=5, h=5)
	
Average Path:
	
	VMG <- boxplot.stats(vmg_diff_topo_ssn$AverageShortestPathLength)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$AverageShortestPathLength)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$AverageShortestPathLength)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","a,b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$AverageShortestPathLength)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_avgpath <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=AverageShortestPathLength)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Average Shortest Path Length",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_avgpath, file="ssn_network_plots/averagepathlength_grey.pdf", w=5, h=5)

Betweenness Centrality:

	VMG <- boxplot.stats(vmg_diff_topo_ssn$BetweennessCentrality)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$BetweennessCentrality)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$BetweennessCentrality)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","a,b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$BetweennessCentrality)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_bc <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=BetweennessCentrality)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Betweenness Centrality",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_bc, file="ssn_network_plots/betweennesscentrality_grey.pdf", w=5, h=5)

Closeness Centrality:

	VMG <- boxplot.stats(vmg_diff_topo_ssn$ClosenessCentrality)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$ClosenessCentrality)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$ClosenessCentrality)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("a","a","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$ClosenessCentrality)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_cc <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=ClosenessCentrality)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Closeness Centrality",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_cc, file="ssn_network_plots/closenesscentrality_grey.pdf", w=5, h=5)

Neighborhood Connectivity:

	VMG <- boxplot.stats(vmg_diff_topo_ssn$NeighborhoodConnectivity)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$NeighborhoodConnectivity)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$NeighborhoodConnectivity)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","a,b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$NeighborhoodConnectivity)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_nc <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=NeighborhoodConnectivity)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Neighborhood Connectivity",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_nc, file="ssn_network_plots/neighborhoodconnectivity_grey.pdf", w=5, h=5)
	
Eccentricity:

	VMG <- boxplot.stats(vmg_diff_topo_ssn$Eccentricity)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$Eccentricity)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$Eccentricity)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$Eccentricity)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_ecc <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=Eccentricity)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Eccentricity",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_ecc, file="ssn_network_plots/eccentricity_grey.pdf", w=5, h=5)		

Clustering Coefficient:

	VMG <- boxplot.stats(vmg_diff_topo_ssn$ClusteringCoefficient)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$ClusteringCoefficient)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$ClusteringCoefficient)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("a","a","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$ClusteringCoefficient)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_clust <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=ClusteringCoefficient)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Clustering Coefficient",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.text.x = element_text(face="bold"), axis.ticks = element_blank())
	ggsave(plot_clust, file="ssn_network_plots/clusteringcoefficient_grey.pdf", w=5, h=5)	

Plot all the above on a 9 panel plot:

	VMG <- boxplot.stats(vmg_diff_topo_ssn$EdgeCount)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$EdgeCount)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$EdgeCount)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("a","a","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$EdgeCount)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_degree <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=EdgeCount)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Total Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.ticks = element_blank())
	
	VMG <- boxplot.stats(vmg_diff_topo_ssn$Indegree)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$Indegree)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$Indegree)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","a,b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$Indegree)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_indegree <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=Indegree)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "In-Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.ticks = element_blank())
	
	VMG <- boxplot.stats(vmg_diff_topo_ssn$Outdegree)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$Outdegree)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$Outdegree)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("a","a","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$Outdegree)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_outdegree <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=Outdegree)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Out-Degree",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.ticks = element_blank())
		
	VMG <- boxplot.stats(vmg_diff_topo_ssn$AverageShortestPathLength)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$AverageShortestPathLength)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$AverageShortestPathLength)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","a,b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$AverageShortestPathLength)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_avgpath <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=AverageShortestPathLength)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Average Shortest Path Length",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.ticks = element_blank())


	VMG <- boxplot.stats(vmg_diff_topo_ssn$BetweennessCentrality)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$BetweennessCentrality)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$BetweennessCentrality)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","a,b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$BetweennessCentrality)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_bc <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=BetweennessCentrality)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Betweenness Centrality",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.ticks = element_blank())


	VMG <- boxplot.stats(vmg_diff_topo_ssn$ClosenessCentrality)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$ClosenessCentrality)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$ClosenessCentrality)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("a","a","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$ClosenessCentrality)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_cc <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=ClosenessCentrality)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Closeness Centrality",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.ticks = element_blank())


	VMG <- boxplot.stats(vmg_diff_topo_ssn$NeighborhoodConnectivity)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$NeighborhoodConnectivity)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$NeighborhoodConnectivity)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","a,b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$NeighborhoodConnectivity)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_nc <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=NeighborhoodConnectivity)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Neighborhood Connectivity",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.ticks = element_blank())


	VMG <- boxplot.stats(vmg_diff_topo_ssn$Eccentricity)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$Eccentricity)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$Eccentricity)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("b","b","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$Eccentricity)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_ecc <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=Eccentricity)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Eccentricity",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.ticks = element_blank())


	VMG <- boxplot.stats(vmg_diff_topo_ssn$ClusteringCoefficient)$stats
	BMG <- boxplot.stats(bmg_diff_topo_ssn$ClusteringCoefficient)$stats
	Non <- boxplot.stats(non_diff_topo_ssn$ClusteringCoefficient)$stats
	merge1 <- data.frame(VMG,BMG,Non)
	merge2 <- t(merge1)
	merge3 <- data.frame(merge2)
	merge3$Dataset <- rownames(merge3)
	merge3$Sign <- c("a","a","a")
	merge3$pos <- (merge3$X5 * 1.05)
	
	sts <- boxplot.stats(vmg_bmg_non_ssn$ClusteringCoefficient)$stats
	order <- c( "Non", "VMG", "BMG")
	#get_n <- function(x){
	#return(data.frame(y = median(x)*1.25, label = paste0("n = ",length(x))))
	#}

	plot_clust <- ggplot(vmg_bmg_non_ssn, aes(x=Dataset, y=ClusteringCoefficient)) +
	geom_boxplot(outlier.colour = NA, aes(fill=Dataset)) +
	scale_fill_grey(start=0.35,end=0.95) +
	labs(x="", y = "Clustering Coefficient",colour="legend" ) +
	guides(fill=FALSE) +
	theme_bw() +
	#stat_summary(fun.data = get_n, geom = "text") +
	geom_text(data = merge3, aes(x = Dataset, y = pos, label = Sign)) +
	coord_cartesian(ylim = c(sts[2]/2,max(merge1)*1.15)) +
	scale_x_discrete(limits=order, labels=c("Non-Differential\n(n=1095)","Differential\nViral MG\n(n=133)","Differential\nTotal MG\n(n=46)")) +
	theme(axis.ticks = element_blank())

	pdf(file="ssn_network_plots/network_topology_panel_grey.pdf",w=12,h=12)
	grid.arrange(plot_degree,plot_indegree,plot_outdegree,plot_bc,plot_cc,plot_clust,plot_avgpath,plot_nc,plot_ecc,legend,ncol=3,nrow = 3)
	dev.off()
	quit()

##Differential KOs mapped to PWs

First, look at how many differential KOs are overlapping between BMG and VMG:

	R-3.2.0/bin/./R
	vmg_lefse <- read.table("viral_blast2tsv/VMG.ko_phast_lefse_output.txt", sep="\t", header=FALSE)
	vmg_lefse$split <- (vmg_lefse$V4 > 1)
	vmg_diff <- split(vmg_lefse, vmg_lefse$split)
	vmg_diff <- as.data.frame(vmg_diff)
	names(vmg_diff)[1] <- "canonicalName"
	names(vmg_diff)[2] <- "KruskallWallis"
	names(vmg_diff)[3] <- "Diet"
	names(vmg_diff)[4] <- "LDA"
	names(vmg_diff)[5] <- "PValue"
	names(vmg_diff)[6] <- "Dataset"
	dim(vmg_diff)

289 diff VMG KOs

	bmg_lefse <- read.table("bmg_blast2tsv/BMG.ko_corrected_abundance_lefse_output.txt", sep="\t", header=FALSE)
	bmg_lefse$split <- (bmg_lefse$V4 > 1)
	bmg_diff <- split(bmg_lefse, bmg_lefse$split)
	bmg_diff <- as.data.frame(bmg_diff)
	names(bmg_diff)[1] <- "canonicalName"
	names(bmg_diff)[2] <- "KruskallWallis"
	names(bmg_diff)[3] <- "Diet"
	names(bmg_diff)[4] <- "LDA"
	names(bmg_diff)[5] <- "PValue"
	names(bmg_diff)[6] <- "Dataset"
	dim(bmg_diff)

118 diff BMG KOs

Intersection:

	overlap <- intersect(bmg_diff$canonicalName, vmg_diff$canonicalName)
	length(overlap)

24 KOs overlap (out of possible 118, ~20.3%)	

What about those specifically involved in metabolism as defined by KEGG:

	ref_topo <- read.table("refnet_attributes.txt", header=TRUE, sep="\t")
	vmg_diff_topo <- merge(x=vmg_diff, y=ref_topo, by="canonicalName")
	bmg_diff_topo <- merge(x=bmg_diff, y=ref_topo, by="canonicalName")
	dim(vmg_diff_topo)
[1] 133  24 - 133 viral KO in metabolism that are diff

	dim(bmg_diff_topo)
[1] 47 24 - 47 BMG KO in metabolism that are diff

	overlap_metab <- intersect(vmg_diff_topo$canonicalName, bmg_diff_topo$canonicalName)
	length(overlap_metab)
[1] 9 - 9 out of a possible 47 overlap

Map KOs to PWs (once again, I can't provide this mapping file because of KEGG liscence):

	write.table(vmg_diff$canonicalName, "viral_diff_ko.txt", row.names=FALSE, col.names=FALSE, quote = FALSE)
	write.table(bmg_diff$canonicalName, "bmg_diff_ko.txt",row.names=FALSE, col.names=FALSE, quote=FALSE)
	write.table(vmg_diff_topo$canonicalName, "viral_diff_metab_ko.txt", row.names=FALSE, col.names=FALSE, quote = FALSE)
	write.table(bmg_diff_topo$canonicalName, "bmg_diff_metab_ko.txt", row.names=FALSE, col.names=FALSE, quote = FALSE)
	q()

Map All Viral KOs:
	
	./ko_to_pw.pl -k_numbers=viral_diff_ko.txt -pw_map=kegg_2/genes/ko/ko_pathway.list
	mv ko_pw.txt viral_ko_pw.txt
	mv pw_list.txt viral_pw_list.txt
	wc -l viral_ko_pw.txt
173 viral_ko_pw.txt - 173 (of289) KOs mapped to PWs
	
	wc -l viral_pw_list.txt 
385 viral_pw_list.txt - 385 total PWs (KOs can map to >1 PW)

Map Metabolic Viral KOs:

	./ko_to_pw.pl -k_numbers=viral_diff_metab_ko.txt -pw_map=kegg_2/genes/ko/ko_pathway.list
	mv ko_pw.txt viral_metab_ko_pw.txt
	mv pw_list.txt viral_metab_pw_list.txt
	wc -l viral_metab_ko_pw.txt
133 viral_metab_ko_pw.txt

	wc -l viral_metab_pw_list.txt
332 viral_metab_pw_list.txt

Map BMG KOs:
	
	./ko_to_pw.pl -k_numbers=bmg_diff_ko.txt -pw_map=kegg_2/genes/ko/ko_pathway.list
	mv ko_pw.txt bmg_ko_pw.txt
	mv pw_list.txt bmg_pw_list.txt
	wc -l bmg_ko_pw.txt 
68 bmg_ko_pw.txt - 68 (of 118) KOs mapped to PWs

	wc -l bmg_pw_list.txt 
153 bmg_pw_list.txt - 153 total PWs

Map Metabolic BMG KOs:

	./ko_to_pw.pl -k_numbers=bmg_diff_metab_ko.txt -pw_map=kegg_2/genes/ko/ko_pathway.list
	mv ko_pw.txt bmg_metab_ko_pw.txt
	mv pw_list.txt bmg_metab_pw_list.txt
	wc -l bmg_metab_ko_pw.txt
47 bmg_metab_ko_pw.txt
      
	wc -l bmg_metab_pw_list.txt
126 bmg_metab_pw_list.txt

Overlap:
	
	R-3.2.0/bin/./R
	viral_pw <- read.table("viral_pw_list.txt")
	bmg_pw <- read.table("bmg_pw_list.txt")
	pw_overlap <- intersect(bmg_pw$V1, viral_pw$V1)
	length(pw_overlap)
[1] 53 - 53 of 153 possible PWs overlap between BMG and VMG (~34.6%)

	viral_metab_pw <- read.table("viral_metab_pw_list.txt")
	bmg_metab_pw <- read.table("bmg_metab_pw_list.txt")
	metab_pw_overlap <- intersect(viral_metab_pw$V1, bmg_metab_pw$V1)
	length(metab_pw_overlap)
[1] 45 - 45 of 126 metabolic pathways overlap (~35.7%)


##Community Stats on all KEGG functions, ORF abundance, metabolic genes


##Viral-Bacterial metagenome read-sharing


##16S alpha and beta-diverstiy





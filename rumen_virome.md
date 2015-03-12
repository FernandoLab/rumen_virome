Impact of dietary change on auxillary viral metabolism and viral-bacterial interactions in the rumen of beef cattle
===============
Author: Christopher L. Anderson (canderson30@unl.edu)


##Introduction
To recreate the analysis from Anderson et al. manuscript.  All commands below were done in a linux enviornment and memory intensive commands were carried out on a high performance cluster at UNL.  Requirements noted so far, but are most likely available to you (other versions likely work):

- java (version 1.8 used for manuscript)
- perl (version?)
- python (version 2.7)


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
    scp canderson3@crane.unl.edu:/work/samodha/canderson3/raw_viral/raw_viral_fastq.tgz raw_viral/
    scp canderson3@crane.unl.edu:/work/samodha/canderson3/raw_total/raw_total_fastq.tgz raw_total/

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
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.1.fastq trimmomatic_output/VMG.1_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:330 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.2.fastq trimmomatic_output/VMG.2_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:200 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.3.fastq trimmomatic_output/VMG.3_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:330 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.4.fastq trimmomatic_output/VMG.4_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:330 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.5.fastq trimmomatic_output/VMG.5_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:150 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.6.fastq trimmomatic_output/VMG.6_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:285 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.7.fastq trimmomatic_output/VMG.7_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:350 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.8.fastq trimmomatic_output/VMG.8_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:330 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.9.fastq trimmomatic_output/VMG.9_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:330 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.10.fastq trimmomatic_output/VMG.10_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:230 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.11.fastq trimmomatic_output/VMG.11_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:330 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.12.fastq trimmomatic_output/VMG.12_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:350 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.13.fastq trimmomatic_output/VMG.13_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:330 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.14.fastq trimmomatic_output/VMG.14_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:230 MINLEN:85
    java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/VMG.15.fastq trimmomatic_output/VMG.15_trimm.fastq ILLUMINACLIP:Trimmomatic-0.33/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:330 MINLEN:85

Alternatively, the outputs from the above command are available at interm/trimmomatic_output.

To remove artifical duplicates we use cd-hit-454.  We noticed early on the transposon library preps with the Ion Torrent Museek kits introduce many duplicates (associated with viral sample?). While originally designed for 454 sequencing, Ion Torrent has shares many of the same error profiles as 454 and cd-hit-454 appears to work fairly well for removing duplicates in our data.

First download cd-hit-454:

    wget https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz
    tar xvf cd-hit-v4.6.1-2012-08-27.tgz
    cd cd-hit-v4.6.1-2012-08-27
    make openmp=yes
    cd ..

Then remove the duplicates from all viral metagenome samples and write output to cd_hit_454_output

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

###Total Metagenome:
TO_DO




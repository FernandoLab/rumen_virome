Impact of dietary change on auxillary viral metabolism and viral-bacterial interactions in the rumen of beef cattle
===============
Author: Christopher L. Anderson


##Introduction
To recreate the analysis from Anderson et al. manuscript.  All commands below were done in a linux enviornment and memory intensive commands were carried out on a high performance cluster at UNL.  Requirements noted so far, but are most likely available to you (other versions likely work):

* java (version 1.8 used for manuscript)
* perl (version?)
* python (version 2.7)

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
In order to skip some long computational steps, you can use the outputs in this directory (too large to just store on github?)

    scp -r canderson3@crane.unl.edu:/work/samodha/canderson3/interm ./

##QC
To ensure the data was properly trimmed by the Torrent Server and to deal with known biases such as artifical duplications created by the transposon library preps used for the viral metagenomes, we performed additional QC steps.

To ensure adaptor trimming, removal of barcode seqeunce, and removal of inserted transposon sequences we used Trimmomatic.  To get the software:

    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
    unzip Trimmomatic-0.33.zip 

The file transp_adapt_remove.cat.txt (in intermediate file directory, iterm/transp_adapt_remove.cat.txt) was supplied to direct removal of adaptors, barcodes, and transposon seqeunces. Further, observed GC and k-mer bias was noted in the 5' and 3' ends.  In order to begin to alleviate those biases, the first 20 bp of the read was trimmed and a sample dependent trimming from the 3' end was done (more 3' trimming later).  To trim all samples:

    mkdir trimmomatic_output
    cd raw_viral/
    for f in *.fastq
    do
        filename=$(basename "$f")
        cd /work/samodha/canderson3
        java -jar Trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 raw_viral/$f trimmomatic_output/"$filename""_trimmomatic.fastq" ILLUMINACLIP:interm/transp_adapt_remove.cat.txt:2:20:10 HEADCROP:20 CROP:300 MINLEN:85
    done

Alternatively, the output from the above command is available at interm/trimmomatic_output.






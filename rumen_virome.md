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
    scp canderson3@crane.unl.edu:/work/samodha/canderson3/raw_viral.tgz raw_viral/
    scp canderson3@crane.unl.edu:/work/samodha/canderson3/raw_total.tgz raw_total/
    tar -zxvf raw_viral.tgz
    tar -zxvf raw_total.tgz
    #set home variable

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


Get khmer software to compare kmer profiles of different metagenomes:

    cd virtualenv-*; python virtualenv.py ../khmerEnv; cd ..
    source khmerEnv/bin/activate
    pip install khmer

Now, load all total metagenomes into counting:

    source khmerEnv/bin/activate
    mkdir /work/samodha/canderson3/khmer_counting
    cd /work/samodha/canderson3/cd_hit_454_output_total
    for f in *.fastq
    do
        filename=$(basename "$f")
        filename="${filename%_*}"
        python /work/samodha/canderson3/khmer-1.3/scripts/load-into-counting.py -k 20 -N 4 -x 8e9 /work/samodha/canderson3/khmer_counting/"$filename""_k20.ct" $f --report-total-kmers -s tsv
    done



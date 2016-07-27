# Setup an Environment and Install Dependencies
# This first part of the anlaysis is to install the dependency 
# software and associated packages in the hopes of creating a reproducible 
# environment to work from over time. 

# ensure pwd is the cloned repository
result=${PWD##*/}
if [ "$result" != "rumen_virome" ]
then
	printf "\nCurrent directory is not the cloned repository.\nSee https://github.com/chrisLanderson/rumen_virome for details.\n\n"
	exit 1
fi

# ensure provided link to usearch download
if [ "$1" = "" ]; then
    printf "\nProvide a link for USEARCH download (from email) as argument.\nGet a license from http://www.drive5.com/usearch/download.html\nSee https://github.com/chrisLanderson/rumen_virome for details.\n\n"
    exit 1
fi

# trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
unzip Trimmomatic-0.33.zip 

# cd-hit
wget https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz
tar -xvf cd-hit-v4.6.1-2012-08-27.tgz
cd cd-hit-v4.6.1-2012-08-27
make openmp=yes
cd ..

# prinseq-lite
wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz
tar -xvf prinseq-lite-0.20.4.tar.gz

# anaconda 
wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda-2.3.0-Linux-x86_64.sh
bash Anaconda-2.3.0-Linux-x86_64.sh
anaconda/bin/conda create -n rumenVirome python=2.7 qiime=1.9.1
source anaconda/bin/activate rumenVirome
#cp pinned anaconda/conda-meta

# r
conda install -c r rpy2=2.5.6 r-devtools=1.9.1 r-curl=0.9.4
conda install -c r r=3.2.2

# cutadapt
pip install cutadapt==1.8.1

# khmer
pip install khmer==1.4.1

# biopython
pip install biopython

# mothur
wget https://github.com/mothur/mothur/releases/download/v1.35.1/Mothur.cen_64.zip
unzip Mothur.cen_64.zip

# usearch
wget -O anaconda/envs/rumenVirome/bin/usearch8.0.1623 $1
chmod 775 anaconda/envs/rumenVirome/bin/usearch8.0.1623

# rRNA prediction
wget http://weizhong-lab.ucsd.edu/meta_rna/rRNA_prediction.tar.bz2
bzip2 -d rRNA_prediction.tar.bz2
tar -xvf rRNA_prediction.tar

# pandoc
conda install -c https://conda.binstar.org/asmeurer pandoc

# spades
wget http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0-Linux.tar.gz
tar -zxvf SPAdes-3.5.0-Linux.tar.gz

# prodigal
wget https://github.com/hyattpd/Prodigal/archive/v2.6.2.zip
unzip v2.6.2.zip
cur=$(pwd)
cd Prodigal-2.6.2/
make install INSTALLDIR=$cur
cd ..

# bowtie2
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip/download -O bowtie-2.2.5.source.zip
unzip bowtie-2.2.5.source.zip

# sra toolkit
wget ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -xzf sratoolkit.current-centos_linux64.tar.gz

# cytoscape utility 
wget https://raw.githubusercontent.com/idekerlab/cy-rest-R/develop/utility/cytoscape_util.R

# R packages
printf "\nInstallling R packages, will take some time...\n"
R CMD BATCH scripts/install_pack.R

conda install bioconductor-genomicranges=1.20.8
wget https://bioc.ism.ac.jp/packages/3.1/bioc/src/contrib/DESeq2_1.8.2.tar.gz
R CMD INSTALL DESeq2_1.8.2.tar.gz

# unpack intermediate results
for f in intermediate_results/*.zip
do
	unzip $f -d intermediate_results/
done

for f in intermediate_results/*.biom
do
	filename=$(basename "$f")
	filename="${filename%.biom}"
	biom convert --table-type="OTU table" --to-tsv -i $f -o intermediate_results/$filename.txt
done

unzip prophage_virus.db.zip

# clean up
rm intermediate_results/*.zip
mv total_adapt_remove.cat.txt Trimmomatic-0.33/
mv transp_adapt_remove.cat.txt Trimmomatic-0.33/
rm Trimmomatic-0.33.zip
mv cd-hit-v4.6.1-2012-08-27/cd-hit-454 anaconda/envs/rumenVirome/bin/
rm cd-hit-v4.6.1-2012-08-27.tgz
rm -r cd-hit-v4.6.1-2012-08-27
chmod 775 prinseq-lite-0.20.4/prinseq-lite.pl
mv prinseq-lite-0.20.4/prinseq-lite.pl anaconda/envs/rumenVirome/bin/
rm -r prinseq-lite-0.20.4
rm prinseq-lite-0.20.4.tar.gz
mv mothur/mothur anaconda/envs/rumenVirome/bin/
rm v2.6.2.zip
rm -r mothur
rm Mothur.cen_64.zip
rm -r __MACOSX
mv usearch8.0.1623 anaconda/envs/rumenVirome/bin/
rm SPAdes-3.5.0-Linux.tar.gz
chmod 775 prodigal
mv prodigal anaconda/envs/rumenVirome/bin/
rm -r Prodigal-2.6.2
chmod 775 bowtie2-2.2.5/*
mv bowtie2-2.2.5/* anaconda/envs/rumenVirome/bin/
rm bowtie-2.2.5.source.zip
rm -r bowtie2-2.2.5
rm sratoolkit.current-centos_linux64.tar.gz
rm rRNA_prediction.tar
rm Anaconda-2.3.0-Linux-x86_64.sh
rm install_pack.Rout

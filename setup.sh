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
mv Trimmomatic-0.33/trimmomatic-0.33.jar anaconda/envs/rumenVirome/bin/

# cd-hit
wget https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz
tar -xvf cd-hit-v4.6.1-2012-08-27.tgz
cd cd-hit-v4.6.1-2012-08-27
make openmp=yes
cd ..
mv cd-hit-v4.6.1-2012-08-27/cd-hit-454 anaconda/envs/rumenVirome/bin/

# prinseq-lite
wget http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz
tar -xvf prinseq-lite-0.20.4.tar.gz
mv prinseq-lite-0.20.4/prinseq-lite.pl anaconda/envs/rumenVirome/bin/

# anaconda 
wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda-2.3.0-Linux-x86_64.sh
bash Anaconda-2.3.0-Linux-x86_64.sh
anaconda/bin/conda create -n rumenVirome python=2.7 pip numpy=1.9.2 matplotlib=1.4.3 scipy=0.16.0 pandas=0.16.2 cython mock=1.1.3 nose=1.3.7 urllib beautifulsoup4
source anaconda/bin/activate rumenVirome

# R
conda install -c r r=3.2.1
conda install -c https://conda.binstar.org/r rpy2=2.5.6 r-devtools r-rcurl

# QIIME
pip install qiime==1.9.1

# cutadapt
pip install cutadapt==1.8.1

# khmer
pip install khmer==1.4.1

# mothur
wget https://github.com/mothur/mothur/releases/download/v1.35.1/Mothur.cen_64.zip
unzip Mothur.cen_64.zip
mv mothur/mothur anaconda/envs/rumenVirome/bin/

# usearch
wget $1 -O usearch8.0.1623
chmod 775 usearch8.0.1623
mv usearch8.0.1623 anaconda/envs/rumenVirome/bin/

# rRNA prediction
wget http://weizhong-lab.ucsd.edu/meta_rna/rRNA_prediction.tar.bz2
bzip2 -d rRNA_prediction.tar.bz2
tar -xvf rRNA_prediction.tar
chmod 777 -R rRNA_prediction/

# pandoc
conda install -c https://conda.anaconda.org/asmeurer pandoc

# spades
wget http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0-Linux.tar.gz
tar -zxvf SPAdes-3.5.0-Linux.tar.gz
mv SPAdes-3.5.0-Linux/bin/* anaconda/envs/rumenVirome/bin/

# prodigal
wget https://github.com/hyattpd/Prodigal/archive/v2.6.2.zip
unzip v2.6.2.zip
cur=$(pwd)
cd Prodigal-2.6.2/
make install INSTALLDIR=$cur
cd ..
mv prodigal anaconda/envs/rumenVirome/bin/

# bowtie2
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip/download -O bowtie-2.2.5.source.zip
unzip bowtie-2.2.5.source.zip
mv bowtie2-2.2.5/bowtie2 anaconda/envs/rumenVirome/bin/
mv bowtie2-2.2.5/bowtie2-build anaconda/envs/rumenVirome/bin/

# R packages
printf "\nInstallling R packages, will take some time...\n"
R CMD BATCH scripts/install_pack.R

# clean up
rm cd-hit-v4.6.1-2012-08-27.tgz
rm -rf cd-hit-v4.6.1-2012-08-27
rm prinseq-lite-0.20.4.tar.gz
rm -rf prinseq-lite-0.20.4
rm Mothur.cen_64.zip
rm -rf __MACOSX
rm Trimmomatic-0.33.zip
rm rRNA_prediction.tar
rm Anaconda-2.3.0-Linux-x86_64.sh
rm SPAdes-3.5.0-Linux.tar.gz
rm -rf SPAdes-3.5.0-Linux
rm bowtie-2.2.5.source.zip
rm -rf bowtie2-2.2.5
rm v2.6.2.zip
rm -rf mothur
rm -rf Trimmomatic-0.33
rm -rf Prodigal-2.6.2
rm install_pack.Rout
mv transp_adapt_remove.cat.txt anaconda/envs/rumenVirome/bin/
mv total_adapt_remove.cat.txt anaconda/envs/rumenVirome/bin/
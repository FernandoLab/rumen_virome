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

# anaconda, R, QIIME
wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda-2.3.0-Linux-x86_64.sh
bash Anaconda-2.3.0-Linux-x86_64.sh
anaconda/bin/conda create -n rumenVirome python=2.7 pip numpy=1.9.2 matplotlib=1.4.3 scipy=0.16.0 pandas=0.16.2 cython mock=1.1.3 nose=1.3.7
source anaconda/bin/activate rumenVirome
conda install -c r r=3.2.1
conda install -c https://conda.binstar.org/r rpy2=2.5.6
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
wget -O anaconda/envs/rumenEnv/bin/usearch8.0.1623 $1
chmod 775 anaconda/envs/rumenEnv/bin/usearch8.0.1623

# rRNA prediction
wget http://weizhong-lab.ucsd.edu/meta_rna/rRNA_prediction.tar.bz2
bzip2 -d rRNA_prediction.tar.bz2
tar -xvf rRNA_prediction.tar
chmod 777 -R rRNA_prediction/

# pandoc
conda install -c https://conda.anaconda.org/asmeurer pandoc

# R packages
printf "\nInstallling R packages, will take some time...\n"
R CMD BATCH scripts/install_pack.R

# Get Data
mkdir raw_viral
mkdir raw_total
scp canderson3@crane.unl.edu:/work/samodha/canderson3/raw_viral.tgz raw_viral/
scp canderson3@crane.unl.edu:/work/samodha/canderson3/raw_total.tgz raw_total/
#get illumina viral data and unzip
tar -zxvf raw_viral.tgz
tar -zxvf raw_total.tgz

# remove
rm cd-hit-v4.6.1-2012-08-27.tgz
rm prinseq-lite-0.20.4.tar.gz
rm Mothur.cen_64.zip
rm -rf __MACOSX
rm Trimmomatic-0.33.zip
rm Anaconda-2.3.0-Linux-x86_64.sh

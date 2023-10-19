# Install Miniconda3
mkdir -p software && cd software
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-py39_23.5.2-0-Linux-x86_64.sh --no-check-certificate


bash Miniconda3-py39_23.5.2-0-Linux-x86_64.sh
## Follow the prompts, press Enter as needed, accept the license agreement, and input 'yes' at the end. Confirm the installation directory and initialize conda with 'yes' as well.

which conda # Check if the installation was successful

# Create an environment for metagenomics analysis
conda create -n "MetaG"
## Add a domestic source
## Open ~/.condarc file, create one if it doesn't exist
touch ~/.condarc
# Open .condarc and paste the following code (excluding single quotes)
'
channels:
#  - bioconda
#  - conda-forge
#  - defaults
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/


auto_activate_base: false
show_channel_urls: true
'


# Install fastp using conda
# https://anaconda.org/anaconda/conda
conda activate MetaG  # Activate the environment
conda install fastp   # Install fastp, refer to the official website https://github.com/OpenGene/fastp
## Note that you need to input 'yes' to confirm or simply press Enter

# Install fastQC
## Typically, fastqc is downloaded from the official website
## But you can try conda
## Let's search to see if there's a suitable version
conda search fastqc
## Found a newer version, install it using conda
conda install fastqc

# Install prodigal mmseqs2 diamond minimap2 kraken2 kaiju humann3 samtools megahit bowtie2 seqkit multiqc axel using conda
conda install prodigal  mmseqs2 diamond minimap2 kraken2 bracken kaiju  samtools megahit  bowtie2 seqkit multiqc axel sortmerna taxonkit vmtouch coverm pysqlite3 pandas


# Manually install bowtie2
## If there are conda conflicts that cannot be resolved, consider downloading the binary file (precompiled, ready to use)
cd ~/software
wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.1/bowtie2-2.5.1-linux-x86_64.zip
unzip bowtie2-2.5.1-linux-x86_64.zip
cd bowtie2-2.5.1-linux-x86_64 
pwd  # Get the path, it should be /PERSONALBIO196/LN00/software/bowtie2-2.5.1-linux-x86_64
## Add this path to ~/.bashrc
echo 'PATH=/PERSONALBIO196/LN00/software/bowtie2-2.5.1-linux-x86_64:$PATH' >> ~/.bashrc
source ~/.bashrc
bowtie2 --version
## If you encounter an error: /usr/bin/env: perl: No such file or directory
## Install perl
conda install perl
## No error means no need to handle it

# Install featureCounts
wget https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz/download
mv download subread-2.0.6-Linux-x86_64.tar.gz
tar -xzvf subread-2.0.6-Linux-x86_64.tar.gz

# Manually reinstall mmseqs2, not usually required but necessary on an older server
cd ~/software
axel -n 10 https://mmseqs.com/latest/mmseqs-linux-sse2.tar.gz
tar xvfz mmseqs-linux-sse2.tar.gz
echo 'PATH=/PERSONALBIO196/LN00/software/mmseqs/bin:$PATH' >> ~/.bashrc

 ~/.bashrc

## Create a new environment and install R 4.3
conda deactivate
conda create -n R4.3 R==4.3
conda activate R4.3
conda install r-data.table r-magrittr

# Create another environment for Humann3
conda deactivate
conda create -n humann3 humann==3.8

# Install awscli, follow the official documentation for installation
## https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html
cd ~/software
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
mkdir -p ~/bin
./aws/install -i $(pwd)/aws-cli -b ${HOME}/bin --update

# Download the kaiju database
cd ~; mkdir -p database/kaijudb;cd database/kaijudb
axel -n 10 https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_nr_euk_2023-05-10.tgz
# If you encounter issues downloading, it's because the resource is on AWS S3
# Check the contents first
aws s3 ls --no-sign-request s3://kaiju-idx/
# Download it
aws s3 cp --no-sign-request s3://kaiju-idx/2023/kaiju_db_refseq_2023-05-23.tgz .

# Download kraken2 database
## Use AWS CLI or axel to download it
## https://benlangmead.github.io/aws-indexes/k2

#axel -n https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz

#awscli
aws s3 sync --no-sign-request s3://genome-idx/kraken/k2_standard_20230605.tar.gz .


# Download sormerna rRNA ref
# https://github.com/sortmerna/sortmerna/tree/master/data/rRNA_databases
cd ~/database
mkdir -p sortmernaDB;cd sortmernaDB
wget https://github.com/sortmerna/sortmerna/archive/refs/heads/master.zip 
unzip master.zip
mv sortmerna-master/data/rRNA_databases/*.fasta .
rm master.zip sortmerna-master/ -r






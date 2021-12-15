#!/bin/sh

#  tutorial.sh
#  
#
<<<<<<< Updated upstream
#  Created by Kim Dill-McFarland on 12/14/21.
#  


nohup snakemake --snakefile Snakefile_step1 --cores 15 > log/SEAsnake_step1.log 2>&1
=======
#  Created by Kim Dill-McFarland on 10/14/21.
#

#AWS snakemake install notes

ssh -i ~/Documents/AWS/keys/SEAsnake_key.cer ec2-user@ec2-44-242-145-12.us-west-2.compute.amazonaws.com

#### Basic AWS ####
sudo yum upgrade -y
sudo yum update -y

#### Setup EBS volume ####
## Get addtl volume name
ebs_name=$(lsblk -o NAME -n -i -r | tail -n 1)
## Format volume
sudo mkfs -t ext4 /dev/$ebs_name
## Attach project directory to volume
sudo mkdir -p ~/project
sudo mount /dev/$ebs_name ~/project/
## Change permissions to read-write
sudo chmod 777 -R ~/project/

#### Install fuse ####
sudo amazon-linux-extras install -y epel
sudo yum install -y s3fs-fuse

#### Install conda ####
## Download conda
sudo mkdir -p ~/apps/anaconda
sudo chmod 777 -R ~/apps
cd ~/apps/anaconda
sudo curl -O https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh

## Compile and install conda
sudo bash Anaconda3-2021.11-Linux-x86_64.sh -b -p /home/ec2-user/apps/anaconda -u
eval "$(/home/ec2-user/apps/anaconda/bin/conda shell.bash hook)"
conda init
sudo chmod 777 -R ~/apps
## Restart instance

## Get addtl conda channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set allow_conda_downgrades true

## Install mamba
conda install -n base -c conda-forge mamba -y

## Create environment and install SEAsnake software with mamba
mamba env create --name SEAsnake --file ~/SEAsnake/environment/Hissss_env.yaml

#### Install SEAsnake ####
## Install git
sudo yum install git -y

## Clone SEAsnake
sudo mkdir ~/project/SEAsnake
sudo chmod 777 ~/project/SEAsnake
git clone https://github.com/BIGslu/SEAsnake ~/project/SEAsnake

####################################

## Configure your account
export AWS_ACCESS_KEY_ID=<AWS_ACCESS_KEY>
export AWS_SECRET_ACCESS_KEY=<AWS_SECRET_ACCESS_KEY>
export AWS_DEFAULT_REGION=us-west-2

## FILL IN WITH YOUR KEYS ##
echo AWS_ACCESS_KEY:AWS_SECRET_ACCESS_KEY > ~/.passwd-s3fs
chmod 600 ~/.passwd-s3fs

### Mount data from S3 using fuse ###

s3fs YOUR_DATA_BUCKET ~/project/SEAsnake/data \
    -o passwd_file=~/.passwd-s3fs \
    -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007

## If genome index already available on S3
## Index must be in exact path on S3
## release###/STARindex/
sudo mkdir ~/project/SEAsnake/ref
sudo chmod 777 ~/project/SEAsnake/ref
s3fs human-ref ~/project/SEAsnake/ref \
    -o passwd_file=~/.passwd-s3fs \
    -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007

#### Run SEAsnake ####
## Activate environment with software installed
conda activate SEAsnake

## Run step 1 to create config file and perform initial sequence quality assessment
cd ~/project/SEAsnake/
git checkout 2c84225536e402887a6c5aea10d21f729c5d2737
rm -R result
snakemake --snakefile Snakefile_step1 --cores 8
nohup snakemake --snakefile Snakefile_step1 --cores 8 >> log/SEAsnake_step1.log 2>&1

## Modify config file as needed

## Run step 2 to filter, align, assess, and count sequences
nohup snakemake --snakefile Snakefile_step2 --cores 15 >> log/SEAsnake_step2.log 2>&1

## Exit environment
conda deactivate


git reset --hard origin/main
git pull
fusermount -u SEAsnake/ref

## Update SEAsnake
cd ~/SEAsnake
git pull
>>>>>>> Stashed changes

#!/bin/bash

: <<'COMMENT'
Install conda/mamba if not yet installed following this instruction.

wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

bash Miniforge3-Linux-x86_64.sh # this will create a directory miniforge3

COMMENT


# create a directory to store your executable files
mkdir -p ~/tools/bin

# create & activate conda environment for the analysis
# You can chose to use conda or mamba.

# Check if mamba is installed
if command -v mamba &> /dev/null;
then
        # create conda environment and install necessary tools including java version 11 or higher
        echo "Installing tools using mamba"
        mamba create \
                -n gatk \
                -y \
                -c bioconda \
                -c conda-forge \
                openjdk \
                python=3.12 \
                fastqc=0.11 \
                fastp=0.23 \
                bwa=0.7 \
                samtools=1.21 \
                bcftools=1.21 
                
        mamba activate gatk        

elif command -v conda &> /dev/null;
then
        echo "Installing tools using conda"
        conda create \
                -n gatk \
                -y \
                -c bioconda \
                -c conda-forge \
                openjdk \
                python=3.12 \
                fastqc=0.11 \
                fastp=0.23 \
                bwa=0.7 \
                samtools=1.21 \
                bcftools=1.21 
                
        conda activate gatk        
else
        echo "OOPs..! miniconda/anaconda/micromamba not installed. Please install miniconda using the following commands and try again!"
        echo '
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && bash Miniforge3-Linux-x86_64.sh'
        
        exit 1
fi


# Download and unzip gatk-4.2.2.0
wget https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip

unzip gatk-4.2.2.0.zip
cp gatk-4.2.2.0/gatk ~/tools/bin
export PATH=$PATH:~/tools/bin # add gatk to your PATH

# Download and unzip snpEff for variant annotation
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

unzip snpEff_latest_core.zip

# If everything in this file is set successfully, go ahead and edit the 'OncoVarDetect.slurm.sh' file and proceed to run
# Ensure both 'OncoVarDetect.slurm.sh' and 'OncoVarDetect.sh' files are in the same directory.

# <======================================================SUCCESS========================================================>

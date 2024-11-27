# Please read this carefully and ensure your environment is fully set then proceed to the submission script. 

# Install conda/mamba if not yet installed
wget https://github.com/conda-\
forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

bash Miniforge3-Linux-x86_64.sh # this will create a directory miniforge3

# create a directory to store your executable files
mkdir -p ~/tools/bin
cp -r miniforge3/bin/* ~/tools/bin
export PATH=$PATH:~/tools/bin


# create & activate conda environment for the analysis
# You can chose to use conda or mamba.
# with conda
conda create -n gatk
conda activate gatk

# OR
mamba create -n gatk
mamba activate gatk

# Install java version 11 or higher in the gatk environment
conda install -c conda-forge openjdk
# OR
mamba install -c conda-forge openjdk
java -version  # check the version of installed java

# Run this command to install the following tools
conda install -y \
       -c bioconda \
       -c conda-forge \
       fastqc=0.11 fastp=0.23 bwa=0.7 samtools=1.21 bcftools=1.21 

# OR
mamba install -y \
       -c bioconda \
       -c conda-forge \
       fastqc=0.11 fastp=0.23 bwa=0.7 samtools=1.21 bcftools=1.21 

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
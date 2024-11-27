#!/bin/bash

# ONLY edit the SUBMISSION and PATH BLOCKS of this script!!!!

# ------------
## SUBMISSION BLOCK:
# ------------
# NOTE; DO NOT UNCOMMENT ANY LINE IN THIS BLOCK.
## edit this section to match your HPC system and your memory availability.

#SBATCH --job-name=OncoVarDetect
#SBATCH --output=var_output.txt
#SBATCH --error=var_errors.txt
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=16

# -----------
## PATH BLOCK:
# -----------
# ==> Provide the path to the following to fit your workspace.

# 1. raw reads directory
sampleDir=/etc/ace-data/CancerGenomicsWG/users/walter/GATK/rptGatk/READS

# 2. reference directory
refDir=/etc/ace-data/CancerGenomicsWG/users/walter/GATK/rptGatk/ref

# 3. output directory
outDir=/etc/ace-data/CancerGenomicsWG/users/walter/GATK/rptGatk/Results    

# 4. tab-delimited file with sample IDs
sampleNamesFile=/etc/ace-data/CancerGenomicsWG/users/walter/GATK/rptGatk/sample.tab

# 5. snpEff.jar file; for annotation
snpEff=/etc/ace-data/CancerGenomicsWG/users/walter/GATK/rptGatk/snpEff/snpEff.jar


PL=ILLUMINA      # edit this if your reads are not sequenced using ILLUMINA
env=gatk         # edit to match your analysis conda env't
exe=~/bin        # edit to match the directory with your executables
readType=fastq   # edit this to match your reads' extension [eg., <fastq>, <fq>]

# ------------------
## DESCRIPTION BLOCK:
# ------------------
# ==> Read the comment below to understand what kind of files you should have before running the script.
________________________________________________________________________________________________________
: << 'comment'
1. <sampleNamesFile>  should look like this;

A tab delimited sampleNamesFile contains sample ids and should look like this;
sample1  tumorID  normalID
sample2  tumorID  normalID
sample3  tumorID  normalID
sample4  tumorID  normalID
...
...
...
sampleN  tumorID  normalID

REPLACE tumorID and normalID with the actual sampleIDs e.g.,
sample1  SRR25434464     SRR25434465
sample2  SRR25434466     SRR25434467

2. <reference directory>  should contain these ref & index files;
you can download from google cloud following this link;
<https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false>

1000G_phase1.snps.high_confidence.hg38.vcf
1000G_phase1.snps.high_confidence.hg38.vcf.idx
af-only-gnomad.hg38.vcf
af-only-gnomad.hg38.vcf.idx
Homo_sapiens_assembly38.dbsnp138.vcf
Homo_sapiens_assembly38.dbsnp138.vcf.idx
Homo_sapiens_assembly38.dict
Homo_sapiens_assembly38.fasta
Homo_sapiens_assembly38.fasta.amb
Homo_sapiens_assembly38.fasta.ann
Homo_sapiens_assembly38.fasta.bwt
Homo_sapiens_assembly38.fasta.fai
Homo_sapiens_assembly38.fasta.pac
Homo_sapiens_assembly38.fasta.sa
Homo_sapiens_assembly38.known_indels.vcf
Homo_sapiens_assembly38.known_indels.vcf.idx
Illumina_Exome_TargetedRegions_v1.2.hg38.bed
Mills_and_1000G_gold_standard.indels.hg38.vcf
Mills_and_1000G_gold_standard.indels.hg38.vcf.idx
Mills_and_1000G_gold_standard.indels.hg38.vcf.tbi

3. <sample directory>
- Your sample directory should contain all your fastq samples (both normal & tumor).
- fastq files come with the following possible file extensions
        > fastq
        > fq
        e.g.,
        > read1.fastq, read2.fastq 
        OR
        > read1.fq,    read2.fq
- Make sure you edit the 'ReadType' above (in the PATH BLOCK) to correctly match the extensions 
  of your reads

comment
# __________________________________________________________________________________________________

# ----------------
## EXECUTION BLOCK:          PLEASE, DO NOT EDIT THIS BLOCK
# ----------------

# Execute your main script and pass the directories and file as arguments
conda activate ${env}
export PATH=$PATH:${exe}
bash OncoVarDetect.sh \
    "${sampleDir}" \
    "${refDir}" \
    "${sampleNamesFile}" \
    "${outDir}" \
    "${snpEff}" \
    "${PL}" \
    "${env}" \
    "${exe}" \
    "${readType}"
    






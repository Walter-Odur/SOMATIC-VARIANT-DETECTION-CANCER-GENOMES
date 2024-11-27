#!/bin/bash

SECONDS=0
# Please, DO NOT edit this script

# Check if the correct number of arguments is provided
if [ "$#" -ne 9 ];
then
            echo "Please Enter the correct number of arguments!"
            echo "Usage: $0 <sampleDir> <refDir> <sampleFile> <outDir> <snpEff.jar> <seqPlatform> <condaEnv> <exeDir> <readType>"
            exit 1
fi

# Get the parameters from the arguments
sampleDir=$1
refDir=$2
sampleNamesFile=$3
outDir=$4
snpEff=$5
PL=$6
env=$7
exe=$8
readType=$9

conda activate ${env}
eport PATH=$PATH:${exe}

# Ensure that sample.tab exists
if [ ! -f "${sampleNamesFile}" ];
then
            echo "${sampleNamesFile} not found! Please provide a file with tumor-normal sample pairs."
            exit 1
fi

# Define reference genome
ref="${refDir}/Homo_sapiens_assembly38.fasta"

# Create necessary directories
mkdir -p ${outDir}/{fastQC,trim,bams,gatk,recal,pileups,mut2,annotations,tempDir}

echo "
sample_Directory: ${sampleDir}
reference_Directory: ${refDir}
sample_IDs: ${sampleNamesFile}
Output_Directory: ${outDir}
Sequencing_Platform: ${PL}
snpEff_Directory: ${snpEff}
conda/run env't: ${env}
reads extension type: ${readType}
"

#---------- Step 1: QUALITY CHECKS ----------
echo "Starting FastQC..."
while read -r sample tumor normal;
do
            echo "Performing FastQC for $tumor and $normal"
            fastqc "${sampleDir}/${tumor}"*".${readType}" -o ${outDir}/fastQC || exit 1
            fastqc "${sampleDir}/${normal}"*".${readType}" -o ${outDir}/fastQC || exit 1

done < "${sampleNamesFile}"
echo "QC Done for all samples"


#---------- Step2: TRIMMING POOR QUALITY READS ----------
echo "Starting trimming..."
while read -r sample tumor normal;
do
            echo "Trimming reads for $tumor and $normal"
            fastp \
                    -i "${sampleDir}/${tumor}"*1.${readType} -I "${sampleDir}/${tumor}"*2.${readType} \
                    -o ${outDir}/trim/${tumor}_1_trim.fastq -O ${outDir}/trim/${tumor}_2_trim.fastq || exit 1
            fastp \
                    -i "${sampleDir}/${normal}"*1.${readType} -I "${sampleDir}/${normal}"*2.${readType} \
                    -o ${outDir}/trim/${normal}_1_trim.fastq -O ${outDir}/trim/${normal}_2_trim.fastq || exit 1
done < "${sampleNamesFile}"
echo "Trimming done"

#---------- Step3: ALIGNMENT TO THE REFERENCE ----------
echo "Starting alignment..."
while read -r sample tumor normal;
do
            echo "Aligning $tumor and $normal to the reference genome"
            bwa mem \
                    -t 16 -aM ${ref} -R "@RG\tID:${tumor}\tSM:${tumor}\tPL:${PL}" \
                    ${outDir}/trim/${tumor}_1_trim.fastq ${outDir}/trim/${tumor}_2_trim.fastq | \
                    samtools view -bS - > ${outDir}/bams/${tumor}.bam || exit 1
            bwa mem \
                    -t 16 -aM ${ref} -R "@RG\tID:${normal}\tSM:${normal}\tPL:${PL}" \
                    ${outDir}/trim/${normal}_1_trim.fastq ${outDir}/trim/${normal}_2_trim.fastq | \
                    samtools view -bS - > ${outDir}/bams/${normal}.bam || exit 1
            samtools sort \
                    -o ${outDir}/bams/${tumor}_sorted.bam ${outDir}/bams/${tumor}.bam && \
                    rm -rf ${outDir}/bams/${tumor}.bam && \
                    mv ${outDir}/bams/${tumor}_sorted.bam ${outDir}/bams/${tumor}.bam || exit 1
            samtools sort \
                    -o ${outDir}/bams/${normal}_sorted.bam ${outDir}/bams/${normal}.bam && \
                    rm -rf ${outDir}/bams/${normal}.bam && \
                    mv ${outDir}/bams/${normal}_sorted.bam ${outDir}/bams/${normal}.bam || exit 1
            samtools index ${outDir}/bams/${tumor}.bam || exit 1
            samtools index ${outDir}/bams/${normal}.bam || exit 1
            

done < "${sampleNamesFile}"
echo "Alignment done"


# #---------- Step4: MARKING DUPLICATES ----------
echo "Marking duplicates..."
while read -r sample tumor normal;
do
            echo "Marking duplicates for $tumor and $normal"
            gatk MarkDuplicates \
                    --INPUT ${outDir}/bams/${tumor}.bam \
                    --OUTPUT ${outDir}/gatk/${tumor}_mkdp.bam \
                    --METRICS_FILE ${outDir}/gatk/${tumor}_metrics.txt \
                    --REMOVE_DUPLICATES false \
                    --CREATE_INDEX true || exit 1
            samtools flagstat ${outDir}/gatk/${tumor}_mkdp.bam > ${outDir}/gatk/${tumor}_mkdp.bam.flagstat || exit 1

            gatk MarkDuplicates \
                    --INPUT ${outDir}/bams/${normal}.bam \
                    --OUTPUT ${outDir}/gatk/${normal}_mkdp.bam \
                    --METRICS_FILE ${outDir}/gatk/${normal}_metrics.txt \
                    --REMOVE_DUPLICATES false \
                    --CREATE_INDEX true || exit 1
            samtools flagstat ${outDir}/gatk/${normal}_mkdp.bam > ${outDir}/gatk/${normal}_mkdp.bam.flagstat || exit 1

done < "${sampleNamesFile}"
echo "Duplicates marking done"


#---------- Step5: BASE QUALITY SCORE RECALIBRATION ----------
echo "Performing BQSR..."
while read -r sample tumor normal;
do
            echo "Performing BQSR for $tumor and $normal"
            gatk BaseRecalibrator \
                    -R ${ref} \
                    -O ${outDir}/recal/${tumor}_recal_data.table \
                    -I ${outDir}/gatk/${tumor}_mkdp.bam \
                    --known-sites ${refDir}/Mills_and_1000G_gold_standard.indels.hg38.vcf \
                    --known-sites ${refDir}/Homo_sapiens_assembly38.known_indels.vcf \
                    --known-sites ${refDir}/1000G_phase1.snps.high_confidence.hg38.vcf || exit 1

            gatk BaseRecalibrator \
                    -R ${ref} \
                    -I ${outDir}/gatk/${normal}_mkdp.bam \
                    -O ${outDir}/recal/${normal}_recal_data.table \
                    --known-sites ${refDir}/Mills_and_1000G_gold_standard.indels.hg38.vcf \
                    --known-sites ${refDir}/Homo_sapiens_assembly38.known_indels.vcf \
                    --known-sites ${refDir}/1000G_phase1.snps.high_confidence.hg38.vcf || exit 1

            gatk ApplyBQSR \
                    -R ${ref} \
                    -I ${outDir}/gatk/${tumor}_mkdp.bam \
                    -O ${outDir}/recal/${tumor}_recal.bam \
                    --bqsr-recal-file ${outDir}/recal/${tumor}_recal_data.table || exit 1

            gatk ApplyBQSR \
                    -R ${ref} \
                    -I ${outDir}/gatk/${normal}_mkdp.bam \
                    -O ${outDir}/recal/${normal}_recal.bam \
                    --bqsr-recal-file ${outDir}/recal/${normal}_recal_data.table || exit 1

done < "${sampleNamesFile}"
echo "BQSR done"


# #---------- Step6: SOMATIC VARIANT CALLING ----------
echo "Calling variants..."
while read -r sample tumor normal;
do
            call_variants() {
                echo "Calling variants for $tumor and $normal"
                gatk Mutect2 \
                    -R ${ref} \
                    -I ${outDir}/recal/${tumor}_recal.bam \
                    -I ${outDir}/recal/${normal}_recal.bam \
                    -O ${outDir}/mut2/${tumor}_${normal}.vcf \
                    --java-options "-Xmx16g -Djava.io.tmpdir=tmp" \
                    --tumor-sample ${tumor} \
                    --normal-sample ${normal} \
                    --germline-resource ${refDir}/af-only-gnomad.hg38.vcf \
                    --tmp-dir ${outDir}/tempDir
            }
            call_variants || exit 1

done < "${sampleNamesFile}"
echo "Variant calling done"


#---------- Step7: VARIANT FILTERING ----------
echo "Filtering variants..."
while read -r sample tumor normal;
do
            filter_variants() {
                local mut2_vcf="${outDir}/mut2/${tumor}_${normal}.vcf"

                echo "Getting pileup summaries for ${sample}_tumor at sites of known mutations"
                gatk GetPileupSummaries \
                    -I ${outDir}/recal/${tumor}_recal.bam \
                    -O ${outDir}/pileups/${tumor}.pileups.table \
                    --variant ${refDir}/af-only-gnomad.hg38.vcf \
                    --intervals ${refDir}/Illumina_Exome_TargetedRegions_v1.2.hg38.bed \
                    --tmp-dir ${outDir}/tempDir \
                    --java-options "-Xmx16g -Djava.io.tmpdir=tmp"

                echo "Getting pileup summaries for ${sample}_normal at sites of known mutations"
                gatk GetPileupSummaries \
                    -I ${outDir}/recal/${normal}_recal.bam \
                    -O ${outDir}/pileups/${normal}.pileups.table \
                    --variant ${refDir}/af-only-gnomad.hg38.vcf \
                    --intervals ${refDir}/Illumina_Exome_TargetedRegions_v1.2.hg38.bed \
                    --tmp-dir ${outDir}/tempDir \
                    --java-options "-Xmx16g -Djava.io.tmpdir=tmp"

                echo "Estimating contamination for $sample"
                gatk CalculateContamination \
                    -I ${outDir}/pileups/${tumor}.pileups.table \
                    -O ${outDir}/pileups/${tumor}_${normal}.contamination.table \
                    --matched-normal ${outDir}/pileups/${normal}.pileups.table \
                    --tmp-dir ${outDir}/tempDir \
                    --java-options "-Xmx16g -Djava.io.tmpdir=tmp"


                echo "Applying filter to variant calls for $sample"
                gatk FilterMutectCalls \
                    -O ${outDir}/mut2/${tumor}_${normal}.Filtered.vcf \
                    --variant ${mut2_vcf} \
                    --contamination-table ${outDir}/pileups/${tumor}_${normal}.contamination.table \
                    --reference $ref \
                    --tmp-dir ${outDir}/tempDir \
                    --java-options "-Xmx16g -Djava.io.tmpdir=tmp"

                awk \
                    '/^#/ {print $0; next} $7=="PASS" {print $0}' \
                    ${outDir}/mut2/${tumor}_${normal}.Filtered.vcf > \
                    ${outDir}/mut2/${tumor}_${normal}.PASS.vcf
            }
            filter_variants || exit 1

done < "${sampleNamesFile}"
echo "Variant filtering done"


#---------- Step8: VARIANT ANNOTATION & SUMMARY STATS ----------
echo "Annotating variants and generating summary statistics..."
while read -r sample tumor normal;
do
            annotate_variants() {
                echo "Annotating variants for $sample"
                java \
                    -jar ${snpEff} hg38 \
                    ${outDir}/mut2/${tumor}_${normal}.PASS.vcf > \
                    ${outDir}/annotations/${tumor}_${normal}_annotated.vcf
                wait
                bcftools stats \
                    ${outDir}/annotations/${tumor}_${normal}_annotated.vcf \
                    > ${outDir}/annotations/${tumor}_${normal}.variantSummary
            }
            annotate_variants

done < "${sampleNamesFile}"
echo "Variant annotation done"

## Report the total duration of the process
duration=${SECONDS}
echo "
The whole process took:
$(($duration / 3600)) hours: $((($duration % 3600) / 60)) minutes: $(($duration % 60)) seconds
"
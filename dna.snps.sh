#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --qos=medium
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16gb
#SBATCH --output=/users/sean.montgomery/logs/gatk.genome.txt
# === end SBATCH directives ===

# === begin ENVIRONMENT SETUP ===

#1. set directory containing sample folders
work_folder=/scratch-cbe/users/sean.montgomery/lab/gatk
#2. Fasta file of genome with only autosomes
index1=/groups/berger/lab/cluster_files/sean/Tak1v5.autosome/Tak1v5.autosome.fasta
#3. File names of fastq files to create - Change depending on if working with Cam2 or Tak2 SNPs relative to Tak1
fastq1=/scratch-cbe/users/sean.montgomery/lab/gatk/Cam2-final_1.fastq
fastq2=/scratch-cbe/users/sean.montgomery/lab/gatk/Cam2-final_2.fastq
# fastq1=/scratch-cbe/users/sean.montgomery/lab/gatk/Tak2_NIBB_1.fastq.gz
# fastq2=/scratch-cbe/users/sean.montgomery/lab/gatk/Tak2_NIBB_2.fastq.gz
#4. Output basename - Change depending on if working with Cam2 or Tak2 SNPs relative to Tak1
output=Cam2.autosome
# output=Takv6.Tak2snp

#Load modules
module load samtools/1.9-foss-2018b
module load gatk/4.0.1.2-java-1.8
module load picard/2.18.27-java-1.8
module load star/2.7.1a-foss-2018b
module load bwa/0.7.17-foss-2018b
module load bowtie2/2.3.4.2-foss-2018b
module load bedtools/2.27.1-foss-2018b

mkdir -p $work_folder
cd $work_folder
export TMPDIR=$work_folder/tmp
# # make sure the directory exists
mkdir -p $TMPDIR

##Copy in and merge bam files - Only if Cam2 SNPs
cp /groups/berger/lab/Raw/demultiplexed/129848* /groups/berger/lab/Raw/demultiplexed/129849* /groups/berger/lab/Raw/demultiplexed/129850* /groups/berger/lab/Raw/demultiplexed/134442* /groups/berger/lab/Raw/demultiplexed/134443* /groups/berger/lab/Raw/demultiplexed/134444* $work_folder
samtools merge ${output}.merged.bam 129848_TGGGTTTCTAGATCGC_000000000-JBK88_1_20201007B_20201007.bam 129848_TGGGTTTCTAGATCGC_000000000-JCBY2_1_20201106B_20201106.bam 129848_TGGGTTTCTAGATCGC_HW3FYDRXX_1_20201120B_20201120.bam 129848_TGGGTTTCTAGATCGC_HW3FYDRXX_2_20201120B_20201120.bam 129848_TGGGTTTCTAGATCGC_H72FGBGXH_1_20201210B_20201210.bam 129849_TTGACCCTTAGATCGC_000000000-JBK88_1_20201007B_20201007.bam 129849_TTGACCCTTAGATCGC_000000000-JCBY2_1_20201106B_20201106.bam 129849_TTGACCCTTAGATCGC_HW3FYDRXX_1_20201120B_20201120.bam 129849_TTGACCCTTAGATCGC_HW3FYDRXX_2_20201120B_20201120.bam 129849_TTGACCCTTAGATCGC_H72FGBGXH_1_20201210B_20201210.bam 129850_CCACTCCTTAGATCGC_000000000-JBK88_1_20201007B_20201007.bam 129850_CCACTCCTTAGATCGC_000000000-JCBY2_1_20201106B_20201106.bam 129850_CCACTCCTTAGATCGC_HW3FYDRXX_1_20201120B_20201120.bam 129850_CCACTCCTTAGATCGC_HW3FYDRXX_2_20201120B_20201120.bam 129850_CCACTCCTTAGATCGC_H72FGBGXH_1_20201210B_20201210.bam 134442_GTCGTGATTAGATCGC_000000000-JCBY2_1_20201106B_20201106.bam 134442_GTCGTGATTAGATCGC_HW3FYDRXX_1_20201120B_20201120.bam 134442_GTCGTGATTAGATCGC_HW3FYDRXX_2_20201120B_20201120.bam 134442_GTCGTGATTAGATCGC_H72FGBGXH_1_20201210B_20201210.bam 134443_ACCACTGTTAGATCGC_000000000-JCBY2_1_20201106B_20201106.bam 134443_ACCACTGTTAGATCGC_HW3FYDRXX_1_20201120B_20201120.bam 134443_ACCACTGTTAGATCGC_HW3FYDRXX_2_20201120B_20201120.bam 134443_ACCACTGTTAGATCGC_H72FGBGXH_1_20201210B_20201210.bam 134444_TGGTCACATAGATCGC_000000000-JCBY2_1_20201106B_20201106.bam 134444_TGGTCACATAGATCGC_HW3FYDRXX_1_20201120B_20201120.bam 134444_TGGTCACATAGATCGC_HW3FYDRXX_2_20201120B_20201120.bam 134444_TGGTCACATAGATCGC_H72FGBGXH_1_20201210B_20201210.bam
samtools sort -n -o ${output}.qsort -O bam ${output}.merged.bam

##Create fastq files - Only if Cam2 SNPs
bedtools bamtofastq -i ${output}.qsort -fq ${fastq1} -fq2 ${fastq2}
rm ${output}.qsort 

## Make sequence dictionary
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
  R=$index1 \
  O=/groups/berger/lab/cluster_files/sean/Tak1v5/Tak1v5.dict

## Convert fastq to bam
java -Xmx8G -jar $EBROOTPICARD/picard.jar FastqToSam \
    FASTQ= $fastq1 \
    FASTQ2= $fastq2 \
    OUTPUT=$output.bam \
    READ_GROUP_NAME=readgroup \
    SAMPLE_NAME=$output \
    LIBRARY_NAME=$output-Lib \
    PLATFORM_UNIT=something \
    PLATFORM=illumina \
    SEQUENCING_CENTER=NGS \
    RUN_DATE=2014-08-20T00:00:00-0400 \
    TMP_DIR=$TMPDIR
samtools index $output.bam

## Mark adaptor sequences
java -Xmx8G -jar $EBROOTPICARD/picard.jar MarkIlluminaAdapters \
    I=$output.bam \
    O=${output}_markilluminaadapters.bam \
    M=${output}_markilluminaadapters_metrics.txt \
    TMP_DIR=$TMPDIR

## Map to reference with BWA, MergeBamAlignments
bwa index -a bwtsw $index1
set -o pipefail
java -Xmx8G -jar $EBROOTPICARD/picard.jar SamToFastq \
    I=${output}_markilluminaadapters.bam \
    FASTQ=/dev/stdout \
    CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
    TMP_DIR=$TMPDIR | \
    bwa mem -M -t 7 -p $index1 /dev/stdin | \
    java -Xmx16G -jar $EBROOTPICARD/picard.jar MergeBamAlignment \
    ALIGNED_BAM=/dev/stdin \
    UNMAPPED_BAM=${output}.bam \
    OUTPUT=${output}_mapped.bam \
    R=$index1 CREATE_INDEX=true ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
    TMP_DIR=$TMPDIR

## Mark duplicates with MarkDuplicates, SortSam
java -Xmx32G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    INPUT=${output}_mapped.bam \
    OUTPUT=${output}_mapped_markduplicates.bam \
    METRICS_FILE=${output}_mapped_markduplicates_metrics.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    CREATE_INDEX=true \
    TMP_DIR=$TMPDIR
java -jar $EBROOTPICARD/picard.jar SortSam \
     INPUT=${output}_mapped_markduplicates.bam \
     OUTPUT=${output}_mapped_sorted.bam \
     SORT_ORDER=coordinate \
     TMP_DIR=$TMPDIR
     CREATE_INDEX=true

## Use HaplotypeCaller from GATK to identify SNPs
samtools index ${output}_mapped_sorted.bam
gatk HaplotypeCaller \
    -R $index1 \
    -I ${output}_mapped_sorted.bam \
    --sample-ploidy 1 \
    --standard-min-confidence-threshold-for-calling 30 \
    -O ${output}.raw.snps.indels.vcf

## Hard filter variants
 gatk VariantFiltration \
   -R $index1 \
   -O ${output}.filtered.vcf \
   -V ${output}.raw.snps.indels.vcf \
   --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
   --filter-name "my_snp_filter" 

gatk SelectVariants \
    -R $index1 \
    -V ${output}.raw.snps.indels.vcf \
    --select-type-to-include INDEL \
    -O ${output}.raw_indels.vcf 

gatk VariantFiltration \
    -R $index1 \
    -V ${output}.raw_indels.vcf  \
    --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filter-name "my_indel_filter" \
    -O ${output}.filtered_indels.vcf 

## Filter high quality SNPs
gatk SelectVariants \
   -R $index1 \
   -V ${output}.filtered.vcf \
   -O ${output}.filtered_snps.vcf \
   --select-type-to-include SNP \
   --selectExpressions 'vc.isNotFiltered()'



#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --qos=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=1-10
#SBATCH --output=/users/sean.montgomery/logs/bedmap.Cam2.PE.%a.txt
# === end SBATCH directives ===

#Script to trim adapters then align CUT&RUN paired end data with Bowtie2
#Script adapted from Michael Borg

# === begin ENVIRONMENT SETUP ===

#1. list of sample names on separate lines (eg copied from excel)
sample_list=/groups/berger/user/sean.montgomery/Documents/Scripts/input_files/cutrun.Cam2.PE.txt #10 Cam2 PE
#2. Genome mapped to
output=Takv6.Cam2snp #size=220358409
chromsize=220358409
#3. set directory containing sample folders
work_folder=/scratch-cbe/users/sean.montgomery/cutrun/${output}
#4. Data type "SE" or "PE"
datatype="PE"
#5. Annotation file
annotation_file=/groups/berger/lab/cluster_files/sean/${output}/annotations/${output}.gene.PCG.bed

F=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`
#4. sample name
NAME=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $2}'`

# Load the required modules
module load bedtools/2.27.1-foss-2018b
module load deeptools/2.5.4-foss-2018b-python-2.7.15
module load samtools/1.9-foss-2018b
# ... and then change to  working directory:
mkdir -p $work_folder/${NAME}
cd $work_folder
# === end ENVIRONMENT SETUP ===

##Load in bam files
cp /groups/berger/user/sean.montgomery/Documents/cutrun/${output}/*/merged/expression/bam/${NAME}.sizednuc150.bam* ${NAME}

##Calculate coverage in bedgraph output
bamCoverage -b ${NAME}/${NAME}.sizednuc150.bam -o ${NAME}/${NAME}.sizednuc150.bedgraph -of bedgraph -bs 10 --normalizeTo1x ${chromsize}

##Map chromatin coverage to genes
bedtools map -o mean -c 4 -a ${annotation_file} -b ${NAME}/${NAME}.sizednuc150.bedgraph > ${NAME}/${NAME}.sizednuc150.PCG.bed


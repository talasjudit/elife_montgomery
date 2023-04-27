#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --qos=short
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --output=/users/sean.montgomery/logs/snpsplit.prep.txt
# === end SBATCH directives ===

# === begin ENVIRONMENT SETUP ===

#1. set directory containing sample folders
work_folder=/scratch-cbe/users/sean.montgomery/lab/gatk
#2. Fasta file of genome
index1=/groups/berger/lab/cluster_files/sean/Tak1v5.autosome/Tak1v5.autosome.fasta
#3. Output basename - Change depending on if working with Cam2 or Tak2 SNPs relative to Tak1
output=Takv6.Cam2snp
# output=Takv6.Tak2snp
#4. Annotation file of Protein Coding Genes (.gff)
annotation_file=/groups/berger/lab/cluster_files/sean/${output}/annotations/${output}.gene.PCG.gff

#Load modules
module load samtools/1.9-foss-2018b
module load gatk/3.8-1-java-1.8
module load picard/2.18.27-java-1.8
module load bowtie2/2.3.4.2-foss-2018b
module load star/2.7.1a-foss-2018b
module load rsem/1.3.2-foss-2018b

cd $work_folder
export TMPDIR=$work_folder/tmp
# # make sure the directory exists
mkdir -p $TMPDIR

## Create N-masked SNPs in new reference genome
mkdir -p /scratch-cbe/users/sean.montgomery/lab/${output}
cp /groups/berger/lab/cluster_files/sean/${output}/${output}.fasta /scratch-cbe/users/sean.montgomery/lab/${output}
java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
   -T FastaAlternateReferenceMaker \
   -R $index1 \
   -o /scratch-cbe/users/sean.montgomery/lab/${output}/${output}.fasta \
   -V ${output}.filtered_snps.vcf \
   --snpmask ${output}.filtered_snps.vcf \
   --snpmaskPriority 
sed 's/.*chr/>chr/g' /scratch-cbe/users/sean.montgomery/lab/${output}/${output}.fasta | sed 's/.*unplaced/>unplaced/g' | sed 's/:1//g' > /scratch-cbe/users/sean.montgomery/lab/${output}/tmp.fasta
mv /scratch-cbe/users/sean.montgomery/lab/${output}/tmp.fasta /scratch-cbe/users/sean.montgomery/lab/${output}/${output}.fasta

##Make samtools index
samtools faidx /scratch-cbe/users/sean.montgomery/lab/${output}/${output}.fasta

##Make bowtie2 index
bowtie2-build /scratch-cbe/users/sean.montgomery/lab/${output}/${output}.fasta /scratch-cbe/users/sean.montgomery/lab/${output}/${output}

##Make .dict
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
	R=/scratch-cbe/users/sean.montgomery/lab/${output}/${output}.fasta \
	O=/scratch-cbe/users/sean.montgomery/lab/${output}/${output}.dict

##Make rsem index
mkdir -p /scratch-cbe/users/sean.montgomery/lab/${output}/rsem-index
rsem-prepare-reference --num-threads 8 --gff3 $annotation_file --star \
  /scratch-cbe/users/sean.montgomery/lab/${output}/${output}.fasta /scratch-cbe/users/sean.montgomery/lab/${output}/rsem-index/rsem-index_${output}

##Make STAR index
mkdir /scratch-cbe/users/sean.montgomery/lab/${output}/star-index_${output}/
STAR --runMode genomeGenerate \
	--genomeDir /scratch-cbe/users/sean.montgomery/lab/${output}/star-index_${output} \
	--genomeFastaFiles /scratch-cbe/users/sean.montgomery/lab/${output}/rsem-index/rsem-index_${output}.idx.fa \

## Format vcf into tab-delimited SNPs_<strain_name>/chr<chromosome>.txt.gz file
## Must do solo because SNPsplit_genome_preparation executable only for mouse genomes
## SNP-ID	Chromosome	Position	Strand	Ref/SNP
grep -v '#' ${output}.filtered_snps.vcf | awk -v OFS='\t' '{print NR,$1,$2,1,$4"/"$5}' > ${output}.SNP.file.txt

##Copy out files
cp ${output}.SNP.file.txt ${output}_mapped_sorted.bam* ${output}.filtered_snps.vcf* /groups/berger/lab/cluster_files/sean/snps/
cp -r ${work_folder}/../${output}/ /groups/berger/lab/cluster_files/sean/

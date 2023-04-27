#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --qos=short
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=1-34
#SBATCH --output=/users/sean.montgomery/logs/mapping_Cam2_PE_%a.txt
# === end SBATCH directives ===

# === begin ENVIRONMENT SETUP ===

#1. list of sample names on separate lines (eg copied from excel)
sample_list=/groups/berger/user/sean.montgomery/Documents/Scripts/input_files/rnaseq.Cam2.PE.txt #34 Cam2 PE
# sample_list=/groups/berger/user/sean.montgomery/Documents/Scripts/input_files/rnaseq.Cam2.SE.txt #6 Cam2 SE
# sample_list=/groups/berger/user/sean.montgomery/Documents/Scripts/input_files/rnaseq.Tak2.SE.txt #3 Tak2 SE
#2. Single or paired end sequences
datatype="PE"
# datatype="SE"
#3. Output basename - Change depending on if working with Cam2 or Tak2 SNPs relative to Tak1
output=Takv6.Cam2snp #size=220358409
chromsize=220358409
# output=Takv6.Tak2snp #size=220358409
# chromsize=220358409
#4. set directory containing sample folders
work_folder=/scratch-cbe/users/sean.montgomery/rna-seq/${output}/PE
# work_folder=/scratch-cbe/users/sean.montgomery/rna-seq/${output}/SE
#6. Fasta file of genome
index1=/groups/berger/lab/cluster_files/sean/${output}/${output}
#7. STAR genome directory
star_dir=/groups/berger/lab/cluster_files/sean/${output}/star-index_${output}
#8. RSEM genome directory
rsem_dir=/groups/berger/lab/cluster_files/sean/${output}/rsem-index/rsem-index_${output}

##Read in .bam file names
F=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`
##Read in sample names
NAME=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $2}'`
FILE1=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $3}'`
FILE2=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $4}'`

#Load modules
module load samtools/1.9-foss-2018b
module load star/2.7.1a-foss-2018b
module load rsem/1.3.2-foss-2018b
module load bedtools/2.27.1-foss-2018b
module load picard/2.18.27-java-1.8
module load deeptools/2.5.4-foss-2018b-python-2.7.15
module load trim_galore/0.6.2-foss-2018b-python-3.6.6

mkdir -p $work_folder/${NAME}
cd $work_folder/${NAME}
export TMPDIR=$work_folder/${NAME}/tmp
# # make sure the directory exists

##Check if input is bam and sort and split into fastq
if [[ ${F} == *bam ]]; then
	cp /groups/berger/lab/Raw/demultiplexed/${F} $work_folder/${NAME}
	cp /groups/berger/lab/Sean/bam/${F} $work_folder/${NAME}
	cp /groups/berger/lab/Sean/For_Elin/${F} $work_folder/${NAME}
        cp /resources/ngs/berger/*/${F} $work_folder/${NAME}/
	samtools sort -n -o ${NAME}.qsort -O bam ${F}
	if [ $datatype == "PE" ]; then
		bedtools bamtofastq -i ${NAME}.qsort -fq ${NAME}.end1.fq -fq2 ${NAME}.end2.fq
	elif [ $datatype == "SE" ]; then
		bedtools bamtofastq -i ${NAME}.qsort -fq ${NAME}_1.fastq
	else
		echo "Datatype??"
	fi
	echo "BAM converted to SAM for sample ... "
	echo ${NAME}
	rm ${NAME}.qsort
else
	if [[ $datatype == "PE" ]]; then
		if [[ "$F" == *-* ]]; then
			##This turned out to be unnecessary, therefore does not work...
			mv ../${FILE1}/*.fastq.gz .
			mv ../${FILE2}/*.fastq.gz .
			zcat *.fastq.gz | paste - - - - | awk '$2~ /1:N/ {print $1,$2"\n"$3"\n"$4"\n"$5}' > my_jgi_read1.fastq
		elif [[ "$F" == *RR* ]]; then
			gunzip ${F}_1.fastq.gz
			gunzip ${F}_2.fastq.gz
			mv ${F}_1.fastq ${NAME}.end1.fq
			mv ${F}_2.fastq ${NAME}.end2.fq
		else
			echo "You're working on thin ice, bud!"
		fi
	elif [ $datatype == "SE" ]; then
		gunzip ${F}.fastq.gz
		mv ${F}.fastq ${NAME}_1.fastq
	fi
fi

##Make directory for trimmed data
mkdir trimmed

##Map reads for PE data
if [ $datatype == "PE" ]; then
        trim_galore --dont_gzip --stringency 4 -o trimmed/ --paired ${NAME}.end1.fq ${NAME}.end2.fq
        rm -r $TMPDIR
		STAR --genomeDir $star_dir \
			--readFilesIn trimmed/${NAME}.end1_val_1.fq trimmed/${NAME}.end2_val_2.fq \
			--outFileNamePrefix ${NAME} \
			--outTmpDir $TMPDIR \
			--runThreadN 8 \
			--outSAMtype BAM Unsorted \
			--alignEndsType EndToEnd --alignIntronMax 1 --alignIntronMin 2 --scoreDelOpen -10000 --scoreInsOpen -10000 \
			--outSAMattributes NH HI AS NM MD 
		samtools view -f2 -h ${NAME}Aligned.out.bam > ${NAME}.Aligned.out.paired.bam
		mkdir -p $TMPDIR
		rsem-calculate-expression \
			--alignments ${NAME}.Aligned.out.paired.bam \
			$rsem_dir \
			${NAME}.${output} \
			--paired-end \
			--num-threads 8 \
			--temporary-folder $TMPDIR \
			--append-names \
			--estimate-rspd \
			--output-genome-bam \
			--seed 12345 \
			--calc-ci \
			--ci-memory 40000 \
			--sort-bam-by-coordinate
		#run the rsem plot function
		rsem-plot-model ${NAME}.${output} ${NAME}.${output}.pdf
		# # generate alignment summary tables
		# ### for aligned.bam (i,e to get raw mapping information)
		echo ".............................................................. "
		echo "---> Summarising alignment statistics for sample ... " ${NAME}
		samtools stats ${NAME}Aligned.out.bam > ${NAME}.Aligned.out_stats.txt
		samtools stats ${NAME}.Aligned.out.paired.bam > ${NAME}.Aligned.out.paired_stats.txt
		for stat in "raw total sequences:" "reads mapped:" "reads mapped and paired:" "reads properly paired:" "reads duplicated:" "average length:" "maximum length:" "average quality:" "insert size average:" "insert size standard deviation:" "inward oriented pairs:" "outward oriented pairs:" "pairs with other orientation:" "pairs on different chromosomes:" ; do
		        grep "^SN" ${NAME}.Aligned.out_stats.txt| grep "$stat" | awk -F $'\t' '{print $2,"\t",$3}' >> ${NAME}.aligned_summary.txt
		        grep "^SN" ${NAME}.Aligned.out.paired_stats.txt | grep "$stat" | awk -F $'\t' '{print $2,"\t",$3}' >> ${NAME}.paired_summary.txt
		done
		rm ${NAME}.Aligned.out_stats.txt
		rm ${NAME}.Aligned.out.paired_stats.txt
		##Filter out multimappers
		samtools view -h ${NAME}.${output}.genome.sorted.bam | awk '(($5 == "100") && ($17 == "ZW:f:1")) || ($1 ~ /^@/)' | samtools view -bS > ${NAME}.${output}.genome.sorted.filtered.bam
		samtools index ${NAME}.${output}.genome.sorted.filtered.bam
		bamCoverage -b ${NAME}.${output}.genome.sorted.filtered.bam -o ${NAME}.${output}.genome.sorted.filtered.bw --normalizeTo1x $chromsize --binSize=10
fi

##Map reads for SE data
if [ $datatype == "SE" ]; then
        trim_galore --dont_gzip --stringency 4 -o trimmed/ ${NAME}_1.fastq 
        rm -r $TMPDIR
		STAR --genomeDir $star_dir \
			--readFilesIn trimmed/${NAME}_1_trimmed.fq \
			--outFileNamePrefix ${NAME} \
			--outTmpDir $TMPDIR \
			--runThreadN 8 \
			--outSAMtype BAM Unsorted \
			--alignEndsType EndToEnd --alignIntronMax 1 --alignIntronMin 2 --scoreDelOpen -10000 --scoreInsOpen -10000 \
			--outSAMattributes NH HI AS NM MD 
		mkdir -p $TMPDIR
		rsem-calculate-expression \
			--alignments ${NAME}Aligned.out.bam \
			$rsem_dir \
			${NAME}.${output} \
			--num-threads 8 \
			--temporary-folder $TMPDIR \
			--append-names \
			--estimate-rspd \
			--output-genome-bam \
			--seed 12345 \
			--calc-ci \
			--ci-memory 40000 \
			--sort-bam-by-coordinate
		#run the rsem plot function
		rsem-plot-model ${NAME}.${output} ${NAME}.${output}.pdf
		# # generate alignment summary tables
		# ### for aligned.bam (i,e to get raw mapping information)
		echo ".............................................................. "
		echo "---> Summarising alignment statistics for sample ... " ${NAME}
		samtools stats ${NAME}Aligned.out.bam > ${NAME}.Aligned.out_stats.txt
		for stat in "raw total sequences:" "reads mapped:" "reads mapped and paired:" "reads properly paired:" "reads duplicated:" "average length:" "maximum length:" "average quality:" "insert size average:" "insert size standard deviation:" "inward oriented pairs:" "outward oriented pairs:" "pairs with other orientation:" "pairs on different chromosomes:" ; do
		        grep "^SN" ${NAME}.Aligned.out_stats.txt| grep "$stat" | awk -F $'\t' '{print $2,"\t",$3}' >> ${NAME}.aligned_summary.txt
		done
		rm ${NAME}.Aligned.out_stats.txt
		##Filter out multimappers
		samtools view -h ${NAME}.${output}.genome.sorted.bam | awk '(($5 == "100") && ($17 == "ZW:f:1")) || ($1 ~ /^@/)' | samtools view -bS > ${NAME}.${output}.genome.sorted.filtered.bam
		samtools index ${NAME}.${output}.genome.sorted.filtered.bam
		bamCoverage -b ${NAME}.${output}.genome.sorted.filtered.bam -o ${NAME}.${output}.genome.sorted.filtered.bw --normalizeTo1x $chromsize --binSize=10
fi

##Remove trimmed .fq files
rm trimmed/*.fq

##Copy out files
cp ${NAME}.${output}.genome.sorted.filtered.bam* /groups/berger/user/sean.montgomery/Documents/rna-seq/${output}/bam/
cp ${NAME}.${output}.genome.sorted.filtered.bw /groups/berger/user/sean.montgomery/Documents/rna-seq/${output}/bw/
cp ${NAME}.${output}.genes.results /groups/berger/user/sean.montgomery/Documents/rna-seq/${output}/genes.results/
cp ${NAME}.${output}.pdf /groups/berger/user/sean.montgomery/Documents/rna-seq/${output}/mapping_stats/pdf/
cp ${NAME}.aligned_summary.txt /groups/berger/user/sean.montgomery/Documents/rna-seq/${output}/mapping_stats/txt/
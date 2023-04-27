#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --qos=short
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1gb
#SBATCH --array=1-10
#SBATCH --output=/users/sean.montgomery/logs/snpsplit.cutrun.%a.txt
# === end SBATCH directives ===

# === begin ENVIRONMENT SETUP ===

#1. list of sample names on separate lines (eg copied from excel)
sample_list=/groups/berger/user/sean.montgomery/Documents/Scripts/input_files/cutrun.Cam2.PE.txt #10 Cam2 PE
# sample_list=/groups/berger/user/sean.montgomery/Documents/Scripts/input_files/rnaseq.Cam2.PE.txt #34 Cam2 PE
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
work_folder=/scratch-cbe/users/sean.montgomery/cutrun/${output}
# work_folder=/scratch-cbe/users/sean.montgomery/rna-seq/${output}/PE
# work_folder=/scratch-cbe/users/sean.montgomery/rna-seq/${output}/SE
#5. SNP file
snps=/groups/berger/lab/cluster_files/sean/snps/${output}.SNP.file.txt

F=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`
#4. sample name
NAME=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $2}'`

#Load modules
module load samtools/1.9-foss-2018b

cd $work_folder
export TMPDIR=$work_folder/tmp
# # make sure the directory exists
mkdir -p $TMPDIR

# Run SNPsplit

##For CUT&RUN data
if [ $datatype == "PE" ]; then
	cd ${NAME}
	/groups/berger/user/sean.montgomery/Documents/SNPsplit-0.3.4/SNPsplit \
		--snp_file $snps \
		${NAME}.sizednuc150.bam \
		--no_sort \
		--paired 
fi

# ## For PE RNA-seq data
# if [ $datatype == "PE" ]; then
# 	runDir2=$work_folder/$NAME
# 	cd $runDir2
# 	/groups/berger/user/sean.montgomery/Documents/SNPsplit-0.3.4/SNPsplit \
# 		--snp_file $snps \
# 		${NAME}.${output}.genome.sorted.filtered.bam \
# 		--paired 
# fi

# ##For SE RNA-seq data
# if [ $datatype == "SE" ]; then
# 	runDir2=$work_folder/$NAME
# 	cd $runDir2
# 	/groups/berger/user/sean.montgomery/Documents/SNPsplit-0.3.4/SNPsplit \
# 		--snp_file $snps \
# 		${NAME}.${output}.genome.sorted.filtered.bam
# fi


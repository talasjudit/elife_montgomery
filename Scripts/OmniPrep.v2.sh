##Collection of all prep scripts in one place

##Add gene_ids to genes in clusters
# output=Cam2.autosome
# output=Tak1v5.chrU
output=Takv6.Cam2snp
# output=Takv6.Tak2snp

##Make file with all TPM values together from rsem output
cd /Volumes/groups/berger/user/sean.montgomery/Documents/rna-seq/${output}/rsem_everything/genes.results/
cut -f1 13daf-3.${output}.genes.results | sed 's/.*_M/M/g' > tmp.0.txt
for i in $(seq 2 42); do NAME=`awk NR==$i'{print $1}' ../samples.Takv6.Cam2snp.paper.txt`; awk -v OFS='\n' 'NR==1{print NAME} NR>1{print $6}' NAME="$NAME" ${NAME}.${output}.genes.results > tmp.${i}.txt; done
paste tmp.*.txt > ../../TPM.txt
rm tmp*
##Collection of all preparation scripts in one place
##################################################################################################################################################################
##Prepare TPM file
##################################################################################################################################################################

##Create single TPM file
tpmsingle <- read.table("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/TPM.txt", header = T, row.names = 1, check.names=FALSE)
tpmsample <- tpmsingle
write.csv(tpmsample, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/TPM_sample.csv", sep = ""))

##Calculate average expression levels
tpm <- as.data.frame(row.names(tpmsample))
colnames(tpm) <- c("gene_id")
for (sample in c("Tak1","Cam2","13daf","ez23-16daf","a-","an","ar")){
	tmpmean <- rowMeans(tpmsample[,c(grep(paste0("^",sample),colnames(tpmsample)))])
	tmpsd <- apply(tpmsample[,c(grep(paste0("^",sample),colnames(tpmsample)))],1, sd, na.rm = TRUE)
	tmpn <- apply(tpmsample[,c(grep(paste0("^",sample),colnames(tpmsample)))], 1, function(x) length(which(!is.na(x))))
	tpm <- cbind(tpm,tmpmean,tmpsd,tmpn)
	colnames(tpm)[(length(tpm)-2):length(tpm)] <- c(paste0(sample,"_tpmMean"), paste0(sample,"_tpmSD"), paste0(sample,"_tpmN"))
}
tpm <- tpm[,-1]
write.csv(tpm, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/TPM.csv", sep = ""))


##################################################################################################################################################################
##Prepare TPM file from Frank2014
##################################################################################################################################################################

##Create single TPM file
tpmsinglefrank <- read.table("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/TPM.txt", header = T, row.names = 1, check.names=FALSE)
tpmsamplefrank <- tpmsinglefrank
write.csv(tpmsamplefrank, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/TPM_sample.csv", sep = ""))

##Calculate average expression levels
tpmfrank <- as.data.frame(row.names(tpmsamplefrank))
colnames(tpmfrank) <- c("gene_id")
for (sample in c("tip","sp")){
    tmpmean <- rowMeans(tpmsamplefrank[,c(grep(paste0("^",sample),colnames(tpmsamplefrank)))])
    tmpsd <- apply(tpmsamplefrank[,c(grep(paste0("^",sample),colnames(tpmsamplefrank)))],1, sd, na.rm = TRUE)
    tmpn <- apply(tpmsamplefrank[,c(grep(paste0("^",sample),colnames(tpmsamplefrank)))], 1, function(x) length(which(!is.na(x))))
    tpmfrank <- cbind(tpmfrank,tmpmean,tmpsd,tmpn)
    colnames(tpmfrank)[(length(tpmfrank)-2):length(tpmfrank)] <- c(paste0(sample,"_tpmMean"), paste0(sample,"_tpmSD"), paste0(sample,"_tpmN"))
}
tpmfrank <- tpmfrank[,-1]
write.csv(tpmfrank, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/TPM.csv", sep = ""))



##################################################################################################################################################################
##Prepare allele counts summary files for RNA-seq
##################################################################################################################################################################

##List of file names in folder genecounts
filenames <- list.files("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/genecounts/")
##Get sample names from file names
samplenames <- gsub(".genes.counts.bed","",filenames)
##Initiate output table with gene ids
sampleGenecounts <- read.table(paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/genecounts/",filenames[1],sep=""))
##Initiate output table with gene ids
matratioSummary <- as.data.frame(sampleGenecounts$V4)
colnames(matratioSummary) <- c("gene_id")
##Initiate output table with gene ids
samplecountsSummary <- as.data.frame(sampleGenecounts$V4)
colnames(samplecountsSummary) <- c("gene_id")
##Initiate output table with gene ids
allelecountsSummary <- as.data.frame(samplecountsSummary$gene_id)
colnames(allelecountsSummary) <- c("gene_id")
##Initiate output table with gene ids
alleleratioSummary <- as.data.frame(samplecountsSummary$gene_id)
colnames(alleleratioSummary) <- c("gene_id")
##Initiate output table with gene ids
allelediffSummary <- as.data.frame(samplecountsSummary$gene_id)
colnames(allelediffSummary) <- c("gene_id")

##Loop through each genecounts file and create matratioSummary
for (i in 1:length(filenames)) {
    ##Read in sample genecounts table from 1 file
    sampleGenecounts <- read.table(paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/genecounts/",filenames[i],sep=""))
    ##Calculate mat ratio in new column with NA if sum(mat+pat)<5
    sampleGenecounts$matratio <- NA
    sampleGenecounts$matratio[(sampleGenecounts$V6+sampleGenecounts$V7)>=5] <- sampleGenecounts$V7[(sampleGenecounts$V6+sampleGenecounts$V7)>=5]/(sampleGenecounts$V6[(sampleGenecounts$V6+sampleGenecounts$V7)>=5]+sampleGenecounts$V7[(sampleGenecounts$V6+sampleGenecounts$V7)>=5])
    ##cbind ratio to matratioSummary
    matratioSummary <- cbind(matratioSummary, sampleGenecounts$matratio)
	##cbind ratio to samplecountsSummary
    samplecountsSummary <- cbind(samplecountsSummary, sampleGenecounts$V6, sampleGenecounts$V7)
    ##colname matratioSummary <- samplename
    colnames(matratioSummary)[i+1] <- paste0(samplenames[i],"_ratio")
    ##colname samplecountsSummary <- samplename
    colnames(samplecountsSummary)[2*i] <- paste0(samplenames[i],"_pat")
    colnames(samplecountsSummary)[2*i+1] <- paste0(samplenames[i],"_mat")
	line  <- as.table(matrix(c(gsub(".genes.counts.bed","",filenames[i]),sum(sampleGenecounts$V6[(sampleGenecounts$V6+sampleGenecounts$V7)>=5]), sum(sampleGenecounts$V6[(sampleGenecounts$V6+sampleGenecounts$V7)>=5])),ncol=3,byrow = T))
	write.table(line, file=paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/genecounts.summary.txt"),append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}

##Summarize counts across each day
for (day in c("Tak1","Cam2","13daf","ez23-16daf")) {
    for (allele in c("pat","mat")) {
        ##cbind ratio to allelecountsSummary
        allelecountsSummary[,paste0(day,"_",allele)] <- rowSums(samplecountsSummary[intersect(grep(paste0("^",day),colnames(samplecountsSummary)), grep(allele,colnames(samplecountsSummary)))])
    }
    alleleratioSummary[,paste0(day,"_ratio")] <- allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))]/(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])
	alleleratioSummary[,paste0(day,"_ratio")][(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])<5] <- NA
	allelediffSummary[,paste0(day,"_diff")] <- (allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))]-allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))])
	allelediffSummary[,paste0(day,"_diff")][(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])<5] <- NA
}

##Write csv table of maternal ratios for all genes and samples
write.csv(matratioSummary, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/matratio.csv", sep = ""))
##Write csv table of allelic counts for all genes and samples
write.csv(samplecountsSummary, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/snpcounts.csv", sep = ""))
##Write csv table of allelic counts for all genes and days
write.csv(allelecountsSummary, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/allelecountsSummary.csv", sep = ""))
##Write csv table of allelic ratios for all genes and days
write.csv(alleleratioSummary, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/alleleratioSummary.csv", sep = ""))
##Write csv table of allelic differences for all genes and days
write.csv(allelediffSummary, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/allelediffSummary.csv", sep = ""))

##Repeat for cutoff of 50 reads
##Initiate output table with gene ids
alleleratioSummary50 <- as.data.frame(sampleGenecounts$V4)
colnames(alleleratioSummary50) <- c("gene_id")
##Initiate output table with gene ids
allelediffSummary50 <- as.data.frame(sampleGenecounts$V4)
colnames(allelediffSummary50) <- c("gene_id")

##Summarize counts across each day
for (day in c("Tak1","Cam2","13daf","ez23-16daf","ez23-21daf","ez23-25daf")) {
    alleleratioSummary50[,paste0(day,"_ratio")] <- allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))]/(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])
    alleleratioSummary50[,paste0(day,"_ratio")][(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])<50] <- NA
    allelediffSummary50[,paste0(day,"_diff")] <- (allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))]-allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))])
    allelediffSummary50[,paste0(day,"_diff")][(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])<50] <- NA
}

##Write csv table of allelic ratios for all genes and days
write.csv(alleleratioSummary50, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/alleleratioSummary50.csv", sep = ""))
##Write csv table of allelic differences for all genes and days
write.csv(allelediffSummary50, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/allelediffSummary50.csv", sep = ""))



##################################################################################################################################################################
##Prepare allele counts summary files for CUT&RUN data
##################################################################################################################################################################

##List of file names in folder genecounts
filenames <- list.files("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/variants/genecounts/")
##Get sample names from file names
samplenames <- gsub("-H",".H",gsub("-Cam2xTak1-",".",gsub(".genes.counts.bed","",filenames)))
##Initiate output table with gene ids
sampleGenecounts <- read.table(paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/variants/genecounts/",filenames[1],sep=""))
##Initiate output table with gene ids
matratioCutrun <- as.data.frame(sampleGenecounts$V4)
colnames(matratioCutrun) <- c("gene_id")
##Initiate output table with gene ids
samplecountsCutrun <- as.data.frame(sampleGenecounts$V4)
colnames(samplecountsCutrun) <- c("gene_id")
##Initiate output table with gene ids
allelecountsCutrun <- as.data.frame(sampleGenecounts$V4)
colnames(allelecountsCutrun) <- c("gene_id")
##Initiate output table with gene ids
alleleratioCutrun <- as.data.frame(sampleGenecounts$V4)
colnames(alleleratioCutrun) <- c("gene_id")
##Initiate output table with gene ids
allelediffCutrun <- as.data.frame(samplecountsCutrun$gene_id)
colnames(allelediffCutrun) <- c("gene_id")

##Loop through each genecounts file
for (i in 1:length(filenames))
{
    ##Read in sample genecounts table from 1 file
    sampleGenecounts <- read.table(paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/variants/genecounts/",filenames[i],sep=""))
    ##Calculate mat ratio in new column with NA if sum(mat+pat)<5
    sampleGenecounts$matratio <- NA
    sampleGenecounts$matratio[(sampleGenecounts$V6+sampleGenecounts$V7)>=5] <- sampleGenecounts$V7[(sampleGenecounts$V6+sampleGenecounts$V7)>=5]/(sampleGenecounts$V6[(sampleGenecounts$V6+sampleGenecounts$V7)>=5]+sampleGenecounts$V7[(sampleGenecounts$V6+sampleGenecounts$V7)>=5])
    ##cbind ratio to matratioCutrun
    matratioCutrun <- cbind(matratioCutrun, sampleGenecounts$matratio)
    ##cbind ratio to samplecountsCutrun
    samplecountsCutrun <- cbind(samplecountsCutrun, sampleGenecounts$V6, sampleGenecounts$V7)
    ##colname matratioCutrun <- samplename
    colnames(matratioCutrun)[i+1] <- paste0(samplenames[i],"_ratio")
    ##colname matratioCutrun <- samplename
    colnames(samplecountsCutrun)[2*i] <- paste0(samplenames[i],"_pat")
    colnames(samplecountsCutrun)[2*i+1] <- paste0(samplenames[i],"_mat")
	line  <- as.table(matrix(c(gsub(".genes.counts.bed","",filenames[i]),sum(sampleGenecounts$V6[(sampleGenecounts$V6+sampleGenecounts$V7)>=5]), sum(sampleGenecounts$V7[(sampleGenecounts$V6+sampleGenecounts$V7)>=5])),ncol=3,byrow = T))
	write.table(line, file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/variants/genecounts.summary.txt"),append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}
##Summarize counts across each day
for (day in gsub("ratio","",colnames(matratioCutrun)[-1])) {
	allelediffCutrun[,paste0(day,"diff")] <- (samplecountsCutrun[intersect(grep(paste0("^",day),colnames(samplecountsCutrun)), grep("mat",colnames(samplecountsCutrun)))]-samplecountsCutrun[intersect(grep(paste0("^",day),colnames(samplecountsCutrun)), grep("pat",colnames(samplecountsCutrun)))])
	allelediffCutrun[,paste0(day,"diff")][(samplecountsCutrun[intersect(grep(paste0("^",day),colnames(samplecountsCutrun)), grep("pat",colnames(samplecountsCutrun)))]+samplecountsCutrun[intersect(grep(paste0("^",day),colnames(samplecountsCutrun)), grep("mat",colnames(samplecountsCutrun)))])<5] <- NA
}
##Write csv table of maternal ratios for all genes and samples
write.csv(matratioCutrun, file = paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/matratio.csv", sep = ""), row.names = F)
##Write csv table of allelic counts for all genes and samples
write.csv(samplecountsCutrun, file = paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/snpcounts.csv", sep = ""), row.names = F)
##Write csv table of allelic differences for all genes and days
write.csv(allelediffCutrun, file = paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/allelediffSummary.csv", sep = ""), row.names = F)

##Create table with mean matratios for C&R replicates
##Initiate output table with gene ids
matratioCutrunmean <- as.data.frame(sampleGenecounts$V4)
colnames(matratioCutrunmean) <- c("gene_id")
for (sample in c(unique(gsub(".[0-9]_ratio","",colnames(matratioCutrun))))[-1]){
    tryCatch({
        matratioCutrunmean <- cbind(matratioCutrunmean,rowMeans(matratioCutrun[,grep(paste0("^",sample,"-"),colnames(matratioCutrun))]))
        colnames(matratioCutrunmean)[length(matratioCutrunmean)] <- paste0(sample,"_ratio")
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
##Write csv table of maternal ratios for all genes and samples
write.csv(matratioCutrunmean, file = paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/matratiomean.csv", sep = ""), row.names = F)









##List of file names in folder genecounts
filenames <- list.files("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/variants/genecounts/")
##Get sample names from file names
samplenames <- gsub("-H",".H",gsub("-Cam2xTak1-",".",gsub(".genes.counts.bed","",filenames)))
##Initiate output table with gene ids
sampleGenecounts <- read.table(paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/variants/genecounts/",filenames[1],sep=""))
##Initiate output table with gene ids
matratioCutrun10 <- as.data.frame(sampleGenecounts$V4)
colnames(matratioCutrun10) <- c("gene_id")
##Initiate output table with gene ids
samplecountsCutrun <- as.data.frame(sampleGenecounts$V4)
colnames(samplecountsCutrun) <- c("gene_id")
##Initiate output table with gene ids
allelediffCutrun10 <- as.data.frame(samplecountsCutrun$gene_id)
colnames(allelediffCutrun10) <- c("gene_id")

##Loop through each genecounts file
for (i in 1:length(filenames))
{
    ##Read in sample genecounts table from 1 file
    sampleGenecounts <- read.table(paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/variants/genecounts/",filenames[i],sep=""))
    ##Calculate mat ratio in new column with NA if sum(mat+pat)<5
    sampleGenecounts$matratio <- NA
    sampleGenecounts$matratio[(sampleGenecounts$V6+sampleGenecounts$V7)>=10] <- sampleGenecounts$V7[(sampleGenecounts$V6+sampleGenecounts$V7)>=10]/(sampleGenecounts$V6[(sampleGenecounts$V6+sampleGenecounts$V7)>=10]+sampleGenecounts$V7[(sampleGenecounts$V6+sampleGenecounts$V7)>=10])
    ##cbind ratio to matratioCutrun10
    matratioCutrun10 <- cbind(matratioCutrun10, sampleGenecounts$matratio)
    ##cbind ratio to samplecountsCutrun
    samplecountsCutrun <- cbind(samplecountsCutrun, sampleGenecounts$V6, sampleGenecounts$V7)
    ##colname matratioCutrun10 <- samplename
    colnames(matratioCutrun10)[i+1] <- paste0(samplenames[i],"_ratio")
    ##colname matratioCutrun10 <- samplename
    colnames(samplecountsCutrun)[2*i] <- paste0(samplenames[i],"_pat")
    colnames(samplecountsCutrun)[2*i+1] <- paste0(samplenames[i],"_mat")
    line  <- as.table(matrix(c(gsub(".genes.counts.bed","",filenames[i]),sum(sampleGenecounts$V6[(sampleGenecounts$V6+sampleGenecounts$V7)>=1]), sum(sampleGenecounts$V7[(sampleGenecounts$V6+sampleGenecounts$V7)>=1])),ncol=3,byrow = T))
    write.table(line, file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/variants/genecounts.summary.txt"),append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}
##Summarize counts across each day
for (day in gsub("ratio","",colnames(matratioCutrun10)[-1])) {
    allelediffCutrun10[,paste0(day,"diff")] <- (samplecountsCutrun[intersect(grep(paste0("^",day),colnames(samplecountsCutrun)), grep("mat",colnames(samplecountsCutrun)))]-samplecountsCutrun[intersect(grep(paste0("^",day),colnames(samplecountsCutrun)), grep("pat",colnames(samplecountsCutrun)))])
    allelediffCutrun10[,paste0(day,"diff")][(samplecountsCutrun[intersect(grep(paste0("^",day),colnames(samplecountsCutrun)), grep("pat",colnames(samplecountsCutrun)))]+samplecountsCutrun[intersect(grep(paste0("^",day),colnames(samplecountsCutrun)), grep("mat",colnames(samplecountsCutrun)))])<10] <- NA
}
##Write csv table of maternal ratios for all genes and samples
write.csv(matratioCutrun10, file = paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/matratio10.csv", sep = ""), row.names = F)
##Write csv table of allelic differences for all genes and days
write.csv(allelediffCutrun10, file = paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/allelediffSummary10.csv", sep = ""), row.names = F)

##Create table with mean matratios for C&R replicates
##Initiate output table with gene ids
matratioCutrun10mean <- as.data.frame(sampleGenecounts$V4)
colnames(matratioCutrun10mean) <- c("gene_id")
for (sample in c(unique(gsub(".[0-9]_ratio","",colnames(matratioCutrun10))))[-1]){
    tryCatch({
        matratioCutrun10mean <- cbind(matratioCutrun10mean,rowMeans(matratioCutrun10[,grep(paste0("^",sample,"-"),colnames(matratioCutrun10))]))
        colnames(matratioCutrun10mean)[length(matratioCutrun10mean)] <- paste0(sample,"_ratio10")
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
##Write csv table of maternal ratios for all genes and samples
write.csv(matratioCutrun10mean, file = paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/matratio10mean.csv", sep = ""), row.names = F)



##################################################################################################################################################################
##Prepare allele counts summary files for RNA-seq Frank2014
##################################################################################################################################################################

##List of file names in folder genecounts
filenames <- list.files("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/genecounts/")
##Get sample names from file names
samplenames <- gsub(".genes.counts.bed","",filenames)
##Initiate output table with gene ids
sampleGenecounts <- read.table(paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/genecounts/",filenames[3],sep=""))
##Initiate output table with gene ids
matratioSummary <- as.data.frame(sampleGenecounts$V4)
colnames(matratioSummary) <- c("gene_id")
##Initiate output table with gene ids
samplecountsSummary <- as.data.frame(sampleGenecounts$V4)
colnames(samplecountsSummary) <- c("gene_id")
##Initiate output table with gene ids
allelecountsSummary <- as.data.frame(samplecountsSummary$gene_id)
colnames(allelecountsSummary) <- c("gene_id")
##Initiate output table with gene ids
alleleratioSummary <- as.data.frame(samplecountsSummary$gene_id)
colnames(alleleratioSummary) <- c("gene_id")
##Initiate output table with gene ids
allelediffSummary <- as.data.frame(samplecountsSummary$gene_id)
colnames(allelediffSummary) <- c("gene_id")

##Loop through each genecounts file and create matratioSummary
for (i in 3:length(filenames)) {
    ##Read in sample genecounts table from 1 file
    sampleGenecounts <- read.table(paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/genecounts/",filenames[i],sep=""))
    ##Calculate mat ratio in new column with NA if sum(mat+pat)<5
    sampleGenecounts$matratio <- NA
    sampleGenecounts$matratio[(sampleGenecounts$V6+sampleGenecounts$V7)>=5] <- sampleGenecounts$V7[(sampleGenecounts$V6+sampleGenecounts$V7)>=5]/(sampleGenecounts$V6[(sampleGenecounts$V6+sampleGenecounts$V7)>=5]+sampleGenecounts$V7[(sampleGenecounts$V6+sampleGenecounts$V7)>=5])
    ##cbind ratio to matratioSummary
    matratioSummary <- cbind(matratioSummary, sampleGenecounts$matratio)
    ##cbind ratio to samplecountsSummary
    samplecountsSummary <- cbind(samplecountsSummary, sampleGenecounts$V6, sampleGenecounts$V7)
    ##colname matratioSummary <- samplename
    colnames(matratioSummary)[i+1-2] <- paste0(samplenames[i],"_ratio")
    ##colname samplecountsSummary <- samplename
    colnames(samplecountsSummary)[2*(i-2)] <- paste0(samplenames[i],"_pat")
    colnames(samplecountsSummary)[2*(i-2)+1] <- paste0(samplenames[i],"_mat")
    line  <- as.table(matrix(c(gsub(".genes.counts.bed","",filenames[i]),sum(sampleGenecounts$V6[(sampleGenecounts$V6+sampleGenecounts$V7)>=5]), sum(sampleGenecounts$V6[(sampleGenecounts$V6+sampleGenecounts$V7)>=5])),ncol=3,byrow = T))
    write.table(line, file=paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/genecounts.summary.txt"),append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}

##Summarize counts across each day
for (day in c("tip","sp")) {
    for (allele in c("pat","mat")) {
        ##cbind ratio to allelecountsSummary
        allelecountsSummary[,paste0(day,"_",allele)] <- rowSums(samplecountsSummary[intersect(grep(paste0("^",day),colnames(samplecountsSummary)), grep(allele,colnames(samplecountsSummary)))])
    }
    alleleratioSummary[,paste0(day,"_ratio")] <- allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))]/(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])
    alleleratioSummary[,paste0(day,"_ratio")][(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])<5] <- NA
    allelediffSummary[,paste0(day,"_diff")] <- (allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))]-allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))])
    allelediffSummary[,paste0(day,"_diff")][(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])<5] <- NA
}

##Write csv table of maternal ratios for all genes and samples
write.csv(matratioSummary, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/matratio.csv", sep = ""), row.names = F)
##Write csv table of allelic counts for all genes and samples
write.csv(samplecountsSummary, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/snpcounts.csv", sep = ""), row.names = F)
##Write csv table of allelic counts for all genes and days
write.csv(allelecountsSummary, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/allelecountsSummary.csv", sep = ""), row.names = F)
##Write csv table of allelic ratios for all genes and days
write.csv(alleleratioSummary, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/alleleratioSummary.csv", sep = ""), row.names = F)
##Write csv table of allelic differences for all genes and days
write.csv(allelediffSummary, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/allelediffSummary.csv", sep = ""), row.names = F)

##Repeat for cutoff of 50 reads
##Initiate output table with gene ids
alleleratioSummary50 <- as.data.frame(sampleGenecounts$V4)
colnames(alleleratioSummary50) <- c("gene_id")
##Initiate output table with gene ids
allelediffSummary50 <- as.data.frame(sampleGenecounts$V4)
colnames(allelediffSummary50) <- c("gene_id")

##Summarize counts across each day
for (day in c("tip","sp")) {
    alleleratioSummary50[,paste0(day,"_ratio")] <- allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))]/(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])
    alleleratioSummary50[,paste0(day,"_ratio")][(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])<50] <- NA
    allelediffSummary50[,paste0(day,"_diff")] <- (allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))]-allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))])
    allelediffSummary50[,paste0(day,"_diff")][(allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("pat",colnames(allelecountsSummary)))]+allelecountsSummary[intersect(grep(paste0("^",day),colnames(allelecountsSummary)), grep("mat",colnames(allelecountsSummary)))])<50] <- NA
}

##Write csv table of allelic ratios for all genes and days
write.csv(alleleratioSummary50, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/alleleratioSummary50.csv", sep = ""), row.names = F)
##Write csv table of allelic differences for all genes and days
write.csv(allelediffSummary50, file = paste("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/allelediffSummary50.csv", sep = ""), row.names = F)



##################################################################################################################################################################
##Prepare chromatin enrichment table for Takv6.Cam2snp
##################################################################################################################################################################

##List of file names in folder genecounts
filenameschrom <- list.files("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/expression/bed/pcg/")
##Initiate output table for chromatin enrichment
chromenrichtmp <- read.table(paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/expression/bed/pcg/",filenameschrom[1],sep=""))
##Initiate output table with gene ids
chromenrich <- as.data.frame(chromenrichtmp$V4)
colnames(chromenrich) <- c("gene_id")
##Get sample names from file names
samplenames <- gsub("-H",".H",gsub("-Cam2xTak1-",".",gsub(".sizednuc150.PCG.bed","",filenameschrom)))

##Loop through each genecounts file
for (i in 1:length(filenameschrom))
{
    ##Read in sample PCG chromatin enrichment from 1 file
    chromenrichtmp <- read.table(paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/expression/bed/pcg/",filenameschrom[i],sep=""))
    ##cbind bed score to chromenrich
    chromenrich <- cbind(chromenrich, chromenrichtmp$V6)
    ##Rename columns
    colnames(chromenrich)[i+1] <- paste0(samplenames[i],"_enrich")
}
##Write csv table of allelic differences for all genes and days
write.csv(chromenrich, file = paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/chromenrich.csv", sep = ""), row.names = F)

##Create table with mean matratios for C&R replicates
##Initiate output table with gene ids
chromenrichmean <- as.data.frame(chromenrichtmp$V4)
colnames(chromenrichmean) <- c("gene_id")
for (sample in c(unique(gsub(".[0-9]_enrich","",colnames(chromenrich))))[-1]){
    tryCatch({
        chromenrichmean <- cbind(chromenrichmean,rowMeans(chromenrich[,grep(paste0("^",sample,"-"),colnames(chromenrich))]))
        colnames(chromenrichmean)[length(chromenrichmean)] <- paste0(sample,"_enrich")
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
##Write csv table of maternal ratios for all genes and samples
write.csv(chromenrichmean, file = paste("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/chromenrichmean.csv", sep = ""), row.names = F)


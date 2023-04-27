##Collection of all plotting scripts in one place

##Load in packages
library("ggalluvial")
library("reshape2")
library(plyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library("ggVennDiagram")
data_summary <- function(data, varname, groupnames){
	require(plyr)
	summary_func <- function(x, col){
	c(mean = mean(x[[col]], na.rm=TRUE),
	sd = sd(x[[col]], na.rm=TRUE))
	}
	data_sum<-ddply(data, groupnames, .fun=summary_func,
	varname)
	# data_sum <- rename(data_sum, c("mean" = varname))
	return(data_sum)
}

##Load in TPM data
tpm <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/TPM.csv", header = T, row.names = 1, check.names=FALSE)
tpmsample <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/TPM_sample.csv", header = T, row.names = 1, check.names=FALSE)

##Load in chromatin enrichment data
chromenrich <- read.csv("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/chromenrich.csv", header = T, row.names = 1, check.names=FALSE)
chromenrichmean <- read.csv("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/chromenrichmean.csv", header = T, row.names = 1, check.names=FALSE)

##Load in summary tables
matratioSummary <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/matratio.csv", header = T, row.names = 2, check.names=FALSE)
samplecountsSummary <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/snpcounts.csv", header = T, row.names = 1, check.names=FALSE)
allelecountsSummary <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/allelecountsSummary.csv", header = T, row.names = 1, check.names=FALSE)
alleleratioSummary <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/alleleratioSummary.csv", header = T, row.names = 1, check.names=FALSE)
allelediffSummary <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/allelediffSummary.csv", header = T, row.names = 1, check.names=FALSE)

matratioCutrunmean <- read.csv("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/matratiomean.csv", header = T, row.names = 1, check.names=FALSE)
cutcountsSummary <- read.csv("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/snpcounts.csv", header = T, row.names = 1, check.names=FALSE)
cutdiffSummary <- read.csv("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/allelediffSummary.csv", header = T, row.names = 1, check.names=FALSE)
matratioCutrun <- read.csv("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/matratio.csv", header = T, row.names = 1, check.names=FALSE)
samplecountsCutrun <- read.csv("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/snpcounts.csv", header = T, row.names = 1, check.names=FALSE)
allelediffCutrun <- read.csv("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/allelediffSummary.csv", header = T, row.names = 1, check.names=FALSE)

##Read in table of DEGS 
ez2316dafvs13daf <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/rsem_ez/output/ez23-16daf_vs_13daf.csv",header = T,row.names = 2)
updeg <- row.names(subset(ez2316dafvs13daf,ez2316dafvs13daf$padj<0.01 & (ez2316dafvs13daf$log2FoldChange>=1)))
downdeg <- row.names(subset(ez2316dafvs13daf,ez2316dafvs13daf$padj<0.01 & ez2316dafvs13daf$log2FoldChange <= -1))

##Read in more stringent cutoffs for matratios
alleleratioSummary50 <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/alleleratioSummary50.csv", header = T, row.names = 1, check.names=FALSE)
matratioCutrun10 <- read.csv("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/matratio10.csv", header = T, row.names = 1, check.names=FALSE)
matratioCutrun10mean <- read.csv("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/matratio10mean.csv", header = T, row.names = 1, check.names=FALSE)
colnames(matratioCutrun10) <- paste0(colnames(matratioCutrun10),10)
colnames(alleleratioSummary50) <- paste0(colnames(alleleratioSummary50),50)

##Create single table with all summaries per genes: TPM, matratio TPM, matratio C&R, ass'n w/ peaks, DEG in ez2/3
omnisummary <- as.data.frame(cbind(tpm, chromenrichmean, alleleratioSummary,alleleratioSummary50, matratioCutrunmean,matratioCutrun10mean))
# omnisummary <- as.data.frame(cbind(tpm, chromenrich, alleleratioSummary50, matratioCutrun10))
for (i in 1:length(colnames(genelist))) {
    omnisummary[,paste0(colnames(genelist)[i],"_peak")] <- "No"
    omnisummary[,paste0(colnames(genelist)[i],"_peak")][row.names(omnisummary) %in% genelist[,i]] <- "Yes"
}
for (samp in c("13daf")){
    tmpsamp <- paste0(samp,"_ratio50")
    tmpsampmut <- paste0("ez23-16daf_ratio50")
    tmptpm <- paste0(samp,"_tpmMean")
    tmptpmmut <- paste0("ez23-16daf_tpmMean")
        omnisummary[,paste0(samp,"_tpmMeanDelta")] <- (omnisummary[,tmptpmmut] - omnisummary[,tmptpm])
        omnisummary[,paste0(samp,"_ratio50Delta")] <- (omnisummary[,tmpsampmut] - omnisummary[,tmpsamp])
    for (mark in c("H3K9me1","H3K36me3","H3K27me3","H3")){
        tryCatch({
            tmpenrich <- paste0(samp,".",mark,"_enrich")
            tmpenrichmut <- paste0("ez23-16daf.",mark,"_enrich")
            tmpmark <- paste0(samp,".",mark,"_ratio10")
            tmpmarkmut <- paste0("ez23-16daf.",mark,"_ratio10")
            omnisummary[,paste0(samp,".",mark,"_enrichDelta")] <- (omnisummary[,tmpenrichmut] - omnisummary[,tmpenrich])
            omnisummary[,paste0(samp,".",mark,"_ratio10Delta")] <- (omnisummary[,tmpmarkmut] - omnisummary[,tmpmark])
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
}

omnisummary$degor <- c("None")
omnisummary$degor[row.names(omnisummary) %in% downdeg] <- c("Down")
omnisummary$degor[row.names(omnisummary) %in% updeg] <- c("Up")

##Add Frank2014 RNA-seq data
tpmfrank <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/TPM.csv", header = T, row.names = 1, check.names=FALSE)
tpmsamplefrank <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/TPM_sample.csv", header = T, row.names = 1, check.names=FALSE)
matratioSummaryfrank <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/matratio.csv", header = T, row.names = 1, check.names=FALSE)
samplecountsSummaryfrank <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/snpcounts.csv", header = T, row.names = 1, check.names=FALSE)
allelecountsSummaryfrank <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/allelecountsSummary.csv", header = T, row.names = 1, check.names=FALSE)
alleleratioSummaryfrank <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/alleleratioSummary.csv", header = T, row.names = 1, check.names=FALSE)
allelediffSummaryfrank <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/allelediffSummary.csv", header = T, row.names = 1, check.names=FALSE)
alleleratioSummary50frank <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/alleleratioSummary50.csv", header = T, row.names = 1, check.names=FALSE)
colnames(alleleratioSummary50frank) <- paste0(colnames(alleleratioSummary50frank),50)

omnisummary$FrankTip_tpmMean <- tpmfrank$tip_tpmMean
omnisummary$FrankSp_tpmMean <- tpmfrank$sp_tpmMean
omnisummary$FrankTip_ratio50 <- alleleratioSummary50frank$tip_ratio50
omnisummary$FrankSp_ratio50 <- alleleratioSummary50frank$sp_ratio50

##Filter out genes with ambiguous SNP calls
omnisummary[unique(c(row.names(subset(omnisummary, omnisummary$Cam2_ratio<0.95)), row.names(subset(omnisummary,omnisummary$Tak1_ratio>0.05)))),grep("ratio",colnames(omnisummary))] <- NA

##Add gene info
snplist <- read.table("/groups/berger/lab/cluster_files/sean/snps/Takv6.Cam2snp.snpcountpergene.bed", row.names = 4)
snplist$snp <- c("No")
snplist[row.names(subset(snplist, V6 > 0)),]$snp <- c("Yes")
omnisummary$snp <- snplist$snp
omnisummary$chr <- snplist$V1
omnisummary$start <- snplist$V2
omnisummary$end <- snplist$V3
omnisummary$strand <- snplist$V5

##Plot rna matratio per sample by size
measurements <- read.table("/groups/berger/user/sean.montgomery/Documents/rna-seq/Dissection_pics/Sporophyte_measurements_repnum.txt",header = T)
measurements$Age <- factor(measurements$Age)
measurements$Sample <- factor(measurements$Sample)
measurements <- measurements[order(measurements$Sample), ]
meltdf <- matratioSummary[,grep("ratio",colnames(matratioSummary))]
meltdf <- melt(meltdf)
meltdfsummary <- data_summary(meltdf, varname="value", groupnames=c("variable"))
meltdfsummary$age <- gsub(".I.*","",gsub(".H.*","",gsub("_.*","",meltdfsummary$variable)))
meltdfsummary$exp <- gsub(".*I","I",gsub(".*H","H",gsub("_.*","",meltdfsummary$variable)))
meltdfsummary$exp[grep("[HI].*",meltdfsummary$exp,invert = T)] <- c("RNASeq")
meltdfsummary$age <- factor(meltdfsummary$age)
meltdfsummary <- meltdfsummary[1:153,]
meltdfsummary$Area <- measurements$Area[which(measurements$Sample %in% meltdfsummary$age)]
meltdfsummary$Stage <- as.numeric(measurements$Stage[which(measurements$Sample %in% meltdfsummary$age)])
meltdfsummary$alter <- measurements$Age[which(measurements$Sample %in% meltdfsummary$age)]
meltdfsummary$oida <- as.numeric(gsub("daf","",meltdfsummary$alter))
agecol <- c("12daf"="#CB181D","13daf"="#FB6A4A","14daf"="#4292C6","15daf"="#cdcd00","16daf"="#74C476")
ggplot(meltdfsummary[1:64,], aes(x=Area, y=mean)) +
    geom_point(aes(color=alter),size=2,position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=alter),position=position_dodge(.9)) +
    stat_cor(method = "spearman",size = 12) +
    scale_colour_manual(name="Legend",values=agecol) +
    theme_classic() +
    geom_hline(yintercept = 0.5) +
    theme(axis.text.x = element_text(color = "black",size = 36, angle = 45,hjust = 1),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)) +
    scale_x_continuous(limits = c(0,0.6)) +
    scale_y_continuous(name="Maternal ratio",limits = c(0,1.2),breaks = c(0,0.25,0.5,0.75,1))

pdf("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/plots/matratio.perday.perrep.embryo.pdf", width = 12, height = 12)
ggplot(meltdfsummary[11:21,], aes(x=age, y=mean)) +
    geom_point(aes(group=c(1:length(meltdfsummary[11:21,1]))),size=2,position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,group=c(1:length(meltdfsummary[11:21,1]))),width = 0.1,position=position_dodge(width = 0.9)) +
    stat_cor(method = "spearman",size = 12) +
    theme_classic() +
    geom_hline(yintercept = 0.5) +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)) +
    scale_x_discrete(name="Replicate number",labels=c("11","12","13","14","15","16","17","18","19","20","3")) +
    scale_y_continuous(name="Maternal ratio",limits = c(0,1.1),breaks = c(0,0.25,0.5,0.75,1))
dev.off()


ggplot(meltdfsummary[1:64,], aes(x=Stage, y=mean)) +
    geom_point(aes(group=c(1:length(meltdfsummary[1:64,1])),color=alter),size=2,position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,group=c(1:length(meltdfsummary[1:64,1])),color=alter),width = 0.1,position=position_dodge(width = 0.9)) +
    stat_cor(method = "spearman",size = 12) +
    theme_classic() +
    geom_hline(yintercept = 0.5) +
    theme(axis.text.x = element_text(color = "black",size = 36, angle = 45,hjust = 1),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)) +
    scale_x_continuous(name="Development stage",limits = c(1,5),breaks = c(2,3,4)) +
    scale_y_continuous(name="Maternal ratio",limits = c(0,1.2),breaks = c(0,0.25,0.5,0.75,1))

##Plot heatmap of RNA sample distances for paper
load("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/rsem_everything/output/dds.RData")
rld <- rlog(dds, blind=FALSE)
colors=colorRampPalette(c("white", "steelblue4"))(100)
sampleDists <- dist(t(assay(rld[,c(14:24,169:174,177:182)])))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld[,c(14:24,169:174,177:182)]$condition, sep="-")
colnames(sampleDistMatrix) <- paste(samples$sample[c(14:24,169:174,177:182)], sep="-")
pdf("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/rsem_everything/output/dist.1nvs2n.pdf", height=10, width=10)
pheatmap(sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors)
dev.off()

##Plot matratio for all samples and all marks
for (day in c("Tak1","Cam2","13daf","ez23-16daf")) {
	tmpday <- sym(paste0(day,"_ratio50"))
	tmpdaydiff <- sym(paste0(day,"_diff"))
	MyBreaks <- seq(-0.01,1.1,0.01)
	pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/graphs/",day,".genes.counts.prop.pdf"),height=12,width=12)
    print(ggplot(omnisummary, aes(!!tmpday, fill=cut(..x.., breaks=get("MyBreaks", envir=.GlobalEnv)))) +
        geom_vline(xintercept=0.5) +
	    geom_histogram(show.legend = FALSE, binwidth = 0.01) +
	    theme_classic() +
	    labs(x = "Transcription maternal ratio", y = "Number of resolved genes") +
        xlim(-0.01,1.01) +
        scale_fill_discrete(h = c(180, 360), c = 150, l = 80) +
        theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)))
	dev.off()	
    ##Plot for uneven breaks
    MyBreaks <- c(0,0.05,0.35,0.65,0.95,1)
    pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/graphs/",day,".genes.counts.prop.uneven.pdf"),height=12,width=12)
    print(ggplot(omnisummary, aes(!!tmpday, fill=cut(..x.., breaks=get("MyBreaks", envir=.GlobalEnv)))) +
        geom_vline(xintercept=0.5) +
        geom_histogram(show.legend = FALSE, breaks=MyBreaks) +
        theme_classic() +
        labs(x = "Transcription maternal ratio", y = "Number of resolved genes") +
        xlim(-0.01,1.01) + 
        scale_fill_discrete(h = c(180, 360), c = 150, l = 80) +
        theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)))
    dev.off()
    tryCatch({
        MyBreaks <- c(-0.0000001,0.0500001,0.3500001,0.6499999,0.9499999,1)
        tmp <- omnisummary %>% mutate(category=cut(omnisummary[,paste0(day,"_ratio50")], breaks=MyBreaks))
        pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/graphs/",day,".genes.counts.prop.uneven.percentage.pdf"),height=12,width=12)
        print(ggplot(tmp[!is.na(tmp$category),], aes(x = category)) +  
            geom_bar(aes(y = (..count..)/sum(..count..)), fill=c("#69c7bd","#5ecbe8","#bca9d2","#c572ae","#f175ad")) + 
            geom_text(stat = "count", aes(label = after_stat(count)), position="fill", size = 12) +
            scale_y_continuous(labels=percent) +
            scale_x_discrete(labels = c("Paternal","Paternal bias","No bias","Maternal bias","Maternal")) +
            theme_classic() +
            labs(x = "Transcription maternal ratio", y = "Percentage of resolved genes") +
            theme(axis.text.x = element_text(color = "black",size = 36,hjust = 1,angle = 45),
                axis.text.y = element_text(color = "black",size = 36),
                axis.title.x = element_blank(),
                axis.title.y = element_text(color = "black",size = 36),
                legend.title = element_text(color = "black",size = 36),
                legend.text = element_text(colour = "black",size = 36)))
        dev.off()
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

##Plot matratio for Frank2014
for (day in c("FrankTip","FrankSp")) {
    tmpday <- sym(paste0(day,"_ratio50"))
    tmpdaydiff <- sym(paste0(day,"_diff"))
    MyBreaks <- seq(-0.01,1.1,0.01)
    pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/graphs/",day,".genes.counts.prop.pdf"),height=12,width=12)
    print(ggplot(omnisummary, aes(!!tmpday, fill=cut(..x.., breaks=get("MyBreaks", envir=.GlobalEnv)))) +
        geom_vline(xintercept=0.5) +
        geom_histogram(show.legend = FALSE, binwidth = 0.01) +
        theme_classic() +
        labs(x = "Transcription maternal ratio", y = "Number of resolved genes") +
        xlim(-0.01,1.01) + 
        scale_fill_discrete(h = c(180, 360), c = 150, l = 80) +
        theme(axis.text.x = element_text(color = "black",size = 36),
                axis.text.y = element_text(color = "black",size = 36),
                axis.title.x = element_text(color = "black",size = 36),
                axis.title.y = element_text(color = "black",size = 36),
                legend.title = element_text(color = "black",size = 36),
                legend.text = element_text(colour = "black",size = 36)))
    dev.off()   
    ##Plot for uneven breaks
    MyBreaks <- c(0,0.05,0.35,0.65,0.95,1)
    pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/graphs/",day,".genes.counts.prop.uneven.pdf"),height=12,width=12)
    print(ggplot(omnisummary, aes(!!tmpday, fill=cut(..x.., breaks=get("MyBreaks", envir=.GlobalEnv)))) +
        geom_vline(xintercept=0.5) +
        geom_histogram(show.legend = FALSE, breaks=MyBreaks) +
        theme_classic() +
        xlim(-0.01,1.01) + 
        scale_fill_discrete(h = c(180, 360), c = 150, l = 80) +
        labs(x = "Transcription maternal ratio", y = "Number of resolved genes") +
        theme(axis.text.x = element_text(color = "black",size = 36),
                axis.text.y = element_text(color = "black",size = 36),
                axis.title.x = element_text(color = "black",size = 36),
                axis.title.y = element_text(color = "black",size = 36),
                legend.title = element_text(color = "black",size = 36),
                legend.text = element_text(colour = "black",size = 36)))
    dev.off()
    tryCatch({
    MyBreaks <- c(-0.0000001,0.0500001,0.3500001,0.6499999,0.9499999,1)
    tmp <- omnisummary %>% mutate(category=cut(omnisummary[,paste0(day,"_ratio50")], breaks=MyBreaks))
    pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Tak2snp/graphs/",day,".genes.counts.prop.uneven.percentage.pdf"),height=12,width=12)
    print(ggplot(tmp[!is.na(tmp$category),], aes(x = category)) +  
            geom_bar(aes(y = (..count..)/sum(..count..)), fill=c("#69c7bd","#5ecbe8","#bca9d2","#c572ae","#f175ad")) + 
            geom_text(stat = "count", aes(label = after_stat(count)), position="fill", size = 12) +
            scale_y_continuous(labels=percent) +
            scale_x_discrete(labels = c("Paternal","Paternal bias","No bias","Maternal bias","Maternal")) +
            theme_classic() +
            labs(x = "Transcription maternal ratio", y = "Percentage of resolved genes") +
            theme(axis.text.x = element_text(color = "black",size = 36,hjust = 1,angle = 45),
                axis.text.y = element_text(color = "black",size = 36),
                axis.title.x = element_blank(),
                axis.title.y = element_text(color = "black",size = 36),
                legend.title = element_text(color = "black",size = 36),
                legend.text = element_text(colour = "black",size = 36)))
    dev.off()
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

##Plot maternal ratios and scatterplots per cluster per mark per day
for (samp in c("13daf","Tak1","Cam2","ez23-16daf")){
	tmpsamp <- sym(paste0(samp,"_ratio50"))
	tmptpm <- sym(paste0(samp,"_tpmMean"))
	for (mark in c("H3K9me1","H3K36me3","H3K27me3","H3")){
		tryCatch({
			tmpenrich <- sym(paste0(samp,".",mark,"_enrich"))
			tmpmark <- sym(paste0(samp,".",mark,"_ratio10"))
            MyBreaks <- seq(-0.01,1.1,0.01)
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/hist/pdf/",samp,".",mark,".pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(!!tmpmark, fill=cut(..x.., breaks=get("MyBreaks", envir=.GlobalEnv)))) +
                    geom_vline(xintercept=0.5) +
                    geom_histogram(show.legend = FALSE, binwidth = 0.01) +
                    theme_classic() +
                    labs(x = paste0(mark," maternal ratio"), y = "Number of resolved genes") +
                    xlim(-0.01,1.01) +
                    scale_fill_discrete(h = c(180, 360), c = 150, l = 80) +
                    theme(axis.text.x = element_text(color = "black",size = 36),
                      axis.text.y = element_text(color = "black",size = 36),
                      axis.title.x = element_text(color = "black",size = 36),
                      axis.title.y = element_text(color = "black",size = 36),
                      legend.title = element_text(color = "black",size = 36),
                      legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            ##Plot for uneven breaks
                MyBreaks <- c(0,0.05,0.35,0.65,0.95,1)
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/hist/pdf/",samp,".",mark,".uneven.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(!!tmpmark, fill=cut(..x.., breaks=get("MyBreaks", envir=.GlobalEnv)))) +
                geom_vline(xintercept=0.5) +
                geom_histogram(show.legend = FALSE, breaks=MyBreaks) +
                theme_classic() +
                labs(x = paste0(mark," maternal ratio"), y = "Number of resolved genes") +
                xlim(-0.01,1.01) + 
                scale_fill_discrete(h = c(180, 360), c = 150, l = 80) +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            tryCatch({
            MyBreaks <- c(-0.0000001,0.0500001,0.3500001,0.6499999,0.9499999,1)
            tmp <- omnisummary %>% mutate(category=cut(omnisummary[,paste0(samp,".",mark,"_ratio10")], breaks=MyBreaks))
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/hist/pdf/",samp,".",mark,".uneven.percentage.pdf"),height=12,width=12)
            print(ggplot(tmp[!is.na(tmp$category),], aes(x = category)) +  
                geom_bar(aes(y = (..count..)/sum(..count..),fill = category)) + 
                geom_text(stat = "count", aes(label = after_stat(count)), position="fill", size = 12) +
                scale_y_continuous(labels=percent) +
                scale_x_discrete(drop=FALSE,labels = c("Paternal","Paternal bias","No bias","Maternal bias","Maternal")) +
                scale_fill_manual(values = c("(-1e-07,0.05]"="#69c7bd","(0.05,0.35]"="#5ecbe8","(0.35,0.65]"="#bca9d2","(0.65,0.95]"="#c572ae","(0.95,1]"="#f175ad"),guide = "none") +
                theme_classic() +
                labs(x = "Transcription maternal ratio", y = "Percentage of resolved genes") +
                theme(axis.text.x = element_text(color = "black",size = 36,hjust = 1,angle = 45),
                      axis.text.y = element_text(color = "black",size = 36),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(color = "black",size = 36),
                      legend.title = element_text(color = "black",size = 36),
                      legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/peak/",samp,".",mark,".chromratio.peak.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmppeak, y=!!tmpmark,color=!!tmppeak)) + 
                geom_hline(yintercept=0.5) + 
                geom_violin(position = position_dodge(0.9)) +
                geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA) +
                geom_point(alpha = 1/5, size = 0.5,position = position_jitterdodge(jitter.width = 0.8,dodge.width = 0.9)) +
                stat_compare_means(aes(group = !!tmppeak), label = "p.format") +
                lims( y = c(0,1.05)) +
                labs(x = paste0(mark," peak presence"), y=paste0(mark," maternal ratio")) +
                theme_classic() +
                scale_color_manual(name=paste0(mark," peak presence"),values=c("Yes"="#4682B4","No"="#B47846")) +
                theme(axis.text.x = element_text(color = "black",size = 36),
                      axis.text.y = element_text(color = "black",size = 36),
                      axis.title.x = element_text(color = "black",size = 36),
                      axis.title.y = element_text(color = "black",size = 36),
                      legend.title = element_text(color = "black",size = 36),
                      legend.text = element_text(colour = "black",size = 36))) 
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/",samp,".",mark,".TPM-chromenrich.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=asinh(!!tmpenrich), y=asinh(!!tmptpm))) + 
                geom_hline(yintercept=asinh(5)) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,6.05), y = c(-0.05, 12.5)) +
                labs(x = paste0("arcsinh(",mark," enrichment)"), y="arcsinh(TPM)") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/",samp,".",mark,".RNA.TPM-ratio.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpsamp, y=asinh(!!tmptpm))) + 
                geom_vline(xintercept=0.5) + 
                geom_hline(yintercept=asinh(5)) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 12.5)) +
                labs(x = "Transcription maternal ratio", y="arcsinh(TPM)") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/",samp,".",mark,".RNA.chromenrich.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpsamp, y=asinh(!!tmpenrich))) + 
                geom_vline(xintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 6.05)) +
                labs(x = "Transcription maternal ratio", y=paste0("arcsinh(",mark," enrichment)")) + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/",samp,".",mark,".TPM-ratio.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpmark, y=asinh(!!tmptpm))) + 
                geom_vline(xintercept=0.5) + 
                geom_hline(yintercept=asinh(5)) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 12.5)) +
                labs(x = paste0(mark," maternal ratio"), y="arcsinh(TPM)") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/",samp,".",mark,".ratio.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpmark, y=!!tmpsamp)) + 
                geom_vline(xintercept=0.5) + 
                geom_hline(yintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12,label.y.npc = 0) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 1.05)) +
                labs(x = paste0(mark," maternal ratio"), y="Transcription maternal ratio") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/",samp,".",mark,".ratio-chromenrich.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpmark, y=asinh(!!tmpenrich))) + 
                geom_vline(xintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 6.05)) +
                labs(x = paste0(mark," maternal ratio"), y=paste0("arcsinh(",mark," enrichment)")) + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
}

##Plot maternal ratios and scatterplots per cluster per mark per day
##Do it again for mutants with their respective stages
for (samp in c("WT13vsez16")){
    tmpsamp <- sym(paste0("13daf_ratio50"))
    tmpsampmut <- sym(paste0("ez23-16daf_ratio50"))
    tmptpm <- sym(paste0("13daf_tpmMean"))
    tmptpmmut <- sym(paste0("ez23-16daf_tpmMean"))
    for (mark in c("H3K9me1","H3K36me3","H3K27me3","H3")){
        tryCatch({
            tmpenrich <- sym(paste0("13daf.",mark,"_enrich"))
            tmpenrichmut <- sym(paste0("ez23-16daf.",mark,"_enrich"))
            tmpmark <- sym(paste0("13daf.",mark,"_ratio10"))
            tmpmarkmut <- sym(paste0("ez23-16daf.",mark,"_ratio10"))
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/enrich/",samp,".",mark,".enrich.pdf"),height=12,width=12)
            print(ggplot(melt(omnisummary[,c(paste0("13daf.",mark,"_enrich"),paste0("ez23-16daf.",mark,"_enrich"))]),aes(x=variable,y=asinh(value))) +
                geom_violin(position = position_dodge(0.9)) +
                geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA) +
                stat_compare_means(aes(group = variable), label = "p.format") +
                scale_y_continuous(name=paste0("asinh(enrichment)"),limits = c(0,6.5),breaks=c(0,2,4,6)) +
                scale_x_discrete(labels = c("13daf","ez23-16daf")) +
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                      axis.text.y = element_text(color = "black",size = 36),
                      axis.title.x = element_text(color = "black",size = 36),
                      axis.title.y = element_text(color = "black",size = 36),
                      legend.title = element_text(color = "black",size = 36),
                      legend.text = element_text(colour = "black",size = 36))) 
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".TPM-chromenrich.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=asinh(!!tmpenrich), y=asinh(!!tmptpmmut))) + 
                geom_hline(yintercept=asinh(5)) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,6.05), y = c(-0.05, 12.5)) +
                labs(x = paste0("arcsinh(WT ",mark," enrichment)"), y="arcsinh(Mutant TPM)") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".TPM-TPM.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=asinh(!!tmptpm), y=asinh(!!tmptpmmut))) + 
                geom_hline(yintercept=asinh(5)) + 
                geom_vline(xintercept=asinh(5)) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                geom_abline(slope = 1,intercept = c(0,0)) +
                geom_smooth() +
                lims(x = c(-0.05, 12.5), y = c(-0.05, 12.5)) +
                labs(x = "arcsinh(WT TPM)", y="arcsinh(Mutant TPM)") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".TPM-ratio.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpmark, y=asinh(!!tmptpmmut))) + 
                geom_vline(xintercept=0.5) + 
                geom_hline(yintercept=asinh(5)) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 12.5)) +
                labs(x = paste0("Wild type ",mark," maternal ratio"), y="arcsinh(Mutant TPM)") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".ratio.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpmark, y=!!tmpsampmut)) + 
                geom_vline(xintercept=0.5) + 
                geom_hline(yintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12,label.y.npc = 0) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 1.05)) +
                labs(x = paste0("Wild type ",mark," maternal ratio"), y="Mutant Transcription maternal ratio") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".WT-",samp,".ratio.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpsamp, y=!!tmpsampmut)) + 
                geom_vline(xintercept=0.5) + 
                geom_hline(yintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 1.05)) +
                labs(x = "WT Transcription maternal ratio", y="Mutant Transcription maternal ratio") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".RNA.chromenrich.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=asinh(!!tmpenrich), y=!!tmpsampmut)) + 
                geom_hline(yintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,6.05), y = c(-0.05, 1.05)) +
                labs(x = paste0("arcsinh(WT ",mark," enrichment)"), y="Mutant Transcription maternal ratio") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".chromenrich-chromenrich.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=asinh(!!tmpenrich), y=asinh(!!tmpenrichmut))) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                geom_abline(slope = 1,intercept = c(0,0)) +
                geom_smooth() +
                lims(x = c(-0.05, 6.05), y = c(-0.05, 6.05)) +
                labs(x = paste0("arcsinh(WT ",mark," enrichment)"), y=paste0("arcsinh(Mutant ",mark," enrichment)")) + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".WT-TPMratio.chromenrich.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpsamp, y=asinh(!!tmpenrichmut))) + 
                geom_vline(xintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 6.05)) +
                labs(x = "WT Transcription maternal ratio", y=paste0("arcsinh(Mutant ",mark," enrichment)")) + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".WT-ratio.chromenrich.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpmark, y=asinh(!!tmpenrichmut))) + 
                geom_vline(xintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 6.05)) +
                labs(x = paste0("Wild type ",mark," maternal ratio"), y=paste0("arcsinh(Mutant ",mark," enrichment)")) + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".WT-TPM.chromenrich.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=asinh(!!tmptpm), y=asinh(!!tmpenrichmut))) + 
                geom_vline(xintercept=asinh(5)) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05, 12.5), y = c(-0.05, 6.05)) +
                labs(x = "arcsinh(WT TPM)", y=paste0("arcsinh(Mutant ",mark," enrichment)")) + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".WT-chromenrich.chromratio.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=asinh(!!tmpenrich), y=!!tmpmarkmut)) + 
                geom_hline(yintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05, 6.05), y = c(-0.05, 1.05)) +
                labs(x = paste0("arcsinh(WT ",mark," enrichment)"), y=paste0("Mutant ",mark," maternal ratio")) + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".WT-TPMratio.chromratio.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpsamp, y=!!tmpmarkmut)) + 
                geom_vline(xintercept=0.5) + 
                geom_hline(yintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 1.05)) +
                labs(x = "WT Transcription maternal ratio", y=paste0("Mutant ",mark," maternal ratio")) + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".WT-ratio.chromratio.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=!!tmpmark, y=!!tmpmarkmut)) + 
                geom_vline(xintercept=0.5) + 
                geom_hline(yintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05,1.05), y = c(-0.05, 1.05)) +
                labs(x = paste0("Wild type ",mark," maternal ratio"), y=paste0("Mutant ",mark," maternal ratio")) + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/pdf/mutant-",samp,".",mark,".WT-TPM.chromratio.pdf"),height=12,width=12)
            print(ggplot(omnisummary, aes(x=asinh(!!tmptpm), y=!!tmpmarkmut)) + 
                geom_vline(xintercept=asinh(5)) + 
                geom_hline(yintercept=0.5) + 
                geom_point(alpha = 1/5, size = 0.5) +
                stat_cor(method = "spearman",size = 12) +
                stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(-0.05, 12.5), y = c(-0.05, 1.05)) +
                labs(x = "arcsinh(WT TPM)", y=paste0("Mutant ",mark," maternal ratio")) + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
                    pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/cluster/pdf/mutant-",samp,".",mark,".TPM-chromenrich.clust.pdf"),height=12,width=12)
                    print(ggplot(omnisummary %>% group_by(!!tmpcluster) %>% filter(n() >= 100), aes(x=asinh(!!tmpenrich), y=asinh(!!tmptpmmut))) + 
                        geom_hline(yintercept=asinh(5)) + 
                        geom_point(alpha = 1/5, size = 0.5) +
                        stat_cor(method = "spearman",size = 12) +
                        stat_density_2d(aes(fill = ..level..), geom="polygon") +
                        facet_wrap(vars(!!tmpcluster),ncol=3) + 
                        scale_fill_continuous(type = "viridis",guide="none") +
                        lims(x = c(-0.05,6.05), y = c(-0.05, 12.5)) +
                        labs(x = paste0("arcsinh(",mark," enrichment)"), y="arcsinh(Mutant TPM)") + 
                        theme_classic() +
                        theme(axis.text.x = element_text(color = "black",size = 36),
                          axis.text.y = element_text(color = "black",size = 36),
                          axis.title.x = element_text(color = "black",size = 36),
                          axis.title.y = element_text(color = "black",size = 36),
                          legend.title = element_text(color = "black",size = 36),
                          legend.text = element_text(colour = "black",size = 36)))
                    dev.off()
                    pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/cluster/pdf/mutant-",samp,".",mark,".TPM-TPM.clust.pdf"),height=12,width=12)
                    print(ggplot(omnisummary %>% group_by(!!tmpcluster) %>% filter(n() >= 100), aes(x=asinh(!!tmptpm), y=asinh(!!tmptpmmut))) + 
                        geom_hline(yintercept=asinh(5)) + 
                        geom_vline(xintercept=asinh(5)) + 
                        geom_point(alpha = 1/5, size = 0.5) +
                        stat_cor(method = "spearman",size = 12) +
                        stat_density_2d(aes(fill = ..level..), geom="polygon") +
                        facet_wrap(vars(!!tmpcluster),ncol=3) + 
                        scale_fill_continuous(type = "viridis",guide="none") +
                        geom_abline(slope = 1,intercept = c(0,0)) +
                        geom_smooth() +
                        lims(x = c(-0.05, 12.5), y = c(-0.05, 12.5)) +
                        labs(x = "arcsinh(WT TPM)", y="arcsinh(Mutant TPM)") + 
                        theme_classic() +
                        theme(axis.text.x = element_text(color = "black",size = 36),
                          axis.text.y = element_text(color = "black",size = 36),
                          axis.title.x = element_text(color = "black",size = 36),
                          axis.title.y = element_text(color = "black",size = 36),
                          legend.title = element_text(color = "black",size = 36),
                          legend.text = element_text(colour = "black",size = 36)))
                    dev.off()
                    pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/cluster/pdf/mutant-",samp,".",mark,".TPM-ratio.clust.pdf"),height=12,width=12)
                    print(ggplot(omnisummary %>% group_by(!!tmpcluster) %>% filter(n() >= 100), aes(x=!!tmpmark, y=asinh(!!tmptpmmut))) + 
                        geom_vline(xintercept=0.5) + 
                        geom_hline(yintercept=asinh(5)) + 
                        geom_point(alpha = 1/5, size = 0.5) +
                        stat_cor(method = "spearman",size = 12) +
                        stat_density_2d(aes(fill = ..level..), geom="polygon") +
                        facet_wrap(vars(!!tmpcluster),ncol=3) + 
                        scale_fill_continuous(type = "viridis",guide="none") +
                        lims(x = c(-0.05,1.05), y = c(-0.05, 12.5)) +
                        labs(x = paste0(mark," maternal ratio"), y="arcsinh(Mutant TPM)") + 
                        theme_classic() +
                        theme(axis.text.x = element_text(color = "black",size = 36),
                          axis.text.y = element_text(color = "black",size = 36),
                          axis.title.x = element_text(color = "black",size = 36),
                          axis.title.y = element_text(color = "black",size = 36),
                          legend.title = element_text(color = "black",size = 36),
                          legend.text = element_text(colour = "black",size = 36)))
                    dev.off()
                    pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/cluster/pdf/mutant-",samp,".",mark,".RNA.TPM-ratio.clust.pdf"),height=12,width=12)
                    print(ggplot(omnisummary %>% group_by(!!tmpcluster) %>% filter(n() >= 100), aes(x=!!tmpsampmut, y=asinh(!!tmptpmmut))) + 
                        geom_vline(xintercept=0.5) + 
                        geom_hline(yintercept=asinh(5)) + 
                        geom_point(alpha = 1/5, size = 0.5) +
                        stat_cor(method = "spearman",size = 12) +
                        stat_density_2d(aes(fill = ..level..), geom="polygon") +
                        facet_wrap(vars(!!tmpcluster),ncol=3) + 
                        scale_fill_continuous(type = "viridis",guide="none") +
                        lims(x = c(-0.05,1.05), y = c(-0.05, 12.5)) +
                        labs(x = "Mutant Transcription maternal ratio", y="arcsinh(Mutant TPM)") + 
                        theme_classic() +
                        theme(axis.text.x = element_text(color = "black",size = 36),
                          axis.text.y = element_text(color = "black",size = 36),
                          axis.title.x = element_text(color = "black",size = 36),
                          axis.title.y = element_text(color = "black",size = 36),
                          legend.title = element_text(color = "black",size = 36),
                          legend.text = element_text(colour = "black",size = 36)))
                    dev.off()
                    pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/cluster/pdf/mutant-",samp,".",mark,".ratio.clust.pdf"),height=12,width=12)
                    print(ggplot(omnisummary %>% group_by(!!tmpcluster) %>% filter(n() >= 100), aes(x=!!tmpmark, y=!!tmpsampmut)) + 
                        geom_vline(xintercept=0.5) + 
                        geom_hline(yintercept=0.5) + 
                        geom_point(alpha = 1/5, size = 0.5) +
                        stat_cor(method = "spearman",size = 12,label.y.npc = 0) +
                        stat_density_2d(aes(fill = ..level..), geom="polygon") +
                        facet_wrap(vars(!!tmpcluster),ncol=3) + 
                        scale_fill_continuous(type = "viridis",guide="none") +
                        lims(x = c(-0.05,1.05), y = c(-0.05, 1.05)) +
                        labs(x = paste0(mark," maternal ratio"), y="Mutant Transcription maternal ratio") + 
                        theme_classic() +
                        theme(axis.text.x = element_text(color = "black",size = 36),
                          axis.text.y = element_text(color = "black",size = 36),
                          axis.title.x = element_text(color = "black",size = 36),
                          axis.title.y = element_text(color = "black",size = 36),
                          legend.title = element_text(color = "black",size = 36),
                          legend.text = element_text(colour = "black",size = 36)))
                    dev.off()
                    pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/cluster/pdf/mutant-",samp,".WT-",samp,".ratio.clust.pdf"),height=12,width=12)
                    print(ggplot(omnisummary %>% group_by(!!tmpcluster) %>% filter(n() >= 100), aes(x=!!tmpsamp, y=!!tmpsampmut)) + 
                        geom_vline(xintercept=0.5) + 
                        geom_hline(yintercept=0.5) + 
                        geom_point(alpha = 1/5, size = 0.5) +
                        stat_cor(method = "spearman",size = 12) +
                        stat_density_2d(aes(fill = ..level..), geom="polygon") +
                        facet_wrap(vars(!!tmpcluster),ncol=3) + 
                        scale_fill_continuous(type = "viridis",guide="none") +
                        lims(x = c(-0.05,1.05), y = c(-0.05, 1.05)) +
                        labs(x = "WT Transcription maternal ratio", y="Mutant Transcription maternal ratio") + 
                        theme_classic() +
                        theme(axis.text.x = element_text(color = "black",size = 36),
                          axis.text.y = element_text(color = "black",size = 36),
                          axis.title.x = element_text(color = "black",size = 36),
                          axis.title.y = element_text(color = "black",size = 36),
                          legend.title = element_text(color = "black",size = 36),
                          legend.text = element_text(colour = "black",size = 36)))
                    dev.off()
                    pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/cluster/pdf/mutant-",samp,".",mark,".RNA.chromenrich.clust.pdf"),height=12,width=12)
                    print(ggplot(omnisummary %>% group_by(!!tmpcluster) %>% filter(n() >= 100), aes(x=!!tmpsampmut, y=asinh(!!tmpenrich))) + 
                        geom_vline(xintercept=0.5) + 
                        geom_point(alpha = 1/5, size = 0.5) +
                        stat_cor(method = "spearman",size = 12) +
                        stat_density_2d(aes(fill = ..level..), geom="polygon") +
                        facet_wrap(vars(!!tmpcluster),ncol=3) + 
                        scale_fill_continuous(type = "viridis",guide="none") +
                        lims(x = c(-0.05,1.05), y = c(-0.05, 6.05)) +
                        labs(x = "Mutant Transcription maternal ratio", y=paste0("arcsinh(",mark," enrichment)")) + 
                        theme_classic() +
                        theme(axis.text.x = element_text(color = "black",size = 36),
                          axis.text.y = element_text(color = "black",size = 36),
                          axis.title.x = element_text(color = "black",size = 36),
                          axis.title.y = element_text(color = "black",size = 36),
                          legend.title = element_text(color = "black",size = 36),
                          legend.text = element_text(colour = "black",size = 36)))
                    dev.off()
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
}



################################################################################################################################################
##Redo scatterplots and violin plots for paper
################################################################################################################################################

for (samp in c("13daf","ez23-16daf")){
    tmpcluster <- sym(paste0(samp,"_cluster"))
    tmpsamp <- sym(paste0(samp,"_ratio50"))
    tmptpm <- sym(paste0(samp,"_tpmMean"))
    for (mark in c("H3K27me3")){
        tryCatch({
            tmpenrich <- sym(paste0(samp,".",mark,"_enrich"))
            tmpmark <- sym(paste0(samp,".",mark,"_ratio10"))
            png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/png/",samp,".",mark,".TPM-ratio.png"), type="cairo",width = 1024, height = 1024, units = "px")
            print(ggplot(omnisummary, aes(x=!!tmpmark, y=asinh(!!tmptpm))) + 
                geom_vline(xintercept=0.5) + 
                geom_hline(yintercept=asinh(5)) + 
                geom_point(alpha = 1/5, size = 2) +
                stat_cor(method = "spearman",size = 12,aes(label = ..r.label..)) +
                # stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                # scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(0,1), y = c(0, 12.5)) +
                labs(x = paste0(mark," maternal ratio"), y="arcsinh(TPM)") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/png/",samp,".",mark,".TPM-chromenrich.png"), type="cairo",width = 1024, height = 1024, units = "px")
            print(omnisummary %>%
                dplyr::mutate(
                    out_of_bounds = asinh(!!tmpenrich) > 2,
                    tmpenrich = ifelse(out_of_bounds, sinh(2), !!tmpenrich)
                    ) %>%
                ggplot() +
                    geom_hline(yintercept=asinh(5)) + 
                    geom_point(aes(x = asinh(!!tmpenrich), y = asinh(!!tmptpm), shape = out_of_bounds),alpha = 1/5, size = 2,show.legend = F) +
                    stat_cor(method = "spearman",size = 12,aes(x = asinh(!!tmpenrich), y = asinh(!!tmptpm),label = ..r.label..)) +
                    # stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                    # scale_fill_continuous(type = "viridis",guide="none") +
                    lims(x = c(0,2.05), y = c(0, 12.5)) +
                    labs(x = paste0("arcsinh(",mark," enrichment)"), y="arcsinh(TPM)") + 
                    theme_classic() +
                    theme(axis.text.x = element_text(color = "black",size = 36),
                        axis.text.y = element_text(color = "black",size = 36),
                        axis.title.x = element_text(color = "black",size = 36),
                        axis.title.y = element_text(color = "black",size = 36),
                        legend.title = element_text(color = "black",size = 36),
                        legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/png/",samp,".",mark,".RNA.TPM-ratio.png"), type="cairo",width = 1024, height = 1024, units = "px")
            print(ggplot(omnisummary, aes(x=!!tmpsamp, y=asinh(!!tmptpm))) + 
                geom_vline(xintercept=0.5) + 
                geom_hline(yintercept=asinh(5)) + 
                geom_point(alpha = 1/5, size = 2) +
                stat_cor(method = "spearman",size = 12,aes(label = ..r.label..)) +
                # stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                # scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(0,1), y = c(0, 12.5)) +
                labs(x = "Transcription maternal ratio", y="arcsinh(TPM)") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/png/",samp,".",mark,".RNA.chromenrich.png"), type="cairo",width = 1024, height = 1024, units = "px")
            print(omnisummary %>%
                dplyr::mutate(
                    out_of_bounds = asinh(!!tmpenrich) > 2,
                    tmpenrich = ifelse(out_of_bounds, sinh(2), !!tmpenrich)
                ) %>%
                ggplot() +
                geom_vline(xintercept=0.5) + 
                geom_point(aes(x = !!tmpsamp, y = asinh(!!tmpenrich), shape = out_of_bounds),alpha = 1/5, size = 2,show.legend = F) +
                stat_cor(method = "spearman",size = 12,aes(x = !!tmpsamp, y = asinh(!!tmpenrich),label = ..r.label..)) +
                # stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                # scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(0,1), y = c(0, 2.05)) +
                labs(x = "Transcription maternal ratio", y=paste0("arcsinh(",mark," enrichment)")) + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                      axis.text.y = element_text(color = "black",size = 36),
                      axis.title.x = element_text(color = "black",size = 36),
                      axis.title.y = element_text(color = "black",size = 36),
                      legend.title = element_text(color = "black",size = 36),
                      legend.text = element_text(colour = "black",size = 36)))
            dev.off()
            png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/scatter/png/",samp,".",mark,".ratio.png"), type="cairo",width = 1024, height = 1024, units = "px")
            print(ggplot(omnisummary, aes(x=!!tmpmark, y=!!tmpsamp)) + 
                geom_vline(xintercept=0.5) + 
                geom_hline(yintercept=0.5) + 
                geom_point(alpha = 1/2, size = 2) +
                stat_cor(method = "spearman",size = 12,label.y.npc = 0,aes(label = ..r.label..)) +
                # stat_density_2d(aes(fill = ..level..), geom="polygon") + 
                # scale_fill_continuous(type = "viridis",guide="none") +
                lims(x = c(0,1), y = c(0,1)) +
                labs(x = paste0(mark," maternal ratio"), y="Transcription maternal ratio") + 
                theme_classic() +
                theme(axis.text.x = element_text(color = "black",size = 36),
                  axis.text.y = element_text(color = "black",size = 36),
                  axis.title.x = element_text(color = "black",size = 36),
                  axis.title.y = element_text(color = "black",size = 36),
                  legend.title = element_text(color = "black",size = 36),
                  legend.text = element_text(colour = "black",size = 36)))
            dev.off()
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
}


pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/matratio-comparison.violin.pdf"),height=12,width=12)
ggplot(melt(omnisummary, measure.vars = c("13daf.H3K27me3_ratio10","ez23-16daf.H3K27me3_ratio10","13daf_ratio50","ez23-16daf_ratio50")),aes(x=variable,y=value)) +
    geom_hline(yintercept = 0.5) +
    geom_violin(position = position_dodge(0.9),size=1) +
    geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA,size=1) +
    scale_y_continuous(name="Maternal ratio",breaks=c(0,0.25,0.5,0.75,1)) +
    scale_x_discrete(labels=c("WT H3K27me3","Mutant H3K27me3","WT Transcription","Mutant Transcription")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36,hjust = 1,angle=45),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()

##Stats and Cohen's D
wilcox.test(na.omit(omnisummary[,c("13daf.H3K27me3_enrich","ez23-16daf.H3K27me3_enrich")])$`13daf.H3K27me3_enrich`, na.omit(omnisummary[,c("13daf.H3K27me3_enrich","ez23-16daf.H3K27me3_enrich")])$`ez23-16daf.H3K27me3_enrich`, paired=TRUE)
cohen.d(omnisummary$`ez23-16daf.H3K27me3_ratio10`,omnisummary$`13daf.H3K27me3_ratio10`,na.rm = T)
wilcox.test(na.omit(omnisummary[,c("13daf.H3K27me3_ratio10","ez23-16daf.H3K27me3_ratio10")])$`13daf.H3K27me3_ratio10`, na.omit(omnisummary[,c("13daf.H3K27me3_ratio10","ez23-16daf.H3K27me3_ratio10")])$`ez23-16daf.H3K27me3_ratio10`, paired=TRUE)
cohen.d(omnisummary$`ez23-16daf_ratio50`,omnisummary$`13daf_ratio50`,na.rm = T)
wilcox.test(na.omit(omnisummary[,c("13daf_ratio50","ez23-16daf_ratio50")])$`13daf_ratio50`, na.omit(omnisummary[,c("13daf_ratio50","ez23-16daf_ratio50")])$`ez23-16daf_ratio50`, paired=TRUE)


################################################################################################################################################
##Next section of plots
################################################################################################################################################

cols <- c("RNASeq"="#949494","H3K9me1"="#CB181D","H3K36me3"="#A1D99B","H2AZ"="#4292C6","H2AK119ub"="#6BAED6","H3K27me3"="#C6DBEF","H3"="#cdcd00","H3K27me1"="#FB6A4A","H3K9ac"="#005A32","H3K14ac"="#238B45","H3K4me1"="#74C476","H3K4me3"="#084594","IgG"="#860e6a")


################################################################################################################################################
##Heatmap of specific gene expression
################################################################################################################################################

c2 <- asinh(tpm[,c("Tak1_tpmMean","Cam2_tpmMean","an_tpmMean","ar_tpmMean","13daf_tpmMean")])
c3 <- c2[row.names(c2) %in% c("Mp4g21020","Mp3g06080","Mp5g12760","Mp5g18040","Mp5g17740","Mp5g22900"),]
my_palette <- colorRampPalette(brewer.pal(3, "YlOrRd"))
pdf(file="/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/plots/PRC2.expression.heatmap.pdf", width=10)
pheatmap(c3, 
         show_rownames=T, 
         color=c("white",my_palette(length(c(0,seq(1,max(c3)+0.01,by=0.01))))), 
         breaks = c(0,seq(1,max(c3)+0.01,by=0.01)), 
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 18, 
         labels_row = c("FIE","MSI1","Su(z)12","E(z)2","E(z)1","E(z)3"),
         labels_col = c("Male vegetative","Female vegetative","Male sexual organ","Female sexual organ","Wild-type embryo"),
         treeheight_row = 0,
         angle_col = "45")
dev.off()


#####################################################################################################################################################################
##Plot maternal ratio depending on if a gene is sporophyte-specific or not
#####################################################################################################################################################################
pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/plots/matratio/embryo-specific.matratio.13daf.pdf"),height=12,width=12)
ggplot(omnisummary,aes(x=(row.names(omnisummary) %in% row.names(subset(omnisummary,Tak1_tpmMean<5&Cam2_tpmMean<5&`13daf_tpmMean`>5&g_tpmMean<5&ar_tpmMean<5&an_tpmMean<5))),y=`13daf_ratio50`)) +
    geom_hline(yintercept=0.5) + 
    geom_violin(position = position_dodge(0.9)) +
    geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA) +
    geom_jitter(alpha=1/20, size = 0.5) +
    stat_compare_means() +
    lims( y = c(0,1.05)) +
    labs(x = paste0("Embryo-specific gene"), y=paste0("Maternal ratio")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)) 
dev.off()

##H3K27me3 maternal ratio of embryo-specific genes at 13daf
ggplot(omnisummary,aes(x=(row.names(omnisummary) %in% row.names(subset(omnisummary,Tak1_tpmMean<5&Cam2_tpmMean<5&`13daf_tpmMean`>5&g_tpmMean<5&ar_tpmMean<5&an_tpmMean<5))),y=`13daf.H3K27me3_ratio10`)) +
    geom_hline(yintercept=0.5) + 
    geom_violin(position = position_dodge(0.9)) +
    geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA) +
    geom_jitter(alpha=1/20, size = 0.5) +
    stat_compare_means() +
    lims( y = c(0,1.05)) +
    labs(x = paste0("Embryo-specific gene"), y=paste0("Maternal ratio")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)) 

##H3K27me3 enrichment of embryo-specific genes at 13daf
ggplot(omnisummary,aes(x=(row.names(omnisummary) %in% row.names(subset(omnisummary,Tak1_tpmMean<5&Cam2_tpmMean<5&`13daf_tpmMean`>5&g_tpmMean<5&ar_tpmMean<5&an_tpmMean<5))),y=asinh(`13daf.H3K27me3_enrich`))) +
    geom_hline(yintercept=0.5) + 
    geom_violin(position = position_dodge(0.9)) +
    geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA) +
    geom_jitter(alpha=1/20, size = 0.5) +
    stat_compare_means() +
    lims( y = c(0,6)) +
    labs(x = paste0("Embryo-specific gene"), y=paste0("Maternal ratio")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)) 


################################################################################################################################################
##Plot ratio of maternal ratio between mutant and WT
################################################################################################################################################

ggplot(melt(data.frame(Transcription=omnisummary$`ez23-16daf_ratio50`/omnisummary$`13daf_ratio50`)),aes(x=variable,y=log2(value))) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    geom_hline(yintercept = log2(1)) +
    scale_x_discrete(name=NA) +
    scale_y_continuous(name="Mutant/WT H3K27me3 maternal ratio",limits = c(-2,2)) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))

ggplot(melt(data.frame(H3K27me3=omnisummary$`ez23-16daf.H3K27me3_ratio10`/omnisummary$`13daf.H3K27me3_ratio10`)),aes(x=variable,y=log2(value))) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    geom_hline(yintercept = log2(1)) +
    scale_x_discrete(name=NA) +
    scale_y_continuous(name="Mutant/WT H3K27me3 maternal ratio",limits = c(-5,5)) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))

pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/Mutant-WT.enrich.ratio.pdf"),height=12,width=12)
ggplot(melt(data.frame(H3K27me3=omnisummary$`ez23-16daf.H3K27me3_enrich`/omnisummary$`13daf.H3K27me3_enrich`)),aes(x=variable,y=log2(value))) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    geom_hline(yintercept = log2(1)) +
    scale_x_discrete(name=NA) +
    scale_y_continuous(name="log2(Mutant/WT H3K27me3 enrichment)",limits = c(-4,4)) +
    theme_classic() +
    theme(axis.text.x = element_blank(),,
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()

pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/Mutant-WT.matratio.ratio.pdf"),height=5,width=5)
ggplot(melt(data.frame(H3K27me3=omnisummary$`ez23-16daf.H3K27me3_ratio10` - omnisummary$`13daf.H3K27me3_ratio10`,Transcription=omnisummary$`ez23-16daf_ratio50` - omnisummary$`13daf_ratio50`)),aes(x=variable,y=value)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    geom_hline(yintercept = 0) +
    scale_x_discrete(name=NA) +
    scale_y_continuous(name="Mutant - WT maternal ratio",limits = c(-0.7,0.7),breaks = c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()




################################################################################################################################################
##Plot correlation matrix of WT vs Mutant traits
################################################################################################################################################

tmp <- round(cor(omnisummary[,c("13daf_ratio50","13daf.H3K27me3_ratio10","13daf.H3K27me3_enrich","ez23-16daf_ratio50","ez23-16daf.H3K27me3_ratio10","ez23-16daf.H3K27me3_enrich")], use = "pairwise.complete.obs",method = "spearman"),2)
# Get lower triangle of the correlation matrix
get_lower_tri<-function(tmp){
    tmp[upper.tri(tmp)] <- NA
    return(tmp)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(tmp){
    tmp[lower.tri(tmp)]<- NA
    return(tmp)
}
upper_tri <- get_upper_tri(tmp)
tmp2 <- melt(upper_tri, na.rm = TRUE)
pdf(paste0(file="/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/correlation.spearman.heatmap.pdf"),height=5.5,width=5.5)
ggplot(tmp2, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") +
    scale_x_discrete(labels = c("WT transcription maternal ratio","WT H3K27me3 maternal ratio","WT H3K27me3 enrichment","Mutant transcription maternal ratio","Mutant H3K27me3 maternal ratio","Mutant H3K27me3 enrichment")) +
    scale_y_discrete(labels = c("WT transcription maternal ratio","WT H3K27me3 maternal ratio","WT H3K27me3 enrichment","Mutant transcription maternal ratio","Mutant H3K27me3 maternal ratio","Mutant H3K27me3 enrichment")) +
    theme_classic()+ # minimal theme
    coord_fixed() + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
    theme(axis.text.x = element_text(color = "black",angle = 45, vjust = 1,size = 12, hjust = 1),
        axis.text.y = element_text(color = "black",size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(color = "black",size = 12),
        legend.text = element_text(colour = "black",size = 12),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.7),
        legend.direction = "horizontal") +
    guides(fill = guide_colorbar(barwidth = 7, 
        barheight = 1,
        title.position = "top", 
        title.hjust = 0.5))
dev.off()

tmp2 <- melt(upper_tri[c(1,2,3),c(4,5,6)], na.rm = TRUE)
pdf(paste0(file="/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/correlation.spearman.heatmap.WT-mutant.pdf"),height=5.5,width=5.5)
ggplot(tmp2, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") +
    scale_x_discrete(labels = c("Mutant transcription maternal ratio","Mutant H3K27me3 maternal ratio","Mutant H3K27me3 enrichment")) +
    scale_y_discrete(labels = c("WT transcription maternal ratio","WT H3K27me3 maternal ratio","WT H3K27me3 enrichment")) +
    theme_classic()+ # minimal theme
    coord_fixed() + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
    theme(axis.text.x = element_text(color = "black",angle = 45, vjust = 1,size = 12, hjust = 1),
          axis.text.y = element_text(color = "black",size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(color = "black",size = 12),
          legend.text = element_text(colour = "black",size = 12),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank()) +
    guides(fill = guide_colorbar(title.position = "top", 
                                 title.hjust = 0.5))
dev.off()

tmp <- round(cor(omnisummary[,c("ez23-16daf_ratio50","ez23-16daf.H3K27me3_ratio10","ez23-16daf.H3K27me3_enrich","ez23-16daf_tpmMean")], use = "pairwise.complete.obs",method = "spearman"),2)
# Get lower triangle of the correlation matrix
get_lower_tri<-function(tmp){
    tmp[upper.tri(tmp)] <- NA
    return(tmp)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(tmp){
    tmp[lower.tri(tmp)]<- NA
    return(tmp)
}
upper_tri <- get_upper_tri(tmp)
tmp2 <- melt(upper_tri, na.rm = TRUE)
pdf(paste0(file="/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/correlation.spearman.heatmap.mutant.pdf"),height=5.5,width=5.5)
ggplot(tmp2, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") +
    scale_x_discrete(labels = c("Mutant transcription maternal ratio","Mutant H3K27me3 maternal ratio","Mutant H3K27me3 enrichment","Mutant TPM")) +
    scale_y_discrete(labels = c("Mutant transcription maternal ratio","Mutant H3K27me3 maternal ratio","Mutant H3K27me3 enrichment","Mutant TPM")) +
    theme_classic()+ # minimal theme
    coord_fixed() + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
    theme(axis.text.x = element_text(color = "black",angle = 45, vjust = 1,size = 12, hjust = 1),
        axis.text.y = element_text(color = "black",size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(color = "black",size = 12),
        legend.text = element_text(colour = "black",size = 12),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.7),
        legend.direction = "horizontal") +
    guides(fill = guide_colorbar(barwidth = 7, 
        barheight = 1,
        title.position = "top", 
        title.hjust = 0.5))
dev.off()


################################################################################################################################################
##Plot maternal ratios of WT and mutant depending on maternal ratio category of WT
################################################################################################################################################
tmp <- omnisummary

tmp$`WT transcription maternal bias category` <- NA
# tmp[row.names(subset(tmp,`13daf_ratio50` <= 0.05)),]$`WT transcription maternal bias category` <- c("Fully paternal")
# tmp[row.names(subset(tmp,`13daf_ratio50` > 0.05 & `13daf_ratio50` <= 0.35)),]$`WT transcription maternal bias category` <- c("Paternal bias")
tmp[row.names(subset(tmp,`13daf_ratio50` < 0.65 & `13daf_ratio50` > 0.35)),]$`WT transcription maternal bias category` <- c("No bias")
tmp[row.names(subset(tmp,`13daf_ratio50` >= 0.65 & `13daf_ratio50` < 0.95)),]$`WT transcription maternal bias category` <- c("Maternal bias")
tmp[row.names(subset(tmp,`13daf_ratio50` >= 0.95)),]$`WT transcription maternal bias category` <- c("Fully maternal")
tmp$`WT transcription maternal bias category` <- factor(tmp$`WT transcription maternal bias category`, levels = c("No bias","Maternal bias","Fully maternal"))

tmp$`WT H3K27me3 maternal bias category` <- NA
tmp[row.names(subset(tmp,`13daf.H3K27me3_ratio10` <= 0.05)),]$`WT H3K27me3 maternal bias category` <- c("Fully paternal")
tmp[row.names(subset(tmp,`13daf.H3K27me3_ratio10` > 0.05 & `13daf.H3K27me3_ratio10` <= 0.35)),]$`WT H3K27me3 maternal bias category` <- c("Paternal bias")
tmp[row.names(subset(tmp,`13daf.H3K27me3_ratio10` < 0.65 & `13daf.H3K27me3_ratio10` > 0.35)),]$`WT H3K27me3 maternal bias category` <- c("No bias")
# tmp[row.names(subset(tmp,`13daf.H3K27me3_ratio10` >= 0.65 & `13daf.H3K27me3_ratio10` < 0.95)),]$`WT H3K27me3 maternal bias category` <- c("Maternal bias")
# tmp[row.names(subset(tmp,`13daf.H3K27me3_ratio10` >= 0.95)),]$`WT H3K27me3 maternal bias category` <- c("Fully maternal")
tmp$`WT H3K27me3 maternal bias category` <- factor(tmp$`WT H3K27me3 maternal bias category`, levels = c("Fully paternal", "Paternal bias", "No bias"))

# pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/Mutant-WT.TPMmatratio.TPMcategories.pdf"),height=12,width=12)
png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/Mutant-WT.TPMmatratio.TPMcategories.png"), type="cairo",width = 1024, height = 1024, units = "px")
ggplot(melt(tmp[!is.na(tmp$`WT transcription maternal bias category`),],id.vars = "WT transcription maternal bias category",measure.vars = c("13daf_ratio50","ez23-16daf_ratio50")), aes(x=`WT transcription maternal bias category`, y=value,color=variable)) + 
    geom_hline(yintercept=0.5) +  
    geom_violin(position = position_dodge(0.9),scale = "width") +
    geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA) +
    geom_point(alpha = 1/10, size = 0.5,position = position_jitterdodge(jitter.width = 0.8,dodge.width = 0.9)) +
    scale_y_continuous(name=paste0("Transcription maternal ratio"),limits = c(0,1.2),breaks=c(0,0.25,0.5,0.75,1)) +
    scale_color_manual(values=c("#000000","#FF0000"),labels=c("WT","Mutant")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36, angle = 45,hjust = 1),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_blank(),
          legend.text = element_text(colour = "black",size = 36))
dev.off()

cohen.d(subset(tmp,`WT transcription maternal bias category`=="Fully maternal")$`ez23-16daf_ratio50`,subset(tmp,`WT transcription maternal bias category`=="Fully maternal")$`13daf_ratio50`,na.rm = T)
cohen.d(subset(tmp,`WT transcription maternal bias category`=="Maternal bias")$`ez23-16daf_ratio50`,subset(tmp,`WT transcription maternal bias category`=="Maternal bias")$`13daf_ratio50`,na.rm = T)
cohen.d(subset(tmp,`WT transcription maternal bias category`=="No bias")$`ez23-16daf_ratio50`,subset(tmp,`WT transcription maternal bias category`=="No bias")$`13daf_ratio50`,na.rm = T)

# pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/Mutant-WT.H3K27me3matratio.TPMcategories.pdf"),height=12,width=12)
png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/Mutant-WT.H3K27me3matratio.TPMcategories.png"), type="cairo",width = 1024, height = 1024, units = "px")
ggplot(melt(tmp[!is.na(tmp$`WT transcription maternal bias category`),],id.vars = "WT transcription maternal bias category",measure.vars = c("13daf.H3K27me3_ratio10","ez23-16daf.H3K27me3_ratio10")), aes(x=`WT transcription maternal bias category`, y=value,color=variable)) + 
    geom_hline(yintercept=0.5) +  
    geom_violin(position = position_dodge(0.9),scale = "width") +
    geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA) +
    geom_point(alpha = 1/10, size = 0.5,position = position_jitterdodge(jitter.width = 0.8,dodge.width = 0.9)) +
    scale_y_continuous(name=paste0("H3K27me3 maternal ratio"),limits = c(0,1.2),breaks=c(0,0.25,0.5,0.75,1)) +
    scale_color_manual(values=c("#000000","#FF0000"),labels=c("WT","Mutant")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36, angle = 45,hjust = 1),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_blank(),
          legend.text = element_text(colour = "black",size = 36))
dev.off()

# pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/Mutant-WT.H3K27me3enrich.TPMcategories.pdf"),height=12,width=12)
png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/Mutant-WT.H3K27me3enrich.TPMcategories.png"), type="cairo",width = 1024, height = 1024, units = "px")
ggplot(melt(tmp[!is.na(tmp$`WT transcription maternal bias category`),],id.vars = "WT transcription maternal bias category",measure.vars = c("13daf.H3K27me3_enrich","ez23-16daf.H3K27me3_enrich")), aes(x=`WT transcription maternal bias category`, y=value,color=variable)) + 
    geom_violin(position = position_dodge(0.9),scale = "width") +
    geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA) +
    geom_point(alpha = 1/10, size = 0.5,position = position_jitterdodge(jitter.width = 0.8,dodge.width = 0.9)) +
    scale_y_continuous(name=paste0("arcsinh(H3K27me3 enrichment"),limits = c(0,6),breaks=c(0,2,4,6)) +
    scale_color_manual(values=c("#000000","#FF0000"),labels=c("WT","Mutant")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36, angle = 45,hjust = 1),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_blank(),
          legend.text = element_text(colour = "black",size = 36))
dev.off()

# pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/Mutant-WT.H3K27me3matratio.H3K27me3categories.pdf"),height=12,width=12)
png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/Mutant-WT.H3K27me3matratio.H3K27me3categories.png"), type="cairo",width = 1024, height = 1024, units = "px")
ggplot(melt(tmp[!is.na(tmp$`WT H3K27me3 maternal bias category`),],id.vars = "WT H3K27me3 maternal bias category",measure.vars = c("13daf.H3K27me3_ratio10","ez23-16daf.H3K27me3_ratio10")), aes(x=`WT H3K27me3 maternal bias category`, y=value,color=variable)) + 
    geom_hline(yintercept=0.5) +  
    geom_violin(position = position_dodge(0.9),scale = "width") +
    geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA) +
    geom_point(alpha = 1/10, size = 0.5,position = position_jitterdodge(jitter.width = 0.8,dodge.width = 0.9)) +
    scale_y_continuous(name=paste0("H3K27me3 maternal ratio"),limits = c(0,1.2),breaks=c(0,0.25,0.5,0.75,1)) +
    scale_color_manual(values=c("#000000","#FF0000"),labels=c("WT","Mutant")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36, angle = 45,hjust = 1),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_blank(),
          legend.text = element_text(colour = "black",size = 36))
dev.off()


################################################################################################################################################
##Plot maternal ratios of each gene along length of chromosome per chromosome
################################################################################################################################################

# pdf(paste0(file="/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/plots/matratio.pergene.chr.pdf"),height=5.5,width=5.5)
png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/plots/matratio.pergene.chr.transcription.png"), type="cairo",width = 4096, height = 2048, units = "px",res = 450)
ggplot(omnisummary, aes(x=(((end - start)/2)+start)/1000000, y=`13daf_ratio50`, color = chr)) +
    geom_point(size=2/3) +
    geom_hline(yintercept = 0.5) +
    geom_vline(data = chromsizes, xintercept = chromsizes[1:8,]$V2/1000000, color = c("#01B4C6","#97ce4c","#FFF874","#BEE5FD","#F675DA","#44281d","#F8D3AC","#E64358")) +
    scale_colour_manual(values=c("chr1"="#01B4C6","chr2"="#97ce4c","chr3"="#FFF874","chr4"="#BEE5FD","chr5"="#F675DA","chr6"="#44281d","chr7"="#F8D3AC","chr8"="#E64358"),name="Chromosome") +
    scale_y_continuous(name = "Transcription maternal ratio",limits = c(0,1)) +
    scale_x_continuous(name = "Position on chromosome (Mb)") +
    guides(colour = guide_legend(override.aes = list(size=5,alpha=1))) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 20),
          axis.text.y = element_text(color = "black",size = 20),
          axis.title.x = element_text(color = "black",size = 20),
          axis.title.y = element_text(color = "black",size = 20),
          legend.title = element_text(color = "black",size = 20),
          legend.text = element_text(colour = "black",size = 20))
dev.off()

# pdf(paste0(file="/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio.pergene.chr.H3K27me3.pdf"),height=5.5,width=5.5)
png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio.pergene.chr.H3K27me3.png"), type="cairo",width = 4096, height = 2048, units = "px",res = 450)
ggplot(omnisummary, aes(x=(((end - start)/2)+start)/1000000, y=`13daf.H3K27me3_ratio10`, color = chr)) +
    geom_point(size=2/3) +
    geom_hline(yintercept = 0.5) +
    geom_vline(data = chromsizes, xintercept = chromsizes[1:8,]$V2/1000000, color = c("#01B4C6","#97ce4c","#FFF874","#BEE5FD","#F675DA","#44281d","#F8D3AC","#E64358")) +
    scale_colour_manual(values=c("chr1"="#01B4C6","chr2"="#97ce4c","chr3"="#FFF874","chr4"="#BEE5FD","chr5"="#F675DA","chr6"="#44281d","chr7"="#F8D3AC","chr8"="#E64358"),name="Chromosome") +
    scale_y_continuous(name = "H3K27me3 maternal ratio",limits = c(0,1)) +
    scale_x_continuous(name = "Position on chromosome (Mb)") +
    guides(colour = guide_legend(override.aes = list(size=5,alpha=1))) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 20),
          axis.text.y = element_text(color = "black",size = 20),
          axis.title.x = element_text(color = "black",size = 20),
          axis.title.y = element_text(color = "black",size = 20),
          legend.title = element_text(color = "black",size = 20),
          legend.text = element_text(colour = "black",size = 20))
dev.off()


################################################################################################################################################
##Scatter plots with enrichment values shown as triangles if too large
##Also show WT and mutant in same plot
################################################################################################################################################
##Enrichment by enrichment with cut off values as triangles
omnisummary %>%
         dplyr::mutate(
             out_of_bounds1 = (asinh(`13daf.H3K27me3_enrich`) > 2 | asinh(`ez23-16daf.H3K27me3_enrich`) > 2),
             out_of_bounds2 = asinh(`13daf.H3K27me3_enrich`) > 2,
             out_of_bounds3 = asinh(`ez23-16daf.H3K27me3_enrich`) > 2,
             `13daf.H3K27me3_enrich` = ifelse(out_of_bounds2, sinh(2), `13daf.H3K27me3_enrich`),
             `ez23-16daf.H3K27me3_enrich` = ifelse(out_of_bounds3, sinh(2), `ez23-16daf.H3K27me3_enrich`)
         ) %>%
    ggplot(aes(x = asinh(`13daf.H3K27me3_enrich`), y = asinh(`ez23-16daf.H3K27me3_enrich`), shape = out_of_bounds1)) +
        geom_point(alpha = 1/5, size = 2,show.legend = F) +
        geom_abline(slope = 1,intercept = c(0,0)) +
        lims(x = c(-0.05, 2.05), y = c(-0.05, 2.05)) +
        labs(x = paste0("arcsinh(WT H3K27me3 enrichment)"), y=paste0("arcsinh(Mutant H3K27me3 enrichment)")) + 
        theme_classic() +
        theme(axis.text.x = element_text(color = "black",size = 36),
              axis.text.y = element_text(color = "black",size = 36),
              axis.title.x = element_text(color = "black",size = 36),
              axis.title.y = element_text(color = "black",size = 36),
              legend.title = element_text(color = "black",size = 36),
              legend.text = element_text(colour = "black",size = 36))

##Enrichment by transcription maternal ratio with cut off values as triangles
omnisummary %>%
         dplyr::mutate(
             out_of_bounds = asinh(`ez23-16daf.H3K27me3_enrich`) > 2,
             `ez23-16daf.H3K27me3_enrich` = ifelse(out_of_bounds, sinh(2), `ez23-16daf.H3K27me3_enrich`)
         ) %>%
    ggplot(aes(x = `ez23-16daf_ratio50`, y = asinh(`ez23-16daf.H3K27me3_enrich`), shape = out_of_bounds)) +
        geom_vline(xintercept = 0.5) +
        geom_point(alpha = 1/5, size = 2,show.legend = F) +
        lims(x = c(-0.05, 1.05), y = c(-0.05, 2.05)) +
        labs(x = paste0("Transcription maternal ratio)"), y=paste0("arcsinh(H3K27me3 enrichment)")) + 
        theme_classic() +
        theme(axis.text.x = element_text(color = "black",size = 36),
              axis.text.y = element_text(color = "black",size = 36),
              axis.title.x = element_text(color = "black",size = 36),
              axis.title.y = element_text(color = "black",size = 36),
              legend.title = element_text(color = "black",size = 36),
              legend.text = element_text(colour = "black",size = 36))

data.frame("H3K27me3_enrichment"=c(omnisummary$`13daf.H3K27me3_enrich`,omnisummary$`ez23-16daf.H3K27me3_enrich`),"Transcription_maternal_ratio"=c(omnisummary$`13daf_ratio50`,omnisummary$`ez23-16daf_ratio50`),"Genotype"=c(rep("Wild type",nrow(omnisummary)),rep("Mutant",nrow(omnisummary)))) %>%
    dplyr::mutate(
        out_of_bounds = asinh(`H3K27me3_enrichment`) > 2,
        `H3K27me3_enrichment` = ifelse(out_of_bounds, sinh(2), `H3K27me3_enrichment`)
    ) %>%
    ggplot(aes(x = `Transcription_maternal_ratio`, y = asinh(`H3K27me3_enrichment`), shape = out_of_bounds,color=Genotype)) +
    geom_vline(xintercept = 0.5) +
    geom_point(alpha = 1/5, size = 2,show.legend = T) +
    guides(shape="none") +
    lims(x = c(-0.05, 1.05), y = c(-0.05, 2.05)) +
    labs(x = paste0("Transcription maternal ratio"), y=paste0("arcsinh(H3K27me3 enrichment)")) +  
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))

##H3K27me3 maternal ratio by transcription maternal ratio with cut off values as triangles
    ggplot(omnisummary,aes(x = `ez23-16daf.H3K27me3_ratio10`, y = `ez23-16daf_ratio50`)) +
        geom_vline(xintercept = 0.5) +
        geom_hline(yintercept = 0.5) +
        geom_point(alpha = 1/5, size = 2,show.legend = F) +
        lims(x = c(0, 1), y = c(0,1)) +
        labs(x = paste0("H3K27me3 maternal ratio"), y=paste0("Transcription maternal ratio)")) + 
        theme_classic() +
        theme(axis.text.x = element_text(color = "black",size = 36),
              axis.text.y = element_text(color = "black",size = 36),
              axis.title.x = element_text(color = "black",size = 36),
              axis.title.y = element_text(color = "black",size = 36),
              legend.title = element_text(color = "black",size = 36),
              legend.text = element_text(colour = "black",size = 36))

data.frame("H3K27me3_maternal_ratio"=c(omnisummary$`13daf.H3K27me3_ratio10`,omnisummary$`ez23-16daf.H3K27me3_ratio10`),"Transcription_maternal_ratio"=c(omnisummary$`13daf_ratio50`,omnisummary$`ez23-16daf_ratio50`),"Genotype"=c(rep("Wild type",nrow(omnisummary)),rep("Mutant",nrow(omnisummary)))) %>%
    ggplot(aes(x = `H3K27me3_maternal_ratio`, y = `Transcription_maternal_ratio`, color=Genotype)) +
        geom_vline(xintercept = 0.5) +
        geom_hline(yintercept = 0.5) +
        geom_point(alpha = 1/5, size = 2,show.legend = T) +
        guides(shape="none") +
        lims(x = c(0, 1), y = c(0,1)) +
        labs(x = paste0("H3K27me3 maternal ratio"), y=paste0("Transcription maternal ratio")) +  
        theme_classic() +
        theme(axis.text.x = element_text(color = "black",size = 36),
              axis.text.y = element_text(color = "black",size = 36),
              axis.title.x = element_text(color = "black",size = 36),
              axis.title.y = element_text(color = "black",size = 36),
              legend.title = element_text(color = "black",size = 36),
              legend.text = element_text(colour = "black",size = 36))

##Enrichment by maternal ratio for WT and mutant with cut off values as triangles
png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/H3K27me3.WT-mutant.enrich-ratio.png"), type="cairo",width = 1024, height = 1024, units = "px")
data.frame("H3K27me3_enrichment"=c(omnisummary$`13daf.H3K27me3_enrich`,omnisummary$`ez23-16daf.H3K27me3_enrich`),"H3K27me3_maternal_ratio"=c(omnisummary$`13daf.H3K27me3_ratio10`,omnisummary$`ez23-16daf.H3K27me3_ratio10`),"Genotype"=c(rep("Wild type",nrow(omnisummary)),rep("Mutant",nrow(omnisummary)))) %>%
    dplyr::mutate(
        out_of_bounds = asinh(`H3K27me3_enrichment`) > 2,
        `H3K27me3_enrichment` = ifelse(out_of_bounds, sinh(2), `H3K27me3_enrichment`)
    ) %>%
    ggplot(aes(x = `H3K27me3_maternal_ratio`, y = asinh(`H3K27me3_enrichment`), shape = out_of_bounds,color=Genotype)) +
    geom_vline(xintercept = 0.5) +
    geom_point(alpha = 1/5, size = 2,show.legend = T) +
    scale_color_manual(values=c("Wild type"="#4682B4","Mutant"="#B47846")) +
    guides(shape="none") +
    lims(x = c(0,1), y = c(0, 2.01)) +
    labs(x = paste0("H3K27me3 maternal ratio"), y=paste0("arcsinh(H3K27me3 enrichment)")) + 
    guides(colour = guide_legend(override.aes = list(size=10,alpha=1))) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()

##PCA plot
library(plotly)
library("DESeq2")


load("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/rsem_everything/output/dds.RData")
vst <- vst(dds, blind = T)
vst_mat <- assay(vst)
row.names(vst_mat) <- gsub("_.*","",row.names(vst_mat))
pca <- prcomp(t(vst_mat))

samples <- read.table("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/rsem_everything/samples.Takv6.Cam2snp.paper.txt",header = T)

df <- cbind(samples,pca$x)

pca <- prcomp(t(vst_mat[,c(1:25)]))
df <- cbind(samples[c(1:25),],pca$x)
pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/rsem_everything/output/pca.1n-2n.pdf"),height=6,width=8)
ggplot(df, aes(x=PC1,y=PC2,color=condition)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(size=5, shape=19) +
    labs(x = "PC1 (62%)", y="PC2 (15%)",color="Sample") +
    theme_minimal() +
    scale_color_manual(values=c("Tak1"='#69C4BA',"Cam2"='#EE75AB',"13daf"='#d9bf6f',"an"="#284A46","ar"="#E2C3D0"),labels=c("Tak1"="Male vegetative","Cam2"="Female vegetative","13daf"="Wild-type embryo","an"="Male sexual organ","ar"="Female sexual organ")) +
    theme(axis.text.x = element_text(color = "black",size = 21),
          axis.text.y = element_text(color = "black",size = 21),
          axis.title.x = element_text(color = "black",size = 21),
          axis.title.y = element_text(color = "black",size = 21),
          legend.title = element_text(color = "black",size = 21),
          legend.text = element_text(colour = "black",size = 21))
dev.off()

###################
##Other plots
######################
##Plot transcript maternal ratio for 1n vs 13daf degs
tmp <- omnisummary
camdeg <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/rsem_everything/output/13daf_vs_Cam2.csv",header = T,row.names = 2,nrows = 18274)
takdeg <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/rsem_everything/output/13daf_vs_Tak1.csv",header = T,row.names = 2,nrows = 18274)
tmp$category <- c("Unchanged")
tmp$category[row.names(omnisummary) %in% intersect(row.names(subset(camdeg,camdeg$padj<0.01 & (camdeg$log2FoldChange>=1))),row.names(subset(takdeg,takdeg$padj<0.01 & (takdeg$log2FoldChange>=1))))] <- c("Up")
tmp$category[row.names(omnisummary) %in% intersect(row.names(subset(camdeg,camdeg$padj<0.01 & (camdeg$log2FoldChange <= -1))),row.names(subset(takdeg,takdeg$padj<0.01 & (takdeg$log2FoldChange <= -1))))] <- c("Down")
png(filename =paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/plots/matratio/1nvs13daf.deg.TPM-matratio.png"), type="cairo",width = 1024, height = 1024, units = "px")            
ggplot(tmp, aes(x=category, y=`13daf_ratio50`,color=category)) + 
    geom_hline(yintercept=0.5) + 
    geom_violin(position = position_dodge(0.9)) +
    geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA) +
    geom_point(alpha = 1/20, size = 0.5,position = position_jitterdodge(jitter.width = 0.8,dodge.width = 0.9)) +
    #stat_compare_means(comparisons = list(c("Up","Unchanged"),c("Up","Down"),c("Down","Unchanged")), label="p.format",tip.length = 0.01, label.y = c(1,1.05,1.1)) +
    scale_x_discrete(name=c("Differential expression")) +
    scale_y_continuous(name=paste0("Transcription maternal ratio"),limits = c(0,1.1),breaks=c(0,0.25,0.5,0.75,1)) +
    labs(x = paste0("DEG")) +
    theme_classic() +
    scale_color_manual(values=c("Unchanged"="#949494","Up"="#4682B4","Down"="#B47846"),guide = FALSE) +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()
pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/rna-seq/Takv6.Cam2snp/OmniAnalysis/plots/matratio/1nvs13daf.deg.TPM-matratio.pdf"),height=12,width=12)
ggplot(tmp, aes(x=category, y=`13daf_ratio50`,color=category)) + 
    geom_hline(yintercept=0.5) + 
    geom_violin(position = position_dodge(0.9),size=1) +
    geom_boxplot(position=position_dodge(0.9),width=0.1,outlier.color = NA,size=1) +
    scale_x_discrete(name=c("Differential expression")) +
    scale_y_continuous(name=paste0("Transcription maternal ratio"),limits = c(0,1.1),breaks=c(0,0.25,0.5,0.75,1)) +
    labs(x = paste0("DEG")) +
    theme_classic() +
    scale_color_manual(values=c("Unchanged"="#949494","Up"="#4682B4","Down"="#B47846"),guide = FALSE) +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()

##Plot maternal ratio per chromosome
tmp <- omnisummary
tmp$chr <- gsub("Mp","chr",gsub("g.*","",row.names(tmp)))
pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/cutrun/Takv6.Cam2snp/OmniAnalysis/plots/matratio/perchr.matratio.H3K27me3.pdf"),height=2.35,width=2.35)
ggplot(melt(tmp[grep("[0-9]",tmp$chr),],id.vars = "chr",measure.vars = "13daf.H3K27me3_ratio10"),aes(x=chr,y=value)) +
    geom_hline(yintercept = 0.5,size=0.25) +
    geom_violin(size=0.25) +
    geom_boxplot(width=0.1,outlier.color = NA,size=0.25) +
    lims( y = c(0,1)) +
    labs(x = "Chromosome", y="H3K27me3 maternal ratio") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",size = 8, hjust=1,angle=45),
          axis.text.y = element_text(color = "black",size = 8),
          axis.title.x = element_text(color = "black",size = 8),
          axis.title.y = element_text(color = "black",size = 8),
          legend.title = element_text(color = "black",size = 8),
          legend.text = element_text(colour = "black",size = 8),
          text=element_text(size=6,  family="Helvetica"),
          axis.line=element_line(size=0.25),
          axis.ticks = element_line (size = 0.25),
          plot.margin = unit(c(0,0,0,0), "cm")) 
dev.off()

##Show sex chromosome gene expression
tpmsex <- read.csv("/groups/berger/user/sean.montgomery/Documents/rna-seq/Tak1v5.chrU/OmniAnalysis/TPM.csv", header = T, row.names = 1, check.names=FALSE)
c2 <- asinh(tpmsex[,c("Tak1","Cam2","13 DAF","16 DAF","21 DAF","25 DAF","spore")])
c3 <- c2[grep("Mp[UV]",row.names(c2)),]
c3$chr <- gsub("g.*","",row.names(c3))
c4 <- subset(c3, !(c3$Tak1>asinh(0) & c3$chr=="MpU") & !(c3$Cam2>asinh(0) & c3$chr=="MpV") & rowSums(c3[,-(ncol(c3))])>0 & !(c3$Tak1==0 & c3$chr=="MpV") & !(c3$Cam2==0 & c3$chr=="MpU"))
c5 <- melt(c4,id.vars = c("chr"))
c5$grouping <- paste0(c5$chr,"_",c5$variable)
ggplot(c5,aes(x=variable,y=value,color=chr)) +
    geom_violin(position=position_dodge(1)) +
    geom_boxplot(position=position_dodge(1),width=0.1,outlier.color = NA) +
    geom_point(alpha = 1/20, size = 0.5,position = position_jitterdodge(jitter.width = 0.3,dodge.width = 1)) +
    stat_compare_means(aes(group = chr), label = "p.format") +
    theme_classic() +
    scale_color_manual(values=c("MpU"="#4682B4","MpV"="#B47846")) +
    theme(axis.text.x = element_text(color = "black",size = 18),
          axis.text.y = element_text(color = "black",size = 18),
          axis.title.x = element_text(color = "black",size = 18),
          axis.title.y = element_text(color = "black",size = 18),
          legend.title = element_text(color = "black",size = 18),
          legend.text = element_text(colour = "black",size = 18)) 

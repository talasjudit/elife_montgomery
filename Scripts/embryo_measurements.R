library(ggplot2)
##Measurement of embryo sizes between WT and mutant
measurements <- read.table("/groups/berger/user/sean.montgomery/Documents/rna-seq/Dissection_pics/Sporophyte_measurements.txt",header = T)
measurements$Age <- factor(measurements$Age)
measurements$Stage <- factor(measurements$Stage)
ezmeasurements <- read.table("/groups/berger/user/sean.montgomery/Documents/rna-seq/Dissection_pics/ez2ez3/ez-Sporophyte_measurements.txt",header = T)
ezmeasurements$Age <- factor(ezmeasurements$Age)
ezmeasurements$Stage <- factor(ezmeasurements$Stage)
ezmeasurements$Cross <- c("Mutant")
measurements$Cross <- c("Wild type")
combomeasure <- rbind(measurements,ezmeasurements)
meltmeasure <- melt(combomeasure,id.vars = c("Cross","Age"),measure.vars = c("Area"))
meltmeasure$Cross <- factor(meltmeasure$Cross, levels = c("Wild type","Mutant"))
pdf("/groups/berger/user/sean.montgomery/Documents/rna-seq/Dissection_pics/Size-comparison.pdf", width = 12, height = 8)
ggboxplot(meltmeasure[c(grep("13",meltmeasure$Age),grep("16",meltmeasure$Age)),], x = "Age", y = "value", color="Cross" ,add = "jitter") +
    theme_classic() +
    scale_x_discrete(name="Age",labels=c("13"="13daf","16"="16daf"))+
    scale_y_continuous(limits=c(0,0.6),name=bquote('Embryo area'~(mm^2))) +
    scale_color_manual(values=c("Wild type"="#4682B4","Mutant"="#B47846")) +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)) +
    stat_compare_means(aes(group = Cross),label = "p.signif",size=8)
dev.off()

##Boxplot with stats comparing fertility of mutant to WT 
library("reshape2")
library(plyr)
library(ggpubr)
fertility <- read.csv("/groups/berger/user/sean.montgomery/Documents/Fertility/Fertility_check.csv",header = T)
fertility$alive <- fertility$X..sporophytes.alive/(fertility$X..with.perianths+0.01)*100
fertility$spore <- fertility$X..yellow.sporophytes/(fertility$X..with.perianths+0.01)*100
fertility$Cross <- gsub("Cam-2 ez2 ez3 #19 x Tak-1","Mutant",gsub("Cam-2 x Tak-1","Wild type",fertility$Cross))
fertmelt <- melt(fertility, id.vars = c("Age","Cross"), measure.vars = c("alive","spore"))
fertmelt$Cross <- factor(fertmelt$Cross, levels = c("Wild type","Mutant"))

pdf("/groups/berger/user/sean.montgomery/Documents/rna-seq/Dissection_pics/fertility.alive.pdf", width = 12, height = 8)
ggboxplot(fertmelt[grep("alive",fertmelt$variable),], x = "Age", y = "value", color="Cross" ,add = "jitter") +
    theme_classic() +
    scale_x_discrete(name="Age")+
    scale_y_continuous(limits = c(0,100),name="Surviving embryos (%)") +
    scale_color_manual(values=c("Wild type"="#4682B4","Mutant"="#B47846")) +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)) +
    stat_compare_means(aes(group = Cross),label = "p.signif",size=8)
dev.off()

pdf("/groups/berger/user/sean.montgomery/Documents/rna-seq/Dissection_pics/fertility.spore.pdf", width = 12, height = 8)
ggboxplot(fertmelt[grep("spore",fertmelt$variable),], x = "Age", y = "value", color="Cross" ,add = "jitter") +
    theme_classic() +
    scale_x_discrete(name="Age")+
    scale_y_continuous(limits = c(0,100),name="Embryos producing spores (%)") +
    scale_color_manual(values=c("Wild type"="#4682B4","Mutant"="#B47846")) +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)) +
    stat_compare_means(aes(group = Cross),label = "p.signif",size=8)
dev.off()



################################################################################################################################################
##Plots of 1n mutant vs WT
################################################################################################################################################
##Measure gemmae area
tmp <- read.table("/groups/berger/user/sean.montgomery/Documents/Fertility/gemmae_areas.xls", header = T)
tmp$Age <- as.numeric(gsub("d","",tmp$Age))
pdf("/groups/berger/user/sean.montgomery/Documents/Fertility/Gemmae_area.pdf", width = 16, height = 12)
ggplot(subset(tmp,Genotype == "Cam2" | Genotype == "Cam2-ez23-19"), aes(x=Age, y=Area, color=Genotype)) +
    theme_classic() +
    scale_x_continuous(name="Age (days)",limits=c(0,14),breaks=c(0,7,14)) +
    scale_y_continuous(name=bquote('Gemma area'~(mm^2))) +
    scale_color_manual(labels=c("Cam2"="Cam-2","Cam2-ez23-19"=expression("Cam-2"~italic("e(z)2/e(z)3"~"#19"))),values=c("Cam2"="#4682B4","Cam2-ez23-19"="#B47846")) +
    geom_smooth() +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()

##Measure first gametangiophore emergence
tmp <- read.csv("/groups/berger/user/sean.montgomery/Documents/Fertility/Gametangiophore_emerge.csv",header = T)
pdf("/groups/berger/user/sean.montgomery/Documents/Fertility/Gametangiophore_emerge.pdf", width = 12, height = 12)
ggboxplot(tmp[grep("Cam",tmp$Line),], x = "Line", y = "Days",add = "jitter") +
    theme_classic() +
    scale_x_discrete(name="Line", label = c("Cam-2",expression("Cam-2"~italic("e(z)2/e(z)3"~"#19")))) +
    scale_y_continuous(limits = c(0,35),name="First sexual organ (days)") +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36)) +
    stat_compare_means(comparisons = list(c("Cam-2","Cam-2 ez2 ez3 #19")),label = "p.signif",size=8)
dev.off()
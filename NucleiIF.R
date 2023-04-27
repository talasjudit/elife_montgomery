##Analysis of immunostaining data from masked images
immuno <- read.csv("/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Immuno_result_MS.csv")

##Number of heterochromatic foci
pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/nCyanMask.pdf"),height=12,width=12)
ggplot(subset(immuno, Path.1 == "Tak2xTak1"), aes(x=Path.2, y=nCyanBlobs)) +
    geom_boxplot(outlier.color = NA) +
    geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.5) +
    theme_minimal() +
    scale_x_discrete(name="Immunostain",labels=c("H3K27me3","H3K36me3","H3K9me1")) +
    scale_y_continuous(name="Number of heterochromatic foci") +
    theme(axis.text.x = element_text(color = "black",size = 36, angle = 45,hjust = 1),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()

##Percentage of nucleus area covered by heterochromatic foci
pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/cyanAreaPerc.pdf"),height=12,width=12)
ggplot(subset(immuno, Path.1 == "Tak2xTak1"), aes(x=Path.2, y=cyanAreaPerc*100)) +
    geom_boxplot(outlier.color = NA) +
    geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.5) +
    theme_minimal() +
    scale_x_discrete(name="Immunostain",labels=c("H3K27me3","H3K36me3","H3K9me1")) +
    scale_y_continuous(name="Percent of nucleus covered by foci", limits=c(0,25)) +
    theme(axis.text.x = element_text(color = "black",size = 36, angle = 45,hjust = 1),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()

##Percentage of IF signal within heterochromatic foci
pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/greenSignalPerc.pdf"),height=12,width=12)
ggplot(subset(immuno, Path.1 == "Tak2xTak1"), aes(x=Path.2, y=greenSignalPerc*100)) +
    geom_boxplot(outlier.color = NA) +
    geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.5) +
    stat_compare_means(comparisons = list(c("H3K27me3_488","H3K36me3_488"),c("H3K27me3_488","H3K9me1_488"),c("H3K9me1_488","H3K36me3_488")), label = "p.signif",tip.length = 0.01, size=8) +
    theme_minimal() +
    scale_x_discrete(name="Immunostain",labels=c("H3K27me3","H3K36me3","H3K9me1")) +
    scale_y_continuous(name="Percent of immunostain signal within foci",limits=c(0,100),breaks = c(0,25,50,75,100)) +
    theme(axis.text.x = element_text(color = "black",size = 36, angle = 45,hjust = 1),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()

##Percentage of overlap between heterochromatic foci and H3K27me3 foci
pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/OverlapGreen.pdf"),height=12,width=12)
ggplot(subset(immuno, Path.1 == "Tak2xTak1" & Path.2 == "H3K27me3_488"), aes(x=Path.2, y=Overlap.Green*100)) +
    geom_boxplot(outlier.color = NA) +
    geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.5) +
    theme_minimal() +
    scale_x_discrete(name="Immunostain",labels=c("H3K27me3")) +
    scale_y_continuous(name="Percent of foci area within immunostain foci",limits=c(0,100),breaks = c(0,25,50,75,100)) +
    theme(axis.text.x = element_text(color = "black",size = 36),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_text(color = "black",size = 36),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()
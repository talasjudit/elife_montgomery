##Read in FISH-IF measurements of the percentage of FISH signal overlap with immuno signal
tmp <- read.table("/groups/berger/user/sean.montgomery/Documents/FISH/16_7_21/U_H3K9/Results.H3K9me1.chrU.xls",header = T)
tmp$sample <- c("H3K9me1.chrU")
fish <- tmp
tmp <- read.table("/groups/berger/user/sean.montgomery/Documents/FISH/16_7_21/V_H3K9/Results.H3K9me1.chrV.xls",header = T)
tmp$sample <- c("H3K9me1.chrV")
fish <- rbind(fish,tmp)
tmp <- read.table("/groups/berger/user/sean.montgomery/Documents/FISH/8_5_21/H3K27me3_U/Results.H3K27me3.chrU.xls",header = T)
tmp$sample <- c("H3K27me3.chrU")
fish <- rbind(fish,tmp)
tmp <- read.table("/groups/berger/user/sean.montgomery/Documents/FISH/8_5_21/H3K27_me3_V/Results.H3K27me3.chrV.xls",header = T)
tmp$sample <- c("H3K27me3.chrV")
fish <- rbind(fish,tmp)

pdf(file=paste0("/groups/berger/user/sean.montgomery/Documents/FISH/immunoFISHquant.pdf"),height=12,width=12)
ggplot(fish,aes(x=sample,y=X.Area)) +
    geom_boxplot() +
    geom_dotplot(binaxis = "y", stackdir = "center",dotsize=0.5,binwidth = 3) +
    stat_compare_means(comparisons = list(c("H3K27me3.chrU","H3K27me3.chrV"),c("H3K9me1.chrU","H3K9me1.chrV"),c("H3K9me1.chrU","H3K27me3.chrU"),c("H3K9me1.chrV","H3K27me3.chrV")), label = "p.signif",tip.length = 0.01, label.y = c(90,95,100,105),size=8) +
    theme_minimal() +
    scale_x_discrete(name="Sample") +
    scale_y_continuous(name="FISH within immunostain signal (%)",limits = c(0,110),breaks=c(0,25,50,75,100)) +
    theme(axis.text.x = element_text(color = "black",size = 36, angle = 45,hjust = 1),
          axis.text.y = element_text(color = "black",size = 36),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black",size = 36),
          legend.title = element_text(color = "black",size = 36),
          legend.text = element_text(colour = "black",size = 36))
dev.off()



library(ggplot2)

x <- read.table("enhancer.enrichment.txt")

c1 <- c(mean(x$V1)/mean(x$V2),mean(x$V3)/mean(x$V4))

c2 <- data.frame(type=c("AtoB","BtoA"),value=c1)
pdf("Enrich.enhancer.pdf",width=3,height=4)
ggplot(c2,aes(x=type,y=log2(value))) + geom_bar(stat="identity",width=0.5) + xlab("") + ylab("log2 Enrichment") + theme_classic()
dev.off()

###############################################
up <- read.table("Upregulated-genes-RNA-seq.txt",header=T)

btoa <- read.table("mm9_TSS_GeneName_BtoA.bed")

library(VennDiagram)
T1<-venn.diagram(x =list(Up = (up$GeneName),BtoA = btoa$V5), filename = NULL, height = 600, width= 600, resolution =150, imagetype="png", fill =c("cornflowerblue","green"),alpha=0.5,cex=2.5,cat.cex=2.5,cat.pos=c(350, 150))
grid.draw(T1)


###############################################
atob_neg <- read.table("repeat_expression50/Neg-AtoB.rpkm")
atob_pos <- read.table("repeat_expression50/Pos-AtoB.rpkm")
btoa_neg <- read.table("repeat_expression50/Neg-BtoA.rpkm")
btoa_pos <- read.table("repeat_expression50/Pos-BtoA.rpkm")

boxplot(btoa_neg$V4,btoa_pos$V4,atob_neg$V4,atob_pos$V4,outline = F,col = c("blue","red","blue","red"))

t.test(atob_neg$V4,atob_pos$V4)
t.test(btoa_neg$V4,btoa_pos$V4)

wilcox.test(atob_neg$V4,atob_pos$V4,paired = T)
wilcox.test(btoa_neg$V4,btoa_pos$V4,paired = T)


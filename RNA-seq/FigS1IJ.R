library(stringr)
up_pathway <- read.table("Up-genes-pathway-1.5.txt",sep="\t",header=T)

up_pathway <- up_pathway[4:1,]
library(ggplot2)
up_pathway$Ingenuity.Canonical.Pathways <- factor(up_pathway$Ingenuity.Canonical.Pathways,levels=up_pathway$Ingenuity.Canonical.Pathways)

up_pathway$pvalue = 10 ^ (-1 * up_pathway$X.log.p.value.)
up_pathway$fdr <- p.adjust(up_pathway$pvalue, method="fdr")
up_pathway$log.fdr <- -log10(up_pathway$fdr)

p <- ggplot(up_pathway,aes(x=Ingenuity.Canonical.Pathways,y=(log.fdr))) + geom_bar(stat = "identity",fill="red4") + coord_flip() + xlab("") + ylab("-log10 FDR")

p + scale_x_discrete(labels=function(x) str_wrap(x,width=30))
ggsave("Up-GO-FDR.pdf",width = 6,height = 4)



up_pathway <- read.table("Down-genes-pathway-1.5.txt",sep="\t",header=T)

up_pathway <- up_pathway[4:1,]
library(ggplot2)
up_pathway$Ingenuity.Canonical.Pathways <- factor(up_pathway$Ingenuity.Canonical.Pathways,levels=up_pathway$Ingenuity.Canonical.Pathways)

up_pathway$pvalue = 10 ^ (-1 * up_pathway$X.log.p.value.)
up_pathway$fdr <- p.adjust(up_pathway$pvalue, method="fdr")
up_pathway$log.fdr <- -log10(up_pathway$fdr)

p <- ggplot(up_pathway,aes(x=Ingenuity.Canonical.Pathways,y=(log.fdr))) + geom_bar(stat = "identity",fill="blue4") + coord_flip() + xlab("") + ylab("-log10 FDR")

p + scale_x_discrete(labels=function(x) str_wrap(x,width=30))
ggsave("Down-GO-FDR.pdf",width = 6,height = 4)


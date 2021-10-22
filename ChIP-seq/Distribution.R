x <- read.csv("Distribution.csv",header = F)
x1 <- data.frame(t(x))
s2 <- as.vector(t(x1[1,]))
colnames(x1) <- s2
x1 <- x1[-1,]

library(ggplot2)
library(reshape2)

x2 <- melt(x1,c("Mod", "Genomic"))
x2$s <- paste(x2$Mod,x2$Genomic)
x2$value <- as.numeric(x2$value)
x2$s <- factor(x2$s,levels = c("H3K27ac ESC","H3K27ac TC","H3K4me3 ESC","H3K4me3 TC","H3K4me1 ESC","H3K4me1 TC","H3K27me3 ESC","H3K27me3 TC","H3K9me3 ESC","H3K9me3 TC"))

#x2$variable <- factor(x2$variable,levels = rev(c("TSS","Intergenic","TES","Genebody")))
ggplot(x2,aes(x=s,y=value,fill=variable)) + geom_bar(stat='identity',width=0.5) + theme_classic()

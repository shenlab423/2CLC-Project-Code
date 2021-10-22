x1 <- read.table("SuperEnhancer-Nearest-genes.txt")
tss <- read.table("mm9.TSS.bed")

x2 <- merge(tss,x1,by.x="V4",by.y="V1")
x3 <- data.frame(GeneName=x2$V5)
x3 <- unique(x3)

allgenes <- read.table("Allgenes.txt")
x4 <- allgenes[rownames(allgenes) %in% x3$GeneName,]

boxplot(log2(x4$Neg_mean),log2(x4$Pos_mean),outline = F,ylim=c(0,11),names=c("ESC","2CLC"),col = c("blue","red"),ylab="log2 TMM",main="Super-enhancer neighboring genes")
wilcox.test(log2(x4$Neg_mean),log2(x4$Pos_mean),paired = T)


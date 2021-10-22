### Oct4 R 3.4.1
### paste RNA-2C-CTCF-1.unique.count.txt RNA-2C-CTCF-2.unique.count.txt RNA-2C-N-1.unique.count.txt RNA-2C-N-2.unique.count.txt RNA-2C-P-1.unique.count.txt RNA-2C-P-2.unique.count.txt > Shenlab.2C.exon.fix.count

library(edgeR)
library(dplyr)

table <- read.table("RNA-IAA-2d-2C.count.txt")
table <- table[,c(1,2,4,6,8)]
colnames(table) <- c("GeneSymbol","N1","N2","P1","P2")

rownames(table) <- table$GeneSymbol
table <- table[,-1]
library(edgeR)

cpm_log <- cpm(table, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
hist(median_log2_cpm)
expr_cutoff <- -1
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(median_log2_cpm > expr_cutoff)
table2 <- table[median_log2_cpm > expr_cutoff, ]

group <- c(1,1,2,2)
y <- DGEList(counts = table2, group = group)
y <- calcNormFactors(y,method = "TMM")


counts.per.m <- cpm(y, normalized.lib.sizes=TRUE)
boxplot(counts.per.m,outline = F)
counts.per.m <- data.frame(counts.per.m)
counts.per.m$GeneSymbol <- rownames(table2)

counts.per.m <- counts.per.m[-(grep("^Mir",counts.per.m[,"GeneSymbol"])),]
counts.per.m <- counts.per.m[-(grep("^Sno",counts.per.m[,"GeneSymbol"])),]
counts.per.m <- counts.per.m[-(grep("^Hist",counts.per.m[,"GeneSymbol"])),]


rownames(counts.per.m) <- counts.per.m$GeneSymbol
counts.per.m <- counts.per.m[,-5]
counts.per.m[counts.per.m<1]=1
counts.per.m$N <- apply(counts.per.m[,1:2], 1, mean)
counts.per.m$P <- apply(counts.per.m[,3:4], 1, mean)


counts.per.m$FC <- log2(counts.per.m$P) - log2(counts.per.m$N)

write.table(data.frame(counts.per.m),"Allgenes_CTCF.txt",col.names = T,row.names = T,sep="\t",quote = F)


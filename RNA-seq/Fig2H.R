### Oct4 R 3.4.1

library(edgeR)
library(dplyr)

table <- read.table("Shenlab.2C.exon.fix.2.count")
colnames(table) <- c("GeneSymbol","2C_Neg1_1","2C_Neg1_2","2C_Neg2_1","2C_Neg2_2","2C_Pos1_1","2C_Pos1_2","2C_Pos2_1","2C_Pos2_2")

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

group <- c(1,1,1,1,2,2,2,2)
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
counts.per.m <- counts.per.m[,-9]
counts.per.m[counts.per.m<1]=1
counts.per.m$Neg_mean <- apply(counts.per.m[,1:4],1,mean)
counts.per.m$Pos_mean <- apply(counts.per.m[,5:8], 1, mean)
counts.per.m$FC <- log2(counts.per.m$Pos_mean) - log2(counts.per.m$Neg_mean )

write.table(data.frame(counts.per.m),"Allgenes.txt",col.names = T,row.names = T,sep="\t",quote = F)



counts.per.m$Genename <- (rownames(counts.per.m))
plurinet <- read.table("plurinet.txt")

x2 <- merge(plurinet,counts.per.m,by.x="V1",by.y="Genename",all.x = T)

library(ggplot2)
pdf("Fig2H.pdf",width = 5,height = 5)
g <- ggplot(x2, aes(x=log2(Neg_mean), y=log2(Pos_mean),colour="type")) + geom_point() + scale_color_manual(values=c("blue")) + theme(legend.position = "none") + scale_y_continuous(limits = c(0,12)) + scale_x_continuous(limits = c(0,12))
g <- g + geom_abline(intercept=log2(1),slope=1,colour="#990000", linetype="dashed") + ylab("log2 TMM (2CLC)") + xlab("log2 TMM (ESC)")
g <- g + geom_text(aes(label=V1), size=4, vjust = 0)
g
dev.off()

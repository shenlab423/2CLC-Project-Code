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

group <- c("A","A","A","A","B","B","B","B")
y <- DGEList(counts = table2, group = group)
y <- calcNormFactors(y,method = "TMM")

y <- estimateCommonDisp(y)
y$common.dispersion
plotBCV(y)
s <- exactTest(y,dispersion = "common")
out <- s$table
out$GeneName <- rownames(out)

counts.per.m <- read.table("Allgenes.txt")
counts.per.m$GeneSymbol <- rownames(counts.per.m)

counts.per.m2 <- merge(counts.per.m,out,by="GeneName")

counts.per.m2$FDR <- p.adjust(counts.per.m2$PValue)
up <- counts.per.m2[counts.per.m2$FC>1.5 & counts.per.m2$FDR<0.05,]
down <- counts.per.m2[counts.per.m2$FC < -1.5 & counts.per.m2$FDR<0.05,]
nochange <- counts.per.m2[!counts.per.m2$GeneName %in% up$GeneName & !counts.per.m2$GeneName %in% down$GeneName,]


up$type <- "Up"
down$type <- "Down"
allplot <- rbind(up,down)
library(ggplot2)
allplot$type <- factor(allplot$type,levels = c("Up","Down"))
ggplot(allplot,aes(x=FC,y=-log10(FDR),col=type)) + geom_point() + coord_cartesian(xlim = c(-5,12),ylim = c(0,250))  + theme_classic() + theme(legend.position = "none") + 
  scale_color_manual(values=c("red","blue"))


nochange$type <- "Nochange"
allplot <- nochange
library(ggplot2)
ggplot(allplot,aes(x=FC,y=-log10(FDR),col=type)) + geom_point() + coord_cartesian(xlim = c(-5,12),ylim = c(0,250))  + theme_classic() + theme(legend.position = "none") + 
  scale_color_manual(values=c("grey"))


allplot <- rbind(up,down,nochange)
highlight <- c("Dub1",
               "Zscan4f",
               "Zscan4c",
               "Zscan4d",
               "Tdpoz3",
               "Tdpoz4",
               "Zfp352",
               "Fgf4",
               "Sox2")

table1_in <- allplot[(allplot$GeneName) %in% highlight,]
table1_in$PNAME <- (table1_in$GeneName)

allplot2 <- (table1_in)

library(ggplot2)
g <- ggplot(allplot2, aes(x=FC, y=-log10(FDR),colour=type)) + geom_point() + scale_color_manual(values=c("black","black")) + theme(legend.position = "none") + scale_x_continuous(limits = c(-5,12)) + scale_y_continuous(limits = c(0,250))
g <- g + geom_text(aes(label=PNAME), size=4, vjust = 0)  #,nudge_x = c(0,1,0,0,0,0,0,-1.5,1,0,1.5),nudge_y = c(-0.5,0,0.5,1,0.5,-0.5,0,0,0,0.5,0)
g <- g + annotate(geom="text", x=5, y=200, label="n=1044",color="red")
g <- g + annotate(geom="text", x=-3, y=200, label="n=212",color="blue")
g + theme_classic()+ theme(legend.position = "none") + geom_vline(xintercept = -1.5) + geom_vline(xintercept = 1.5) + geom_hline(yintercept = -log10(0.05))



write.table(data.frame((up$GeneName)),"Upgenes_1.5.txt",col.names = F,row.names = F,sep="\t",quote = F)
write.table(data.frame((down$GeneName)),"Downgenes_1.5.txt",col.names = F,row.names = F,sep="\t",quote = F)
write.table(nochange,"Nochangegenes.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(counts.per.m2,"Allgenes.pvalue.txt",col.names = T,row.names = T,sep="\t",quote = F)


######################

up <- read.table("Upgenes_1.5.txt")
down <- read.table("Downgenes_1.5.txt")
genes <- read.table("GSE71434_FPKM_stage.txt",header=T)

genes_up <- genes[genes$Gene %in% up$V1,]
genes_down <- genes[genes$Gene %in% down$V1,]

par(las=2)
genes_up <- genes_up[,c(4,5,6,7,8,9,10,11)]
genes_up <- genes_up[,c(5,8)]
pdf("2cell.up.pdf",width = 5,height = 5)
boxplot(log2(genes_up+1),outline = F)
dev.off()

genes_down <- genes_down[,c(4,5,6,7,8,9,10,11)]
genes_down <- genes_down[,c(5,8)]
pdf("2cell.down.pdf",width = 5,height = 5)
boxplot(log2(genes_down+1),outline = F)
dev.off()

######################


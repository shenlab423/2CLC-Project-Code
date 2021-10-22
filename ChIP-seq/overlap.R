esc <- read.table("H3K4me3_ESC_peaks.bed")
tc <- read.table("H3K4me3_2C_peaks.bed")

esc$gene <- paste(esc$V1,esc$V2,esc$V3,sep="-")
tc$gene <- paste(tc$V1,tc$V2,tc$V3,sep="-")

library(VennDiagram)

#pdf("OverlapDown.pdf",width = 5,height = 5)
T1<-venn.diagram(x =list(ESC = esc$gene, TCLC= tc$gene), filename = NULL, height = 600, width= 600, resolution =150, imagetype="png", fill =c("cornflowerblue","red"),alpha=0.5,cex=2.5,cat.cex=2.5,cat.pos=c(0,0))
grid.draw(T1)


esc <- read.table("H3K9me3_ESC_peaks.bed")
tc <- read.table("H3K9me3_2C_peaks.bed")

esc$gene <- paste(esc$V1,esc$V2,esc$V3,sep="-")
tc$gene <- paste(tc$V1,tc$V2,tc$V3,sep="-")

library(VennDiagram)

#pdf("OverlapDown.pdf",width = 5,height = 5)
T1<-venn.diagram(x =list(ESC = esc$gene, TCLC= tc$gene), filename = NULL, height = 600, width= 600, resolution =150, imagetype="png", fill =c("cornflowerblue","red"),alpha=0.5,cex=2.5,cat.cex=2.5,cat.pos=c(0,0))
grid.draw(T1)



esc <- read.table("H3K27me3_ESC_peaks.bed")
tc <- read.table("H3K27me3_2C_peaks.bed")

esc$gene <- paste(esc$V1,esc$V2,esc$V3,sep="-")
tc$gene <- paste(tc$V1,tc$V2,tc$V3,sep="-")

library(VennDiagram)

#pdf("OverlapDown.pdf",width = 5,height = 5)
T1<-venn.diagram(x =list(ESC = esc$gene, TCLC= tc$gene), filename = NULL, height = 600, width= 600, resolution =150, imagetype="png", fill =c("cornflowerblue","red"),alpha=0.5,cex=2.5,cat.cex=2.5,cat.pos=c(0,0))
grid.draw(T1)



esc <- read.table("H3K27ac_ESC_peaks.bed")
tc <- read.table("H3K27ac_2C_peaks.bed")

esc$gene <- paste(esc$V1,esc$V2,esc$V3,sep="-")
tc$gene <- paste(tc$V1,tc$V2,tc$V3,sep="-")

library(VennDiagram)

#pdf("OverlapDown.pdf",width = 5,height = 5)
T1<-venn.diagram(x =list(ESC = esc$gene, TCLC= tc$gene), filename = NULL, height = 600, width= 600, resolution =150, imagetype="png", fill =c("cornflowerblue","red"),alpha=0.5,cex=2.5,cat.cex=2.5,cat.pos=c(0,0))
grid.draw(T1)



esc <- read.table("H3K4me1_ESC_peaks.bed")
tc <- read.table("H3K4me1_2C_peaks.bed")

esc$gene <- paste(esc$V1,esc$V2,esc$V3,sep="-")
tc$gene <- paste(tc$V1,tc$V2,tc$V3,sep="-")

library(VennDiagram)

#pdf("OverlapDown.pdf",width = 5,height = 5)
T1<-venn.diagram(x =list(ESC = esc$gene, TCLC= tc$gene), filename = NULL, height = 600, width= 600, resolution =150, imagetype="png", fill =c("cornflowerblue","red"),alpha=0.5,cex=2.5,cat.cex=2.5,cat.pos=c(0,0))
grid.draw(T1)

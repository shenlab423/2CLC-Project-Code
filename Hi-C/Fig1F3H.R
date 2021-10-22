tclc_enhancer <- read.table("ESC.Fig1F.txt",sep=" ")
minvalue <- min(tclc_enhancer,na.rm = T)
tclc_enhancer[is.na(tclc_enhancer)] <- 0
tclc_enhancer <- tclc_enhancer + minvalue
library(gplots)
rampCol1 <- colorRampPalette(c("#00004c","#0000ff","#ffffff","#ff0000","#7f0000"))(30)

pdf("ESC.1.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1.5,1.5,0.1))
dev.off()

tclc_enhancer <- read.table("2CLC.Fig1F.txt",sep=" ")
minvalue <- min(tclc_enhancer,na.rm = T)
tclc_enhancer[is.na(tclc_enhancer)] <- 0
tclc_enhancer <- tclc_enhancer + minvalue
library(gplots)
rampCol1 <- colorRampPalette(c("#00004c","#0000ff","#ffffff","#ff0000","#7f0000"))(30)

pdf("2CLC.1.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1.5,1.5,0.1))
dev.off()


tclc_enhancer <- read.table("ESC.Fig3H.txt",sep=" ")
minvalue <- min(tclc_enhancer,na.rm = T)
tclc_enhancer[is.na(tclc_enhancer)] <- 0
tclc_enhancer <- tclc_enhancer + minvalue
library(gplots)
rampCol1 <- colorRampPalette(c("#00004c","#0000ff","#ffffff","#ff0000","#7f0000"))(30)

pdf("ESC.3.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1.5,1.5,0.1))
dev.off()

tclc_enhancer <- read.table("2CLC.Fig3H.txt",sep=" ")
minvalue <- min(tclc_enhancer,na.rm = T)
tclc_enhancer[is.na(tclc_enhancer)] <- 0
tclc_enhancer <- tclc_enhancer + minvalue
library(gplots)
rampCol1 <- colorRampPalette(c("#00004c","#0000ff","#ffffff","#ff0000","#7f0000"))(30)

pdf("2CLC.3.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1.5,1.5,0.1))
dev.off()



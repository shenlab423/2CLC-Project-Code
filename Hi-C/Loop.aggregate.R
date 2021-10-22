
tclc_enhancer <- read.table("ESC.Loop.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("ESC.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()


tclc_enhancer <- read.table("2CLC.Loop.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("2CLC.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()


tclc_enhancer <- read.table("Untreat.Loop.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("Untreat.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()


tclc_enhancer <- read.table("auxin-2days.Loop.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("Auxin.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()

tclc_enhancer <- read.table("washoff-2days.Loop.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("Washoff.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()




tclc_enhancer <- read.table("early-2Cell.Loop.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("early2Cell.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()

tclc_enhancer <- read.table("ICM.Loop.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("ICM.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()


tclc_enhancer <- read.table("ESC.EPI.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("ESC.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()


tclc_enhancer <- read.table("2CLC.EPI.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("2CLC.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()




tclc_enhancer <- read.table("early-2Cell.EPI.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("early2Cell.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()


tclc_enhancer <- read.table("ICM.EPI.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("ICM.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()



tclc_enhancer <- read.table("Untreat.EPI.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("Untreat.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()


tclc_enhancer <- read.table("auxin-2days.EPI.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("auxin-2days.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()


tclc_enhancer <- read.table("washoff-2days.EPI.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(20)

pdf("washoff-2days.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.1))
dev.off()

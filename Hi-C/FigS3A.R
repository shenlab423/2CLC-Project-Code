
tclc_enhancer <- read.table("ESC.Nanog.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("white","orange","red"))(50)

pdf("ESC.Nanog.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(0,5,0.1))
dev.off()

tclc_enhancer <- read.table("2CLC.Nanog.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("white","orange","red"))(50)

pdf("2CLC.Nanog.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(0,5,0.1))
dev.off()



tclc_enhancer <- read.table("ESC.Klf4.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("white","orange","red"))(50)

pdf("ESC.Klf4.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(0,5,0.1))
dev.off()

tclc_enhancer <- read.table("2CLC.Klf4.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("white","orange","red"))(50)

pdf("2CLC.Klf4.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(0,5,0.1))
dev.off()



tclc_enhancer <- read.table("ESC.Sox2.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("white","orange","red"))(50)

pdf("ESC.sox2.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(0,5,0.1))
dev.off()

tclc_enhancer <- read.table("2CLC.Sox2.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("white","orange","red"))(50)

pdf("2CLC.sox2.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(0,5,0.1))
dev.off()
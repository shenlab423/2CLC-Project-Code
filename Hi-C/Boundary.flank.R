
tclc_enhancer <- read.table("ESC.boundary.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(60)

pdf("ESC.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-0.6,0.6,0.02))
dev.off()


tclc_enhancer <- read.table("2CLC.boundary.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(60)

pdf("2CLC.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-0.6,0.6,0.02))
dev.off()



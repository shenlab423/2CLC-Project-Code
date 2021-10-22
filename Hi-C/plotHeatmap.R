

tclc_enhancer <- read.table("2CLC.compartment.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(100)

pdf("2CP-2CLC.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.02))
dev.off()



tclc_enhancer <- read.table("ESC.compartment.txt",sep=" ")
library(gplots)
rampCol1 <- colorRampPalette(c("blue","white","red"))(100)

pdf("2CN-ESC.pdf",width=6,height=6)
heatmap.2(as.matrix(log2(tclc_enhancer)),trace = "none",col = rampCol1,Rowv = F,Colv = F,density.info = "none",breaks = seq(-1,1,0.02))
dev.off()


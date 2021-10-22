################## putative 2C enhancers
x <- read.table("H3K27ac_2C_UP_rpkm.bed.annotated.txt",sep="\t")
x2 <- x[x$V8=="Intergenic",]
x3 <- x2[abs(x2$V10)>3000,]
write.table(x3[,c(2,3,4)],"H3K27ac_2C_UP_rpkm.intergenic.bed",col.names=F,row.names=F,sep="\t",quote=F)

################## putative 2C enhancers (fold change > 4)
x <- read.table("H3K27ac_2C_UP_rpkm_FC4.bed.annotated.txt",sep="\t")
x2 <- x[x$V8=="Intergenic",]
x3 <- x2[abs(x2$V10)>3000,]
write.table(x3[,c(2,3,4)],"H3K27ac_2C_UP_rpkm_FC4.intergenic.bed",col.names=F,row.names=F,sep="\t",quote=F)

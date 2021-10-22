x <- read.table("ESC.PC1.bedGraph")
x <- na.omit(x)
schr <- paste("chr",c(1:19,"X","Y"),sep="")
out <- c()
for(chr in schr){
	tmp <- x[x$V1==chr,]
	tmp <- tmp[order(tmp$V4,decreasing=T),]
    out <- rbind(out,tmp)
}
out <- data.frame(out)
write.table(out[,c(1,2,3,4)],"Order.Compartment.ESC.bed",sep="\t",quote=F,row.names=F,col.names=F)



x <- read.table("2CLC.PC1.bedGraph")
x <- na.omit(x)
schr <- paste("chr",c(1:19,"X","Y"),sep="")
out <- c()
for(chr in schr){
	tmp <- x[x$V1==chr,]
	tmp <- tmp[order(tmp$V4,decreasing=T),]
    out <- rbind(out,tmp)
}
out <- data.frame(out)
write.table(out[,c(1,2,3,4)],"Order.Compartment.2CLC.bed",sep="\t",quote=F,row.names=F,col.names=F)


x <- read.table("ESC.PC1.bedGraph")
y <- read.table("2CLC.PC1.bedGraph")

s <- cbind(x,y)
s <- na.omit(s)
s <- s[,c(1,2,3,4,8)]
colnames(s) <- c("chr","start","end","ESC","2CLC")
atob <- s[s$ESC>0 & s$`2CLC`<0,]
btoa <- s[s$ESC<0 & s$`2CLC`>0,]

dim(atob)
dim(btoa)
write.table(atob[,c(1:3)],"AtoB.bed",col.names = F,row.names = F,sep="\t",quote=F)
write.table(btoa[,c(1:3)],"BtoA.bed",col.names = F,row.names = F,sep="\t",quote=F)


#################################

x <- read.table("10wESC.PC1.bedGraph")
x <- na.omit(x)
schr <- paste("chr",c(1:19,"X","Y"),sep="")
out <- c()
for(chr in schr){
  tmp <- x[x$V1==chr,]
  tmp <- tmp[order(tmp$V4,decreasing=T),]
  out <- rbind(out,tmp)
}
out <- data.frame(out)
write.table(out[,c(1,2,3,4)],"Order.Compartment.10wESC.bed",sep="\t",quote=F,row.names=F,col.names=F)



x <- read.table("10w2CLC.PC1.bedGraph")
x <- na.omit(x)
schr <- paste("chr",c(1:19,"X","Y"),sep="")
out <- c()
for(chr in schr){
  tmp <- x[x$V1==chr,]
  tmp <- tmp[order(tmp$V4,decreasing=T),]
  out <- rbind(out,tmp)
}
out <- data.frame(out)
write.table(out[,c(1,2,3,4)],"Order.Compartment.10w2CLC.bed",sep="\t",quote=F,row.names=F,col.names=F)


#################################

x <- read.table("4w-2Cn.PC1.bedGraph")
x <- na.omit(x)
schr <- paste("chr",c(1:19,"X","Y"),sep="")
out <- c()
for(chr in schr){
  tmp <- x[x$V1==chr,]
  tmp <- tmp[order(tmp$V4,decreasing=T),]
  out <- rbind(out,tmp)
}
out <- data.frame(out)
write.table(out[,c(1,2,3,4)],"Order.Compartment.4w-2Cn.bed",sep="\t",quote=F,row.names=F,col.names=F)



x <- read.table("4w-2Cp.PC1.bedGraph")
x <- na.omit(x)
schr <- paste("chr",c(1:19,"X","Y"),sep="")
out <- c()
for(chr in schr){
  tmp <- x[x$V1==chr,]
  tmp <- tmp[order(tmp$V4,decreasing=T),]
  out <- rbind(out,tmp)
}
out <- data.frame(out)
write.table(out[,c(1,2,3,4)],"Order.Compartment.4w-2Cp.bed",sep="\t",quote=F,row.names=F,col.names=F)

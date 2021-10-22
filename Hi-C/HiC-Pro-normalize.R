args <- commandArgs(T)
filename <- args[1]
filename2 <- paste(filename,".norm.matrix",sep='')
#print(filename2)
x <- read.table(filename)
#print(sum(x$V3))
x$V3 <- x$V3 * 100000000 / sum(x$V3)
#print(sum(x$V3))
x <- x[x$V3>0,]
write.table(x,filename2,sep='\t',col.names=F,row.names=F,quote=F)

files <- list.files(".",pattern="*.insulation$")

for(file in files){
x <- read.table(file,header=T)


a <- c()
for (i in 1:dim(x)[1]) {
  if (i>=100 & i+100 <= dim(x)[1]) {
    splices <- x[(i-100):(i+100),]
  }
  else{
    if (i+100 > dim(x)[1]){
      splices <- x[(i-100):dim(x)[1],]
    }
    else{
      splices <- x[1:(i+100),]
      
    }
  }
  #print(i)
  num = mean(splices$smoothedInsulaton,na.rm=T)
  a[i] <- num
}

x$mean <- a
x$IS <- log2(x$smoothedInsulaton / x$mean)
chrom <- sub(".matrix--is500001--nt0--ids260001--ss1--immean.insulation","",file)
out <- data.frame(chrom=chrom,start=x$start,end=x$end,IS=x$IS)
write.table(out,paste(file,".out",sep=""),col.names = F,row.names = F,sep="\t",quote=F)
}

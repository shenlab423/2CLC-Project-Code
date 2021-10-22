 m <- read.table("early2Cell.Loopstrength.txt")
 n <- read.table("ICM.Loopstrength.txt")
 x <- read.table("ESC.Loopstrength.txt")
 z <- read.table("2CLC.Loopstrength.txt")

 s <- cbind(m,n,z,x)
 pdf("Loopstrength.pdf",width = 5,height = 5)
 boxplot(log2(s),outline=F,names=c("2Cell","ICM","2CLC","ESC"))
 dev.off()

#t.test(x$V1,z$V1)
wilcox.test(log2(x$V1),log2(z$V1),alternative = "greater") # where y and x are numeric
wilcox.test(log2(n$V1),log2(m$V1),alternative = "greater") # where y and x are numeric


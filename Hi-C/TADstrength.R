x <- read.table("ESC.tadstrength.txt")
y <- read.table("2CLC.tadstrength.txt")

#m <- read.table("/oldnest/yezhang_zhu/Project/2C-nethome/HiC/embryo/RAWDATA/Aggregate/TAD2/2Cell.tadstrength.txt")
#n <- read.table("/oldnest/yezhang_zhu/Project/2C-nethome/HiC/embryo/RAWDATA/Aggregate/TAD2/ICM.tadstrength.txt")

pdf("TAD.strength.pdf",width = 5,height = 5)
boxplot(log2(y$V1),log2(x$V1),names=c("2CLC","ESC"),outline = F)
dev.off()

wilcox.test(log2(x$V1),log2(y$V1),paired = T)




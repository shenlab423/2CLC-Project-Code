

x <- read.table("ESC.EPIstrength.txt")
y <- read.table("2CLC.EPIstrength.txt")
boxplot(log2(y$V1),log2(x$V1),outline = F,names = c("2CLC","ESC"))


m <- read.table("early2Cell.EPIstrength.txt")
n <- read.table("ICM.EPIstrength.txt")

pdf("Enhancer.Loop.strength.pdf",width = 5,height = 5)
boxplot(log2(m$V1),log2(n$V1),log2(y$V1),log2(x$V1),outline = F,names=c("early2Cell","ICM","2CLC","ESC"))
dev.off()

wilcox.test(log2(m$V1),log2(n$V1),paired = T)
wilcox.test(log2(x$V1),log2(y$V1),paired = T)


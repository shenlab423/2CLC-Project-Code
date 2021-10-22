esc <- read.table("ESC.interTAD.txt")
tc <- read.table("2CLC.interTAD.txt")

esc <- esc[esc$V2 < 10000000,]
tc <- tc[tc$V2 < 10000000,]

pdf("Inter.TAD.pdf",width=5,height=5)
boxplot(esc$V1,tc$V1,names=c("ESC","2CLC"),outline=F)
dev.off()

wilcox.test(esc$V1,tc$V1,paired = T)


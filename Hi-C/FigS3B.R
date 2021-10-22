esc <- c(11.7872528151109, 23.4541957293347, 16.15600518)
tc <- c(6.72286868881754,11.63693734976075, 7.482516957)

rate <- tc/esc

x <- data.frame(Group=c("ESC","2CLC","ESC","2CLC","ESC","2CLC"),gene=c("Nanog","Nanog","Sox2","Sox2","Klf4","Klf4"),value=c(1,rate[1],1,rate[2],1,rate[3]))

x$Group <- factor(x$Group,levels=c("ESC","2CLC"))
x$gene <- factor(x$gene,levels=c("Nanog","Sox2","Klf4"))
library(ggplot2)
pdf("FigS3B.pdf",width = 5,height = 5)
ggplot(x,aes(x=gene,y=value,fill=Group)) + geom_bar(stat="identity",position=position_dodge(0.7),width=0.5) + theme_classic()
dev.off()

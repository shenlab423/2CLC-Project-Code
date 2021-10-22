x1 <- read.table("Allgenes.txt")
x2 <- read.table("Allgenes_CTCF.txt")

x1 <- x1[,c(1,2,3,4,5,6,7,8)]
x2 <- x2[,c(1,2,3,4)]
x1$Gene <- rownames(x1)
x2$Gene <- rownames(x2)

table <- merge(x1,x2,by="Gene")

genes <- c("Ctcf","Smc1a","Smc3","Rad21","Yy1")
table2 <- table[table$Gene %in% genes,]

rownames(table2) <- table2$Gene
table2 <- table2[,-1]

x <- c("Neg1","Neg1","Neg1","Neg1","Pos1","Pos1","Pos1","Pos1","Neg2","Neg2","Pos2","Pos2")

table3 <- rbind(table2,x)

library(reshape2)
table4 <- t(table3)
table5 <- melt(table4)

table6 <- table5[c(1:60),]
table7 <- table5[c(61:72),]
table7 <- table7[,c(1,3)]

table8 <- merge(table6,table7,by="Var1")
table8$value.x <- log2(as.numeric(table8$value.x))
library(dplyr)

table9 <- table8 %>% group_by(value.y,Var2) %>% summarise(meanValue=mean(value.x),sdValue=sd(value.x))


library(ggplot2)

table9$value.y <- factor(table9$value.y,levels = c("Neg1","Pos1","Neg2","Pos2"))
ggplot(table9,aes(x=Var2,y=meanValue,fill=value.y)) + geom_bar(stat = "identity", color = "black", width = 0.6, position = position_dodge()) + 
  geom_errorbar(aes(ymin=meanValue-sdValue, ymax=meanValue+sdValue),width=0.2,position=position_dodge(0.6)) + theme_classic() + ylab("Log2 Expression") + xlab("")
  

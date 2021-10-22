esc <- read.table("4w-2Cn.compartment.5.txt")
tc <- read.table("4w-2Cp.compartment.5.txt")

esc1 <- (esc[1,1])
tc1 <- (tc[1,1] )

esc <- read.table("10wESC.compartment.5.txt")
tc <- read.table("10w2CLC.compartment.5.txt")

esc2 <- (esc[1,1])
tc2 <- (tc[1,1])

x1 <- data.frame(type=c("rep1","rep1","rep2","rep2"),group=c("ESC","2CLC","ESC","2CLC"),values=c(esc1,tc1,esc2,tc2),type2=c("A-A","A-A","A-A","A-A"))

esc <- read.table("4w-2Cn.compartment.5.txt")
tc <- read.table("4w-2Cp.compartment.5.txt")

esc1 <- (esc[5,5])
tc1 <- (tc[5,5])

esc <- read.table("10wESC.compartment.5.txt")
tc <- read.table("10w2CLC.compartment.5.txt")

esc2 <- (esc[5,5])
tc2 <- (tc[5,5])

x2 <- data.frame(type=c("rep1","rep1","rep2","rep2"),group=c("ESC","2CLC","ESC","2CLC"),values=c(esc1,tc1,esc2,tc2),type2=c("B-B","B-B","B-B","B-B"))


esc <- read.table("4w-2Cn.compartment.5.txt")
tc <- read.table("4w-2Cp.compartment.5.txt")

esc1 <- (esc[5,1])
tc1 <- (tc[5,1])

esc <- read.table("10wESC.compartment.5.txt")
tc <- read.table("10w2CLC.compartment.5.txt")

esc2 <- (esc[5,1])
tc2 <- (tc[5,1])

x3 <- data.frame(type=c("rep1","rep1","rep2","rep2"),group=c("ESC","2CLC","ESC","2CLC"),values=c(esc1,tc1,esc2,tc2),type2=c("A-B","A-B","A-B","A-B"))

x <- rbind(x1,x2,x3)

library(ggplot2)
library(dplyr)

x2 <- x %>% group_by(group,type2) %>% summarise(mean=mean(values),minvalue=min(values),maxvalue=max(values))
x2$group <- factor(x2$group,levels=c("ESC","2CLC"))
x2$type2 <- factor(x2$type2,levels=c("A-A","B-B","A-B"))
ggplot(x2,aes(x=type2,y=log2(mean),fill=group,col=group)) + geom_errorbar(aes(ymin=log2(minvalue), ymax=log2(maxvalue)), width=0.1) + coord_cartesian(ylim = c(-1,1)) + theme_classic()


ggplot(x2, aes(x=type2, weight=log2(mean), ymin=log2(minvalue), ymax=log2(maxvalue), fill=group)) +
  geom_bar      (position=position_dodge(0.7), aes(y=log2(mean)), stat="identity",width=0.5) +
  geom_errorbar (position=position_dodge(width=0.7), colour="black",width=0.3) +
  theme_classic()



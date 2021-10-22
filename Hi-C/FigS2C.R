x <- read.table("ESC.PC1.bedGraph")
y <- read.table("2CLC.PC1.bedGraph")

mydata <- data.frame(x=x$V4,y=y$V4)
mydata <- na.omit(mydata)


atob <- mydata[mydata$x>0 & mydata$y<0,]
btoa <- mydata[mydata$x<0 & mydata$y>0,]

s <- c(55,368)/5321
data <- data.frame(type=c("AtoB","BtoA"),value=s)
library(ggplot2)
ggplot(data,aes(x=type,y=value,fill=type)) + geom_bar(stat="identity",width = 0.5) + theme_classic()

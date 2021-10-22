#############

x <- read.table("mESC_Dixon2012-raw_TADs")
x2 <- x[,7:206]
x3 <- x[,1:6]
insulation2c <- x2[,150]
insulationesc <- x2[,50]

#plot(insulationesc,insulation2c,pch=19,xlim=c(-0.3,1.2),ylim=c(-0.3,1.2),xlab = "Insulation score(ESC)",ylab="Insulation score(2CLC)")
#abline(a = 0 , b = 1, col = "black")
library(LSD)

heatscatter(insulationesc,insulation2c,xlim=c(-0.3,1.2),ylim=c(-0.3,1.2))
abline(a = 0 , b = 1, col = "black")

boxplot(insulationesc,insulation2c,outline = F,names=c("ESC","2CLC"),ylab="Insulation Score")
wilcox.test(insulationesc,insulation2c,paired = T)


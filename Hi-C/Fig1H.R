esc <- read.table("10wESC.PS.2.txt")
tc <- read.table("10w2CLC.PS.2.txt")

esc <- esc[-1,]
tc <- tc[-1,]

esc <- esc[-1,]
tc <- tc[-1,]

esc$V4 <- esc$V3 / sum(esc$V3) 
tc$V4 <- tc$V3 / sum(tc$V3) 
esc$x <- c(1:14)
tc$x <- c(1:14)

esc$type <- "ESC"
tc$type <- "2CLC"


data <- rbind(esc,tc)


esc <- read.table("4w-2Cn-facs.PS.2.txt")
tc <- read.table("4w-2Cp-facs.PS.2.txt")

esc <- esc[-1,]
tc <- tc[-1,]

esc <- esc[-1,]
tc <- tc[-1,]

esc$V4 <- esc$V3 / sum(esc$V3) 
tc$V4 <- tc$V3 / sum(tc$V3) 
esc$x <- c(1:14)
tc$x <- c(1:14)

esc$type <- "ESC"
tc$type <- "2CLC"


data2 <- rbind(esc,tc)

data <- rbind(data,data2)

library(ggplot2)
#ggplot(data,aes(x=x,y=(V4),col=type)) + geom_smooth(method = loess, formula = y ~ x, span=0.05) +
# scale_x_continuous(breaks= seq(0,14,by=1), labels = c(0:14),limits = c(0,14), expand = c(0,0))  + theme_classic()

#ggplot(data,aes(x=x,y=(V4),col=type)) + geom_line(size=0.5) +
#  scale_x_continuous(breaks= seq(0,14,by=1), labels = c(0:14),limits = c(0,14), expand = c(0,0))  + theme_classic()

library(dplyr)
data2 <- data %>% group_by(x,type) %>% summarise(lower=min(V4),upper=max(V4),ave=mean(V4))


p<-ggplot(data=data2, aes(x=x, y=ave, colour=type)) + geom_line() #+ scale_color_manual(values = c("white","white"))
p<-p+geom_ribbon(aes(ymin=data2$lower, ymax=data2$upper), fill = "grey70",alpha=0.5) + scale_x_continuous(breaks= seq(0,14,by=1), labels = c(0:14),limits = c(0,14), expand = c(0,0))  + theme_classic()

p
#+ coord_cartesian(ylim = c(-10,-5),xlim=c(5,8))
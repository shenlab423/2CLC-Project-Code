data <- read.table("enrich.txt")
data$enrich <- data$V4 / data$V1
data$type <- c("Down","Up","Nochange","Top10","Top30")

data <- data[!data$type %in% c("Down10"),]
library(ggplot2)
data$type <- factor(data$type,levels=c("Nochange","Down","Up","Top30","Top10"))
ggplot(data,aes(x=type,y=enrich)) + geom_bar(stat="identity",width = 0.5) + theme_classic()

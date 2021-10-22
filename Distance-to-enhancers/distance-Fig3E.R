
down <- read.table("Down_promoter_nearest_2Cenhancer.txt",stringsAsFactors=F)
down <- down[down$V1 == down$V7,]

up <- read.table("Up_promoter_nearest_2Cenhancer.txt",stringsAsFactors=F)
up <- up[up$V1 == up$V7,]

nochange <- read.table("Nochange_promoter_nearest_2Cenhancer.txt",stringsAsFactors=F)
nochange <- nochange[nochange$V1 == nochange$V7,]


library(dplyr)
down$type <- "Down"
down <- down[,c(10,11)]

up$type <- "Up"
up <- up[,c(10,11)]

nochange$type <- "Nochange"
nochange <- nochange[,c(10,11)]

#all <- rbind(up,down,nochange)

up$distance <- log10(abs(up$V10)+1)
down$distance <- log10(abs(down$V10)+1)
nochange$distance <- log10(abs(nochange$V10)+1)

svalue <- function(x){
  s <- hist(x,breaks = seq(0,8,0.1),xlim = c(0,8))
  return((s$breaks[1:length(s$breaks)-1]))
}
svalue2 <- function(x){
  s <- hist(x,breaks = seq(0,8,0.1),xlim = c(0,8))
  return(s$counts)
}

up2 <- data.frame(type="Up",Group=svalue(up$distance),n=svalue2(up$distance))
down2 <- data.frame(type="Down",Group=svalue(down$distance),n=svalue2(down$distance))
nochange2 <- data.frame(type="Nochange",Group=svalue(nochange$distance),n=svalue2(nochange$distance))

s <- merge(up2,down2,by="Group",all = T)
all <- merge(s,nochange2,by="Group",all=T)
all <- data.frame(all$Group)

up2 <- merge(all,up2,by.x="all.Group",by.y = "Group",all.x = T)
down2 <- merge(all,down2,by.x="all.Group",by.y = "Group",all.x = T)
nochange2 <- merge(all,nochange2,by.x="all.Group",by.y = "Group",all.x = T)

up2$type <- "Up"
down2$type <- "Down"
nochange2$type <- "Nochange"


all <- rbind(up2,down2,nochange2)
colnames(all) <- c("distance","type","n")
all[is.na(all)] <- 0

all2 <- all %>% group_by(type) %>% mutate(n2 = cumsum(n))
all3 <- all2 %>% group_by(type) %>% summarise(total=max(n2))

all4 <- all2 %>% left_join(all3,by="type")
all4 <- all4 %>% mutate(fre=n2/total)

all4 <- data.frame(all4[,c("type","distance","fre")])

#max(all4$distance)
#systemadd <- data.frame(type=c("Down","Up","Nochange"),
#                        distance = c(8,8,8),
#                        fre = c(1,1,1),stringsAsFactors=F) 

#all4 <- data.frame(rbind(all4,systemadd))
all4$type <- factor(all4$type,levels=c("Up","Nochange","Down"))
library(ggplot2)
ggplot(all4,aes(x=distance,y=fre,fill=type,col=type)) + geom_line(size=1) + coord_cartesian(xlim = c(0,8))+ theme_classic() +
  xlab("log10 distance from TSS to the closest expressed Enhancer(bp)") + ylab("cumulative portion of genes")

#all4 <- all4[all4$distance<6.3,]

ks.test(jitter(all4$fre[all4$type=="Up"]),jitter(all4$fre[all4$type=="Down"]))     #4.271e-08
ks.test(jitter(all4$fre[all4$type=="Nochange"]),jitter(all4$fre[all4$type=="Down"]))  #2.928e-07
ks.test(jitter(all4$fre[all4$type=="Nochange"]),jitter(all4$fre[all4$type=="Up"]))  #2.928e-07






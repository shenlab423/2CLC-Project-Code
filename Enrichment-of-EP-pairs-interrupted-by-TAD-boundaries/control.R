down <- read.table("Down.close.2Cenhancer.len",stringsAsFactors=F)
down <- down[down$V1 == down$V7,]

up <- read.table("Up.close.2Cenhancer.len",stringsAsFactors=F)
up <- up[up$V1 == up$V7,]

nochange <- read.table("Nochange.close.2Cenhancer.len",stringsAsFactors=F)
nochange <- nochange[nochange$V1 == nochange$V7,]


top <- read.table("Up10.close.2Cenhancer.len",stringsAsFactors=F)
top <- top[top$V1 == top$V7,]

top30 <- read.table("Up30.close.2Cenhancer.len",stringsAsFactors=F)
top30 <- top30[top30$V1 == top30$V7,]



down <- down[abs(down$V10)<500000,]
library(dplyr)
down <- down %>% mutate(start2=ifelse(V10>0,V2-V10,V2),end2=ifelse(V10>0,V2,V2+abs(V10)))
down_out <- down[,c(1,11,12,5)]
down <- down %>% mutate(start=ifelse(V9>V2,V2,V9),end=ifelse(V9>V2,V9,V2))
down_out2 <- down[,c(1,13,14,5)]

down_out[down_out<0] <- 1
write.table(down_out,paste("Down.close.","Dixon",".control.txt",sep=""),col.names = F,row.names = F,sep = "\t",quote=F)
write.table(down_out2,paste("Down.close.","Dixon",".txt",sep=""),col.names = F,row.names = F,sep = "\t",quote=F)



up <- up[abs(up$V10)<500000,]
library(dplyr)
up <- up %>% mutate(start2=ifelse(V10>0,V2-V10,V2),end2=ifelse(V10>0,V2,V2+abs(V10)))
up_out <- up[,c(1,11,12,5)]
up <- up %>% mutate(start=ifelse(V9>V2,V2,V9),end=ifelse(V9>V2,V9,V2))
up_out2 <- up[,c(1,13,14,5)]

up_out[up_out<0] <- 1
write.table(up_out,paste("Up.close.","Dixon",".control.txt",sep=""),col.names = F,row.names = F,sep = "\t",quote=F)
write.table(up_out2,paste("Up.close.","Dixon",".txt",sep=""),col.names = F,row.names = F,sep = "\t",quote=F)


nochange <- nochange[abs(nochange$V10)<500000,]
library(dplyr)
nochange <- nochange %>% mutate(start2=ifelse(V10>0,V2-V10,V2),end2=ifelse(V10>0,V2,V2+abs(V10)))
nochange_out <- nochange[,c(1,11,12,5)]
nochange <- nochange %>% mutate(start=ifelse(V9>V2,V2,V9),end=ifelse(V9>V2,V9,V2))
nochange_out2 <- nochange[,c(1,13,14,5)]

nochange_out[nochange_out<0] <- 1
write.table(nochange_out,paste("Nochange.close.","Dixon",".control.txt",sep=""),col.names = F,row.names = F,sep = "\t",quote=F)
write.table(nochange_out2,paste("Nochange.close.","Dixon",".txt",sep=""),col.names = F,row.names = F,sep = "\t",quote=F)



top <- top[abs(top$V10)<500000,]
library(dplyr)
top <- top %>% mutate(start2=ifelse(V10>0,V2-V10,V2),end2=ifelse(V10>0,V2,V2+abs(V10)))
top_out <- top[,c(1,11,12,5)]
top <- top %>% mutate(start=ifelse(V9>V2,V2,V9),end=ifelse(V9>V2,V9,V2))
top_out2 <- top[,c(1,13,14,5)]

top_out[top_out<0] <- 1
write.table(top_out,paste("Top.close.","Dixon",".control.txt",sep=""),col.names = F,row.names = F,sep = "\t",quote=F)
write.table(top_out2,paste("Top.close.","Dixon",".txt",sep=""),col.names = F,row.names = F,sep = "\t",quote=F)


top30 <- top30[abs(top30$V10)<500000,]
library(dplyr)
top30 <- top30 %>% mutate(start2=ifelse(V10>0,V2-V10,V2),end2=ifelse(V10>0,V2,V2+abs(V10)))
top30_out <- top30[,c(1,11,12,5)]
top30 <- top30 %>% mutate(start=ifelse(V9>V2,V2,V9),end=ifelse(V9>V2,V9,V2))
top30_out2 <- top30[,c(1,13,14,5)]

top30_out[top30_out<0] <- 1
write.table(top30_out,paste("Top30.close.","Dixon",".control.txt",sep=""),col.names = F,row.names = F,sep = "\t",quote=F)
write.table(top30_out2,paste("Top30.close.","Dixon",".txt",sep=""),col.names = F,row.names = F,sep = "\t",quote=F)


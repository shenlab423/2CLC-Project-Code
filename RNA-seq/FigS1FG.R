
repeats <- read.table("mm9_repeats_simple.bed")
rnaesc <- read.table("Neg.RNA.repeats.txt")
rna2c <- read.table("Pos.RNA.repeats.txt")

library(dplyr)

rnaesc <- rnaesc %>% mutate(ESC=(V4+V8+V12+V16)/4)
rnaesc <- data.frame(rnaesc[,17])

rna2c <- rna2c %>% mutate(`2CLC`=(V4+V8+V12+V16)/4)
rna2c <- data.frame(rna2c[,17])

repeats_rna <- cbind(repeats,rnaesc,rna2c)

ervl <- repeats_rna[repeats_rna$V7=="ERVL",]

mervl <- ervl[grep("^M",ervl$V5),]

mervl_int <- ervl[grep("MERVL-int",ervl$V5),]

mervl_mm <- ervl[grep("MT2_Mm",ervl$V5),]

mervl_out <- rbind(mervl_int,mervl_mm)
write.table(mervl_out[,c(1,2,3,5,6,4)],"mm9-MERVL-int-MT2-Mm.bed",col.names = F,row.names = F,sep="\t",quote=F)


mervl <- mervl %>% mutate(type = "MERVL(n=9987)") %>% filter(rnaesc...17.>0 | rna2c...17.>0)
mervl_int <- mervl_int %>% mutate(type = "MERVL-int(n=1206)") %>% filter(rnaesc...17.>0 | rna2c...17.>0)
mervl_mm <- mervl_mm %>% mutate(type = "MT2_Mm(n=1477)") %>% filter(rnaesc...17.>0 | rna2c...17.>0)


all <- rbind(mervl,mervl_int,mervl_mm)
all <- all[,c("type","rnaesc...17.","rna2c...17.")]
colnames(all) <- c("type","ESC","2CLC")

library(reshape2)

all2 <- melt(all)

library(ggplot2)
all2$type <- factor(all2$type,levels=c("MERVL(n=9987)","MERVL-int(n=1206)","MT2_Mm(n=1477)"))
ggplot(all2,aes(x=type,y=value,fill=variable)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0,4)) + ylab("FPKM") + xlab("") + 
  theme_classic() + scale_fill_manual(values = c("blue","red"))

wilcox.test(mervl$rnaesc...17.,mervl$rna2c...17.,paired = T)
wilcox.test(mervl_int$rnaesc...17.,mervl_int$rna2c...17.,paired = T)
wilcox.test(mervl_mm$rnaesc...17.,mervl_mm$rna2c...17.,paired = T)


###################################################################

repeats <- read.table("mm9_repeats_simple.bed")
rnaesc <- read.table("Neg.RNA.repeats.txt")
rna2c <- read.table("Pos.RNA.repeats.txt")

library(dplyr)

rnaesc <- rnaesc %>% mutate(ESC=(V4+V8+V12+V16)/4)
rnaesc <- data.frame(rnaesc[,17])

rna2c <- rna2c %>% mutate(`2CLC`=(V4+V8+V12+V16)/4)
rna2c <- data.frame(rna2c[,17])

repeats_rna <- cbind(repeats,rnaesc,rna2c)


repeats_all2 <- repeats_rna %>% group_by(V5) %>% summarise(ESC_RNA=mean(rnaesc...17.),TCLC_RNA=mean(rna2c...17.))

repeats_all2 <- repeats_all2 %>% filter(ESC_RNA>1 | TCLC_RNA>1)
write.table(repeats_all2,"repeats_all2.txt",col.names = F,row.names = F,sep="\t",quote = F)
library(ggplot2)

g <- ggplot(repeats_all2, aes(x=log2(ESC_RNA+1), y=log2(TCLC_RNA+1),colour="type")) + geom_point() + scale_color_manual(values=c("black")) + theme(legend.position = "none") + scale_y_continuous(limits = c(0,4)) + scale_x_continuous(limits = c(0,4))
g <- g  + ylab("2CLC expression (log2RPKM)") + xlab("ESC expression(log2RPKM)") 

g <- g + geom_text(aes(label=V5), size=4, vjust = 0)  
g + theme_classic()+ theme(legend.position = "none")


###################################################################


repeats <- read.table("mm9_repeats_simple.bed")
rnaesc <- read.table("Neg.RNA.repeats.txt")
rna2c <- read.table("Pos.RNA.repeats.txt")

library(dplyr)

rnaesc <- rnaesc %>% mutate(ESC=(V4+V8+V12+V16)/4)
rnaesc <- data.frame(rnaesc[,17])

rna2c <- rna2c %>% mutate(`2CLC`=(V4+V8+V12+V16)/4)
rna2c <- data.frame(rna2c[,17])

repeats_rna <- cbind(repeats,rnaesc,rna2c)

ervl <- repeats_rna[repeats_rna$V7=="ERVL",]


mervl_int <- ervl[grep("MERVL-int",ervl$V5),]

mervl_mm <- ervl[grep("MT2_Mm",ervl$V5),]

mervl <- rbind(mervl_int,mervl_mm)


expressed3 <- mervl[mervl$rna2c...17.>10*mervl$rnaesc...17. & mervl$rna2c...17.>5,]
write.table(expressed3[,c(1,2,3,5,6,4)],"MERVL-int-MT2-Mm-expressed-FC10-2C5.bed",col.names = F,row.names = F,sep="\t",quote=F)

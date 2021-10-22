###############

all <- read.table("Allgenes.pvalue.txt")
all <- all[order(all$FC,decreasing = T),]
all <- all[all$FDR<0.05 & all$FC > 1.5,]

up30 <- all[1:313,]
up10 <- all[1:104,]

promoter <- read.table("mm9_promoter_2k2k_strand.bed")
up30 <- promoter[promoter$V5 %in% (up30$GeneName),]
up10 <- promoter[promoter$V5 %in% (up10$GeneName),]


up <- read.table("Upgenes_1.5.txt")
down <- read.table("Downgenes_1.5.txt")
nochange <- read.table("Nochangegenes.txt")

up <- promoter[promoter$V5 %in% up$V1,]
down <- promoter[promoter$V5 %in% down$V1,]
nochange <- promoter[promoter$V5 %in% (nochange$GeneName),]


write.table(up30,"Up30_promoter_enrich.bed",col.names = F,row.names = F,sep="\t",quote = F)
write.table(up10,"Up10_promoter_enrich.bed",col.names = F,row.names = F,sep="\t",quote = F)

write.table(up,"Up_promoter_enrich.bed",col.names = F,row.names = F,sep="\t",quote = F)
write.table(down,"Down_promoter_enrich.bed",col.names = F,row.names = F,sep="\t",quote = F)
write.table(nochange,"Nochange_promoter_enrich.bed",col.names = F,row.names = F,sep="\t",quote = F)

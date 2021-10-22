# get Up,Down,Non-regulated gene promoters

up <- read.table("Upgenes_1.5.txt")
down <- read.table("Downgenes_1.5.txt")
nochange <- read.table("Nochangegenes.txt")

promoter <- read.table("mm9_promoter_2k2k_strand.bed")
up2 <- promoter[promoter$V5 %in% up$V1,]
down2 <- promoter[promoter$V5 %in% down$V1,]
nochange2 <- promoter[promoter$V5 %in% (nochange$GeneName),]

up2$V2 <- up2$V2 + 1999
up2$V3 <- up2$V3 - 2000


down2$V2 <- down2$V2 + 1999
down2$V3 <- down2$V3 - 2000


nochange2$V2 <- nochange2$V2 + 1999
nochange2$V3 <- nochange2$V3 - 2000

write.table(up2,"Up_promoter.bed",col.names = F,row.names = F,sep="\t",quote = F)
write.table(down2,"Down_promoter.bed",col.names = F,row.names = F,sep="\t",quote = F)
write.table(nochange2,"Nochange_promoter.bed",col.names = F,row.names = F,sep="\t",quote = F)
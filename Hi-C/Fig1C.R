esc <- read.table("4w-2Cn.compartment.strength.txt")
tc <- read.table("4w-2Cp.compartment.strength.txt")

esc2 <- read.table("10wESC.compartment.strength.txt")
tc2 <- read.table("10w2CLC.compartment.strength.txt")

esc <- rbind(esc,esc2)
tc <- rbind(tc,tc2)

esc <- esc[!esc$V1%in%c("chrX","chrY"),]
tc <- tc[!tc$V1%in%c("chrX","chrY"),]

esc$value <- (esc$V2 + esc$V5) / (esc$V3 + esc$V4)
tc$value <- (tc$V2 + tc$V5) / (tc$V3 + tc$V4)

boxplot(esc$value,tc$value,outline = F,names = c("ESC","2CLC"))

wilcox.test(esc$value,tc$value,paired = T)

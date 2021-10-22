x1 <- read.table("Order.Compartment.ESC.bed")
x1 <- x1[x1$V1 == "chr1",]
barplot(x1$V4,ylim = c(-0.1,0.1))



x1 <- read.table("Order.Compartment.2CLC.bed")
x1 <- x1[x1$V1 == "chr1",]
barplot(x1$V4,ylim = c(-0.1,0.1))

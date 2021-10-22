up <- read.table("Upgenes_1.5.txt")
up <- unique(up)
rownames(up) <- up$V1
############################################
library(dplyr)
library(gplots)

table <- read.table("genes_CTCF.fpkm_table",header=T)
table <- table[,c(1,2,3,4,5,6)]
colnames(table) <- c("GeneSymbol","UT","1d","2d","4d","wash")

rownames(table) <- table$GeneSymbol


table2 <- read.table("genes_rad21.fpkm_table",header=T)
colnames(table2) <- c("GeneSymbol","12h","18h","1h","24h-1","24h-2","3h","6h","9h","ctrl1","ctrl2")

rownames(table2) <- table2$GeneSymbol
table2$ctrl <- 0.5 * (table2$ctrl1 + table2$ctrl2)
table2$`24h` <- 0.5 * (table2$`24h-1` + table2$`24h-2`)

table2 <- table2[,c("GeneSymbol","ctrl","1h","3h","6h","9h","12h","18h","24h")]


table_up <- table[table$GeneSymbol %in% rownames(up),]
table2_up <- table2[table2$GeneSymbol %in% rownames(up),]

table_up <- table_up[,-1]
table2_up <- table2_up[,-1]

plot1 <- table_up[,c(1,3,4)]
plot1$Change2d <- (plot1$`2d`+1) / (plot1$UT+1)
plot1$Change4d <- (plot1$`4d`+1) / (plot1$UT+1)

boxplot(log2(plot1[,c(4,5)]),outline = F)


out1 <- plot1[,c(4,5)]
out1$type <- "2Cspecific"

wilcox.test(log2(plot1$Change2d))
wilcox.test(log2(plot1$Change4d))


plot2 <- table2_up[,c(1,6,8)]
plot2$Change12h <- (plot2$`12h`+1) / (plot2$ctrl+1)
plot2$Change24h <- (plot2$`24h`+1) / (plot2$ctrl+1)

boxplot(log2(plot2[,c(4,5)]),outline = F)

out4 <- plot2[,c(4,5)]
out4$type <- "2Cspecific"

wilcox.test(log2(plot2$Change12h))
wilcox.test(log2(plot2$Change24h))


###################################################

all <- read.table("Allgenes.txt")
random <- all 

############################################
library(dplyr)
library(gplots)

table <- read.table("genes_CTCF.fpkm_table",header=T)
table <- table[,c(1,2,3,4,5,6)]
colnames(table) <- c("GeneSymbol","UT","1d","2d","4d","wash")

rownames(table) <- table$GeneSymbol


table2 <- read.table("genes_rad21.fpkm_table",header=T)
colnames(table2) <- c("GeneSymbol","12h","18h","1h","24h-1","24h-2","3h","6h","9h","ctrl1","ctrl2")

rownames(table2) <- table2$GeneSymbol
table2$ctrl <- 0.5 * (table2$ctrl1 + table2$ctrl2)
table2$`24h` <- 0.5 * (table2$`24h-1` + table2$`24h-2`)

table2 <- table2[,c("GeneSymbol","ctrl","1h","3h","6h","9h","12h","18h","24h")]


table_random <- table[table$GeneSymbol %in% rownames(random),]
table2_random <- table2[table2$GeneSymbol %in% rownames(random),]

table_random <- table_random[,-1]
table2_random <- table2_random[,-1]


plot1 <- table_random[,c(1,3,4)]
plot1$Change2d <- (plot1$`2d`+1) / (plot1$UT+1)
plot1$Change4d <- (plot1$`4d`+1) / (plot1$UT+1)

boxplot(log2(plot1[,c(4,5)]),outline = F)

out2 <- plot1[,c(4,5)]
out2$type <- "All"

wilcox.test(log2(plot1$Change2d))
wilcox.test(log2(plot1$Change4d))


plot2 <- table2_random[,c(1,6,8)]
plot2$Change12h <- (plot2$`12h`+1) / (plot2$ctrl+1)
plot2$Change24h <- (plot2$`24h`+1) / (plot2$ctrl+1)

boxplot(log2(plot2[,c(4,5)]),outline = F)

wilcox.test(log2(plot2$Change12h))
wilcox.test(log2(plot2$Change24h))

out5 <- plot2[,c(4,5)]
out5$type <- "All"


############################################

up <- read.table("plurinet")
rownames(up) <- up$V1

############################################
library(dplyr)
library(gplots)

table <- read.table("genes_CTCF.fpkm_table",header=T)
table <- table[,c(1,2,3,4,5,6)]
colnames(table) <- c("GeneSymbol","UT","1d","2d","4d","wash")

rownames(table) <- table$GeneSymbol


table2 <- read.table("genes_rad21.fpkm_table",header=T)
colnames(table2) <- c("GeneSymbol","12h","18h","1h","24h-1","24h-2","3h","6h","9h","ctrl1","ctrl2")

rownames(table2) <- table2$GeneSymbol
table2$ctrl <- 0.5 * (table2$ctrl1 + table2$ctrl2)
table2$`24h` <- 0.5 * (table2$`24h-1` + table2$`24h-2`)

table2 <- table2[,c("GeneSymbol","ctrl","1h","3h","6h","9h","12h","18h","24h")]


table_up <- table[table$GeneSymbol %in% rownames(up),]
table2_up <- table2[table2$GeneSymbol %in% rownames(up),]

table_up <- table_up[,-1]
table2_up <- table2_up[,-1]

plot1 <- table_up[,c(1,3,4)]
plot1$Change2d <- (plot1$`2d`+1) / (plot1$UT+1)
plot1$Change4d <- (plot1$`4d`+1) / (plot1$UT+1)

boxplot(log2(plot1[,c(4,5)]),outline = F)

out3 <- plot1[,c(4,5)]
out3$type <- "pluri"

wilcox.test(log2(plot1$Change2d))
wilcox.test(log2(plot1$Change4d))


plot2 <- table2_up[,c(1,6,8)]
plot2$Change12h <- (plot2$`12h`+1) / (plot2$ctrl+1)
plot2$Change24h <- (plot2$`24h`+1) / (plot2$ctrl+1)

boxplot(log2(plot2[,c(4,5)]),outline = F)

out6 <- plot2[,c(4,5)]
out6$type <- "pluri"

wilcox.test(log2(plot2$Change12h))
wilcox.test(log2(plot2$Change24h))


###########

out <- rbind(out1,out2,out3)
library(reshape2)
out <- melt(out)
library(ggplot2)
ggplot(out,aes(x=type,y=log2(value),fill=variable)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim=c(-1.5,2.5)) + theme_classic()


out <- rbind(out4,out5,out6)
library(reshape2)
out <- melt(out)
library(ggplot2)
ggplot(out,aes(x=type,y=log2(value),fill=variable)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim=c(-3,2.5)) + theme_classic()

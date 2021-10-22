x1 <- read.table("Allgenes.txt")
x2 <- read.table("Allgenes_CTCF.txt")

x1$Gene <- rownames(x1)
x2$Gene <- rownames(x2)
x3 <- merge(x1[,c("Gene","FC")],x2[,c("Gene","FC")],by="Gene")

library(ggplot2)
g <- ggplot(x3, aes(x=(FC.x), y=(FC.y))) + geom_point() + theme(legend.position = "none") + scale_y_continuous(limits = c(-5,12)) + scale_x_continuous(limits = c(-5,12))
g <- g  + ylab("2CLC/ESC (CTCF-deplete)") + xlab("2CLC/ESC (WT)") #+ geom_abline(intercept=log2(1),slope=1,colour="#990000", linetype="dashed")
g + theme_classic()+ theme(legend.position = "none")

cor.test(x3$FC.x,x3$FC.y,method = "pearson")   #cor = 0.7465889   p-value < 2.2e-16



highlight <- c("Dub1",
               "Zscan4f",
               "Zscan4c",
               "Zscan4d",
               "Tdpoz3",
               "Tdpoz4",
               "Zfp352",
               "Pou5f1","Nanog","Myc","Zfp42")

table1_in <- x3[(x3$Gene) %in% highlight,]
table1_in$PNAME <- (table1_in$Gene)

table1_nin <- x3[!(x3$Gene) %in% highlight,]
table1_nin$PNAME <- ""

x3 <- rbind(table1_in,table1_nin)

library(ggplot2)
g <- ggplot(table1_in, aes(x=(FC.x), y=(FC.y))) + geom_point()  + theme(legend.position = "none") + scale_y_continuous(limits = c(-5,12)) + scale_x_continuous(limits = c(-5,12))
g <- g  + ylab("2CLC/ESC (CTCF-deplete)") + xlab("2CLC/ESC (WT)") 
g <- g + geom_text(aes(label=PNAME), size=4, vjust = 0)  
g + theme_classic()+ theme(legend.position = "none")



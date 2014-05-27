# remake the cluster plots with the newly organized cluster data, sigh.  
# 052614
library(ggplot2)
library(reshape2)
library(scales)


setwd("~/Desktop/ShadePaper/clustering/")

somm <- read.table("kohonen/ILSunSh_shortDEshade_no622noSPE_LogFC_SOM8x8_melt_052614.txt", sep = " ")

levels(somm$chr)
levels(somm$geno)
somm$chr <- factor(somm$chr, levels(somm$chr)[c(13, 1, 5, 6, 7, 8, 9, 10, 11, 12, 2, 3, 4)])
somm$geno <- factor(somm$geno, levels(somm$geno)[c(73, 1:7, 25:72, 8:24)])
head(somm)

somo <- somm[order(c(somm$cluster, somm$number), decreasing = F), ]
head(somo)
tail(somo)
summary(somo)
somo <- na.omit(somo)
head(somo)
tail(somo)
dim(somo)


####plot of gene expression for each gene, listed by cluster####
clusterz <- NULL
for (i in 1:64) {
  clusterz$cluster[i] <- i
  clusterz$length[i] <- length(unique(somm$itag[somm$cluster == i]))
}
clusterz <- as.data.frame(clusterz)
head(clusterz) # how long is the first cluster? 222
clusterz$total[1] <- 222 
for (i in 2:64) clusterz$total[i] <- clusterz$length[i] + clusterz$total[i-1] 
for (i in 1:64) clusterz$breaks <- round(clusterz$total - (0.5*clusterz$length))
head(clusterz)
tail(clusterz)

write.csv(clusterz, "ILSunSh_KohonenClusterAxisBreaks.csv")
clusterz <- read.csv("ILSunSh_KohonenClusterAxisBreaks.csv", row.names = 1)


library(ggplot2) 


##### plot gene expression in order of cluster by IL #####
  p <- ggplot(somo, aes(x = geno, y = number, fill = logFC, color = logFC)) +
    geom_tile() +
    theme_bw() +
    theme(aspect.ratio = 2) + 
    scale_fill_gradientn(colours = c("green", "white", "darkorchid"), 
                         limits = c(-1, 1), oob = squish) +
    scale_colour_gradientn(colours = c("green", "white", "darkorchid")) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0.8, size = 6)) +
    #theme(axis.text.y = element_blank()) +
    theme(axis.text.y = element_text(size = 6)) +
    #theme(panel.grid.major = element_blank()) +
    scale_y_discrete(breaks = clusterz$breaks, labels = clusterz$cluster) +
    ylab("Clustered Gene Expression") +
    xlab("Introgression Lines") +
    labs(fill = "Log2 Fold Change\n(Shade/Sun)\nExpression", 
         title = "Shade-responsive Gene Expression in each IL") +
    guides(colour = "none")
  ggsave("orderClusterTest_white.pdf")


####PLOT OF GENE EXPRESSION BY CLUSTER BY IL#####
p <- ggplot(ccmo, aes(x = geno, y = cluster, fill = logFC, color = logFC)) +
  geom_tile() +
  theme_bw() +
  theme(aspect.ratio = 2) + 
  scale_fill_gradientn(colours = c("blue", "black", "yellow"), limits = c(-2.5, 2.5), oob = squish) +
  scale_colour_gradientn(colours = c("blue", "black", "yellow")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0.8, size = 6)) +
  theme(axis.text.y = element_text(size = 7)) +
  ylab("Clustered Gene Expression") +
  xlab("Introgression Lines") +
  labs(fill = "Log2 Fold Change\n(Shade/Sun)\nExpression", 
       title = "Clustered Shade-responsive Gene Expression in each IL") +
  guides(colour = "none")
ggsave("clusterTestBiggerScale_black.pdf")
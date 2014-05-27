# self-organizing maps!  using HOPACH package from Bioconductor
# April 2014
  #using kohonen:
  # start with IL*trt logFC genes (>17700)
  # use the log(ILsha/ILsun) value (aka ILsha + M82sha - ILsun)
  # reduce gene number by getting rid of genes with low variance across ILs (unlog value first)
  # reduce gene number by including 8770 genes DE in trtsha (no geno factor)
# exclude PENNELLII!!!!

#hopach

library(hopach)

log.data.cv.sh <- read.csv("ILSunShDE_log2Diff_041714.csv", row.names = 1)

head(log.data.cv.sh)
summary(log.data.cv.sh)

log.data.cv.sh <- log.data.cv.sh[, names(log.data.cv.sh) != "SPE"]

gene.dist <- distancematrix(log.data.cv.sh, "cor") ## shockingly fast! what about "cor"???
dim(gene.dist)
head(gene.dist)[, 1:5]

hobj <- hopach(log.data.cv.sh, dmat = gene.dist) 
names(hobj)

table(hobj$clust$labels)
hobj$clust$k # 182 clusters, sigh  #only 46 with cor distance matrix!!

setwd("~/Desktop/ShadePaper/clustering/hopach/")

png("hopach_COR_clusters_noSPE_052014.png")
dplot(dist = gene.dist, hopachobj = hobj, ord = "final", 
      main = "HOPACH cor clusters, log2(ILsha/ILsun)", showclusters = T, )
dev.off()

#####

gd <- as.data.frame(as.matrix(gene.dist))
head(gd)[, 1:5]
names(gd) <- rownames(log.data.cv.sh)
rownames(gd) <- rownames(log.data.cv.sh)
head(gd)[, 1:5]

gd$order <- hobj$final$order
#gd$numord <- seq(1:8352) # don't need this
gd$cl.labs <- hobj$final$labels
dim(gd)
gd <- gd[, c(8354, 8353, 1:8352)]
gd[1:5, 1:5]
gd <- gd[order(gd$cl.labs), ]
gd[8349:8352, 1:5]
#unique(substring(as.character(gd$cl.labs), 1, 6)) # so the first 6 digits of label = 82 clusters!
#gd$cluster <- factor(substring(as.character(gd$cl.labs), 1, 6))

names(table(hobj$clust$labels))
unique(round(gd$cl.labs/10e10)) #1001, 2101, 2311, 2401, 2631, 2654, 2801, 2901, 6201, 6331, 6341, 8211, 8411
gd$cluster <- factor(round(gd$cl.labs/10e10))
gd$cluster[gd$cluster == 1001] <- 1000
gd$cluster[gd$cluster == 2101] <- 2100
gd$cluster[gd$cluster == 2311] <- 2310
gd$cluster[gd$cluster == 2401] <- 2400
gd$cluster[gd$cluster == 2631] <- 2630
gd$cluster[gd$cluster == 2654] <- 2653
gd$cluster[gd$cluster == 2801] <- 2800
gd$cluster[gd$cluster == 2901] <- 2900
gd$cluster[gd$cluster == 6201] <- 6200
gd$cluster[gd$cluster == 6331] <- 6330
gd$cluster[gd$cluster == 6341] <- 6340
gd$cluster[gd$cluster == 8211] <- 8210
gd$cluster[gd$cluster == 8411] <- 8410


gd <- gd[, c(1:2, 8355, 3:8354)]
gd[1:10, 1:6]
#gd <- gd[, -1] #get rid of number order, no need

gd <- gd[order(gd$cluster), ]
gd[1:10, 1:6]

#write.csv(gd, file = "ILSunSh_GeneDist_HopachCor_noSPE_052014.csv")

gd <- read.csv("~/Desktop/ShadePaper/clustering/hopach/ILSunSh_GeneDist_HopachCor_noSPE_052014.csv", row.names = 1)
dim(gd)

#write.csv(gd[, c(1:3)], file = "ILSunSh_GeneDistSHORT_HopachCor_noSPE_052014.csv")

####

bobj <- boothopach(log.data.cv.sh, hobj, B = 1000)
colnames(bobj)

makeoutput(log.data.cv.sh, hobj, bobj, file = "Hopach_Clustering_Output_noSPE_051714", gene.names = rownames(log.data.cv.sh))

####


head(log.data.cv.sh)[, 1:5]
dim(log.data.cv.sh)
names(log.data.cv.sh)
gd[1:4, 1:6]

log.data.cv.sh$itag <- factor(rownames(log.data.cv.sh))
lg <- merge(gd[, c(1,3)], log.data.cv.sh, by.x = 0, by.y = 0, all.x = T, all.y = T)
head(lg)[, 1:5]

#write.csv(lg, file = "ILSunSh_LogFC_hopachCOR_noSPE_052014.csv")

lg <- read.csv("~/Desktop/ShadePaper/clustering/hopach/ILSunSh_LogFC_hopachCOR_noSPE_052014.csv", row.names = 1)
dim(lg)

lg$number <- 1:8352

#####

library(reshape2)

lgm <- melt(lg, id.vars = c("Row.names", "cluster", "cl.labs", "number"))
head(lgm)
summary(lgm)
lgm$cl.labs <- as.factor(lgm$cl.labs)
lgm$cluster <- as.factor(lgm$cluster)
lgm$value <- as.numeric(lgm$value)
lgm$Row.names <- as.factor(lgm$Row.names)

head(lgm[is.na(lgm$value), ], 50)
lgm <- lgm[!is.na(lgm$value), ]
lgm$cluster <- droplevels(lgm$cluster)
summary(lgm)
names(lgm) <- c("itag", "cluster", "cl.labs", "number", "geno", "logFC")

#write.table(lgm, file = "~/Desktop/ShadePaper/clustering/hopach/ILSunSh_LogFC_hopachCOR_noSPE_melt_052414.txt")
lgm <- read.table("~/Desktop/ShadePaper/clustering/hopach/ILSunSh_LogFC_hopachCOR_noSPE_melt_052014.txt", header = T)
dim(lgm)
lgm$cl.labs <- as.factor(lgm$cl.labs)
lgm$cluster <- as.factor(lgm$cluster)
head(lgm)


library(ggplot2)
setwd("~/Desktop/ShadePaper/clustering/hopach/")

p <- ggplot(lgm, aes(x = geno, y = logFC)) +
  geom_point(position = "jitter", size=1.5, alpha=0.6, aes(color = cluster)) + 
  geom_boxplot(outlier.size = 0, alpha = 0.8) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, size = 5)) +
  facet_grid(cluster ~ ., scales = "free", space = "free")
ggsave("BoxplotHopachCOR_noSPE_052014.pdf", p, height = 100, width = 10.5, limitsize = F)

head(lgm)
lgm$chr <- sub("IL", "", lgm$geno)
lgm$chr <- factor(sub("\\.\\d+|\\.\\d+\\.\\d+", "", lgm$chr))
lgm$chr <- factor(lgm$chr, levels(lgm$chr)[c(13, 1, 5, 6, 7, 8, 9, 10, 11, 12, 2, 3, 4)])
lgm$geno <- factor(lgm$geno, levels(lgm$geno)[c(73, 1:7, 25:72, 8:24)])
levels(lgm$chr)
levels(lgm$geno)
head(lgm)

lgo <- lgm[order(c(lgm$cluster, lgm$number), decreasing = F), ]

head(lgo)

for (i in 1:46) { # 46 for correlation dist matrix
  levs <- levels(lgm$cluster)[i]
  plotters <- lgm[lgm$cluster == levs, ]
  number <- length(unique(plotters$itag))
  
  p <- ggplot(plotters, aes(x = geno, y = logFC)) +
    geom_point(position = "jitter", size=2, alpha=0.6, aes(color = chr)) + 
    geom_boxplot(outlier.size = 0, alpha = 0.8) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90, size = 10, vjust = 0.5)) +
    facet_grid(cluster ~ ., scales = "free", space = "free") +
    labs(title = paste("There are", number, "genes in cluster", levels(lgm$cluster)[i], "!"))
  ggsave(paste0("hopach_cluster_boxplots_noSPE/BoxplotHopachCOR_cluster_", i, ".pdf"), p, height = 8, width = 10.5, limitsize = F)
}

clusterz <- NULL
for (i in 1:46) {
  bam <- unique(lgm$cluster)[i]
  clusterz$cluster[i] <- unique(lgm$cluster)[i]
  clusterz$length[i] <- length(unique(lgm$itag[lgm$cluster == bam]))
}
clusterz <- as.data.frame(clusterz)
clusterz$total[1] <- 922 
for (i in 2:46) clusterz$total[i] <- clusterz$length[i] + clusterz$total[i-1] 
for (i in 1:46) clusterz$breaks <- ceiling(clusterz$total - (0.5*clusterz$length))
head(clusterz)
tail(clusterz)

####The Hopach Clusters Just Suck, is all ####
##### plot gene expression in order of cluster by IL #####
library(ggplot2)

p <- ggplot(lgo, aes(x = geno, y = number, fill = logFC, color = logFC)) +
  geom_tile() +
  theme_bw() +
  theme(aspect.ratio = 2) + 
  scale_fill_gradientn(colours = c("green", "white", "magenta")) +
  scale_colour_gradientn(colours = c("green", "white", "magenta")) +
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
ggsave("HopachOrderClusterTest.pdf", p)


####PLOT OF GENE EXPRESSION BY CLUSTER BY IL#####
p <- ggplot(lgo, aes(x = geno, y = cluster, fill = logFC, color = logFC)) +
  geom_tile() +
  theme_bw() +
  theme(aspect.ratio = 2) + 
  scale_fill_gradientn(colours = c("green", "white", "magenta")) +
  scale_colour_gradientn(colours = c("green", "gray", "magenta")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0.8, size = 6)) +
  theme(axis.text.y = element_text(size = 7)) +
  ylab("Clustered Gene Expression") +
  xlab("Introgression Lines") +
  labs(fill = "Log2 Fold Change\n(Shade/Sun)\nExpression", 
       title = "Clustered Shade-responsive Gene Expression in each IL") +
  guides(colour = "none")
ggsave("hopachClusterTest.pdf")



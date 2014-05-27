# find SOM clusters to zoom in on 
5, 9, 12, 13, 14, 21, 23, 27, 30, 31, 49, 54, 57, 64  
# phenos of interest #30 is the one with IL9.3.2; #21 is 4.3.2
# 

# hopach
10, 30, 33/34/35, 38

# plot IL GO based on the GO enrichments in that cluster of interest
# relate to phenos
#052614

setwd("~/Desktop/ShadePaper/")


thirty <- read.csv("kohonen/KohonenClusters_noSPE/SOM_cluster30_GeneList_052014.csv", row.names = 1)
twentyone <- read.csv("kohonen/KohonenClusters_noSPE/SOM_cluster21_GeneList_052014.csv", row.names = 1)

# IL9.3.2
dim(thirty)
names(thirty)
thirty[, c(1, 73, 74, 78, 79)] #trans - there's stuff from all the chr #on ave IL9.3.2 is down reg here

# IL4.3.2
dim(twentyone)
names(twentyone)
twentyone[, c(1, 43, 74, 78)] #trans #

# what is the deal? the clusters all have a ton of genes in them now, instead of the few genes they are reported to have....
# or at least the few genes they had originally.  What changed??  Now nothing is cis-regulating. Can I figure this out?  
# initial analysis: the difference may well be clustering by IL (old way) vs. clustering by gene expression (right way).  This might lead to 
# diffs in cis vs trans plots.  Right?  

#FUCKING FUCK IS WHAT I HAVE TO SAY TO THIS, I JUST GD MADE ALL THE FUCKING PLOTS

somm <- read.table("clustering/kohonen/ILSunSh_shortDEshade_no622noSPE_LogFC_SOM8x8_melt_052614.txt", sep = " ")
head(somm)
summary(somm)
somm$cluster <- as.factor(somm$cluster)
levels(somm$chr)
levels(somm$geno)
somm$chr <- factor(somm$chr, levels(somm$chr)[c(13, 1, 5, 6, 7, 8, 9, 10, 11, 12, 2, 3, 4)])
somm$geno <- factor(somm$geno, levels(somm$geno)[c(73, 1:7, 25:72, 8:24)])

length(unique(somm$itag[somm$cluster == 12])) # plot says 106, data says 48. Hmmm.
length(unique(somm$itag[somm$cluster == 13])) # plot says 243, data says 148.  Augh!

# next work on why the plots don't look like the gene lists...
# but I think we want to plot the log.data.cv.sh data, which is relative counts, not scaled...humm...

log.data.cv.sh <- read.csv("data_tables/ILSunShDE_log2Diff_041714.csv", row.names = 1)
names(log.data.cv.sh)
log.data.cv.sh <- log.data.cv.sh[, -grep("SPE", names(log.data.cv.sh))]
head(log.data.cv.sh)[, 1:5]
dim(log.data.cv.sh)

somz <- read.table("ILSunSh_shortDEshade_no622noSPE_LogFC_SOM8x8_052614.txt", sep = " ")
dim(somz)
head(somz)
names(somz)

clust.counts <- merge(log.data.cv.sh, somz[, 74:76], by.x = 0, by.y = 0)
dim(clust.counts)
head(clust.counts)
rownames(clust.counts) <- clust.counts[, 1]
clust.counts <- clust.counts[, -1]
clust.counts$cluster <- as.factor(clust.counts$cluster)

write.csv(clust.counts, "ILSunSh_Log2Diff_ClusteredSOM8x8_052614.csv")
clust.counts <- read.csv("ILSunSh_Log2Diff_ClusteredSOM8x8_052614.csv", row.names = 1)
summary(clust.counts)
clust.counts$cluster <- as.factor(clust.counts$cluster)

library(reshape2)
ccm <- melt(clust.counts, id.vars = c("itag", "number", "cluster"))
head(ccm)

levels(ccm$variable)
names(ccm)[4] <- "geno"
levels(ccm$geno)
ccm$geno <- factor(ccm$geno, levels = levels(ccm$geno)[c(73, 1:7, 25:72, 8:24)])
ccm$geno <- relevel(ccm$geno, ref = "SLY")

head(ccm)
names(ccm)[5] <- "logFC"
summary(ccm)
#ccm$cluster <- as.factor(ccm$cluster)
head(ccm)

ccm$chr <- sub("IL", "", ccm$geno)
ccm$chr <- factor(sub("\\.\\d+|\\.\\d+\\.\\d+", "", ccm$chr))
ccm$chr <- factor(ccm$chr, levels(ccm$chr)[c(13, 1, 5, 6, 7, 8, 9, 10, 11, 12, 2, 3, 4)])
ccm$geno <- factor(ccm$geno, levels(ccm$geno)[c(73, 1:7, 25:72, 8:24)])
levels(ccm$chr)
levels(ccm$geno)

write.csv(ccm, "ILSunSh_Log2Diff_ClusteredSOM8x8_melt_052614.csv")

ccm <- read.csv("ILSunSh_Log2Diff_ClusteredSOM8x8_melt_052614.csv", row.names = 1)
ccm$chr <- factor(ccm$chr, levels(ccm$chr)[c(13, 1, 5, 6, 7, 8, 9, 10, 11, 12, 2, 3, 4)])
ccm$geno <- factor(ccm$geno, levels(ccm$geno)[c(73, 1:7, 25:72, 8:24)])
summary(ccm)
ccm$cluster <- as.factor(ccm$cluster)

library(ggplot2)

length(unique(ccm$itag[ccm$cluster == 1])) #222
length(clust.counts$itag[clust.counts$cluster == 1]) #222

for (i in 41:64) {
  plotters <- ccm[ccm$cluster == i, ]
  #dim(plotters)
  #head(plotters)
  #summary(plotters)
  number <- length(unique(plotters$itag))
  
  p <- ggplot(plotters, aes(x = geno, y = logFC)) +
    geom_point(position = "jitter", size=2, alpha=0.6, aes(color = chr)) + 
    geom_boxplot(outlier.size = 0, alpha = 0.8) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90, size = 10, vjust = 0.5)) +
    facet_grid(cluster ~ ., scales = "free", space = "free") +
    labs(title = paste("There are", number, "genes in this cluster!"))
  ggsave(paste0("clustering/kohonen/8x8cluster_boxplots_noSPE/BoxplotSOM8x8_cluster_", i, ".pdf"), p, height = 8, width = 10.5, limitsize = F)
}

# look at these plots: awesome.  Still see the trends, now the genes match.  Re-write the gene lists to be sure.
# thank goodness!

for (i in 1:64) {
  writers <- clust.counts[clust.counts$cluster == i, ]
  dim(writers) 
  tail(writers)
  filename <- paste0("clustering/kohonen/KohonenClusters_noSPE/SOM_cluster", i, "_GeneList_052614.csv")
  write.csv(writers, filename)
}

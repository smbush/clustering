# self-organizing maps!  using Kohonen package
# April 2014
# start with IL*trt logFC genes (>17700)
# use the log(ILsha/ILsun) value (aka ILsha + M82sha - ILsun)
# reduce gene number by getting rid of genes with low variance across ILs (unlog value first)
# reduce gene number by including 8770 genes DE in trtsha (no geno factor)
# exclude PENNELLII!!!!

setwd("~/Desktop/ShadePaper/")

#install.packages("kohonen")

library("kohonen")

data <- read.csv("data_tables/ILSunShRenamed_adjDE_032914.csv", row.names = 1)

names(data)
dim(data)
head(data)[, 1:5]
head(data)[, c(1:3, 75:78)]
rownames(data) <- data[, 1]
data <- data[, -1]

log.data <- data
log.data[, 1:75] <- log.data[77:151] - log.data[, 1:75]
head(data)[, 1:5]
head(log.data)[, 1:5]
log.data <- log.data[, 1:76]
head(log.data)[, 1:5]
head(log.data)[, 72:76]

names(log.data) <- sub("logFC.geno", "", names(log.data))
names(log.data)[76] <- "SLY"
head(log.data)[, 1:5]
log.data <- log.data[, -grep("IL6.2", names(log.data))]

unlog.data <- log.data
unlog.data <- apply(unlog.data, c(1,2), function (x) 2^(x))
unlog.data <- as.data.frame(unlog.data)
head(log.data)[, 1:5]
head(unlog.data)[, 1:5]

unlog.data$CV <- apply(unlog.data, 1, function(x) sd(x)/mean(x))
head(unlog.data)

hist(apply(unlog.data[,1:74],1,mean))
hist(apply(unlog.data[,1:74],1,sd), breaks=150)
hist(apply(unlog.data[,1:74],1,function(x) log2(sd(x))), breaks = 150)

summary(abs(unlog.data$CV))
unlog.data.short <- unlog.data[abs(unlog.data$CV) >= 0.1, ] # this knocks off about 1000 genes
dim(unlog.data.short)
dim(unlog.data)

log.data.cv <- log.data[rownames(log.data) %in% row.names(unlog.data.short), ]
dim(log.data.cv)
summary(log.data.cv)

shade <- read.csv("~/Documents/ILSunShadePaper/edgeR_DEoutputs/ILSunShRenamed_shaDE_010414.csv", row.names = 1)
dim(shade) #8770 genes
head(shade)[, 1:5] 

log.data.cv.sh <- log.data.cv[rownames(log.data) %in% shade$Row.names, ]
dim(log.data.cv.sh)
summary(log.data.cv.sh)
head(log.data.cv.sh)[, 1:5]
log.data.cv.sh <- log.data.cv.sh[!is.na(log.data.cv.sh$IL1.1), ]
summary(log.data.cv.sh[, 1:5])

write.csv(log.data.cv.sh, file = "ILSunShDE_log2Diff_041714.csv")

log.data.cv.sh <- read.csv("data_tables/ILSunShDE_log2Diff_041714.csv", row.names = 1)

# transpose the data bc you want to SCALE BY GENE NOT BY IL
data.sc <- t(scale(t(log.data.cv.sh[, names(log.data.cv.sh) != "SPE"]), center = F))

# cluster rows
hcGenes <- hclust(dist(data.sc))
#plot(hcGenes)

# transpose the matrix and cluster columns
hcILs <- hclust(dist(t(data.sc)))
pdf(paste0("GenexIL_", namelist[i], "_ILsCluster.pdf"))
plot(hcILs, cex = 0.6)
dev.off()

library(kohonen)
data.som <- som(data = data.sc, grid = somgrid(20,10,"hexagonal"))
data.som <- som(data = data.sc, grid = somgrid(10,10,"hexagonal"))
data.som <- som(data = data.sc, grid = somgrid(8, 8, "hexagonal")) # I like this one arbitrarily
data.som <- som(data = data.sc, grid = somgrid(5,5,"hexagonal"))

setwd("clustering/kohonen/")

pdf("SOM8x8_052614.pdf")
plot(data.som, type="changes")
plot(data.som, type="codes")
plot(data.som, type="counts")
plot(data.som, type="quality")
dev.off()

cutree(hclust(dist(data.som$codes)), 20)
cutree(hclust(dist(data.som$codes)), 50)

somz <- as.data.frame(cbind(data.som$data, data.som$unit.classif))
dim(somz)
names(somz)[74] <- "cluster"
somz$itag <- factor(substring(rownames(somz), 1, 14))
head(somz)


library(reshape2)

somz$number <- 1:8352
somm <- melt(somz, id.vars = c("itag", "cluster", "number"))
head(somm)
tail(somm)
dim(somm)

levels(somm$variable)
names(somm)[4] <- "geno"
levels(somm$geno)
somm$geno <- factor(somm$geno, levels = levels(somm$geno)[c(1:7, 25:73, 8:24)])
somm$geno <- relevel(somm$geno, ref = "SLY")

head(somm)
names(somm)[5] <- "logFC"
summary(somm)
somm$cluster <- as.factor(somm$cluster)
head(somm)

write.table(somz, file="ILSunSh_shortDEshade_no622noSPE_LogFC_SOM8x8_052614.txt")
write.table(somm, file = "ILSunSh_shortDEshade_no622noSPE_LogFC_SOM8x8_melt_052414.txt")

somm <- read.table("clustering/kohonen/ILSunSh_shortDEshade_no622noSPE_LogFC_SOM8x8_melt_052414.txt", sep = " ")
head(somm)

library(ggplot2)

p <- ggplot(somm, aes(x = geno, y = logFC)) +
  geom_point(position = "jitter", size=1.5, alpha=0.6, aes(color = cluster)) + 
  geom_boxplot(outlier.size = 0, alpha = 0.8) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, size = 5)) +
  facet_grid(cluster ~ ., scales = "free", space = "free")
ggsave("BoxplotSOM8x8_noSPE_052014.pdf", p, height = 100, width = 10.5, limitsize = F)

head(somm)
somm$chr <- sub("IL", "", somm$geno)
somm$chr <- factor(sub("\\.\\d+|\\.\\d+\\.\\d+", "", somm$chr))
somm$chr <- factor(somm$chr, levels(somm$chr)[c(13, 1, 5, 6, 7, 8, 9, 10, 11, 12, 2, 3, 4)])
somm$geno <- factor(somm$geno, levels(somm$geno)[c(73, 1:7, 25:72, 8:24)])
levels(somm$chr)
levels(somm$geno)

write.table(somm, file = "clustering/kohonen/ILSunSh_shortDEshade_no622noSPE_LogFC_SOM8x8_melt_052614.txt")


for (i in 1:64) {
  plotters <- somm[somm$cluster == i, ]
  number <- length(unique(plotters$itag))
  
  p <- ggplot(plotters, aes(x = geno, y = logFC)) +
    geom_point(position = "jitter", size=2, alpha=0.6, aes(color = chr)) + 
    geom_boxplot(outlier.size = 0, alpha = 0.8) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90, size = 10, vjust = 0.5)) +
    facet_grid(cluster ~ ., scales = "free", space = "free") +
    labs(title = paste("There are", number, "genes in this cluster!"))
  ggsave(paste0("8x8cluster_boxplots_noSPE/BoxplotSOM8x8_cluster_", i, ".pdf"), p, height = 8, width = 10.5, limitsize = F)
}

somo <- somm[order(c(somm$cluster, somm$number), decreasing = F), ]
head(somo)
tail(somo)
summary(somo)
somo <- na.omit(somo)
head(somo)
tail(somo)
dim(somo)

write.csv(somo, "clustering/kohonen/ILSunSh_LogFC_Som8x8_melt_052414.csv")

####plot of gene expression for each gene, listed by cluster####
clusterz <- NULL
for (i in 1:64) {
clusterz$cluster[i] <- i
clusterz$length[i] <- length(unique(somm$itag[somm$cluster == i]))
}
clusterz <- as.data.frame(clusterz)
clusterz$total[1] <- 100 
for (i in 2:64) clusterz$total[i] <- clusterz$length[i] + clusterz$total[i-1] 
for (i in 1:64) clusterz$breaks <- round(clusterz$total - (0.5*clusterz$length))
head(clusterz)
tail(clusterz)

write.csv(clusterz, "ILSunSh_KohonenClusterAxisBreaks.csv")

library(ggplot2) 

##### plot gene expression in order of cluster by IL #####
  p <- ggplot(somo, aes(x = geno, y = number, fill = logFC, color = logFC)) +
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
  ggsave("orderClusterTest.pdf")


####PLOT OF GENE EXPRESSION BY CLUSTER BY IL#####
p <- ggplot(somo, aes(x = geno, y = cluster, fill = logFC, color = logFC)) +
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
ggsave("clusterTest.pdf")

# GO analysis on SOM clusters
# 041514
setwd("~/Desktop/ShadePaper/clustering/")

somm <- read.table("kohonen/ILSunSh_shortDEshade_no622noSPE_LogFC_SOM8x8_052614.txt", header = T)
#somm <- read.csv("hopach/ILSunSh_GeneDistSHORT_HopachCor_noSPE_052014.csv", row.names = 1)
dim(somm)
head(somm)

somm$cluster <- as.factor(somm$cluster)

annot1 <- read.delim("~/Documents/Daniele'sdata_script/ITAG2.3_all_Arabidopsis_annotated.tsv", header = T)[, 1:4]
dim(annot1)
head(annot1)
annot2 <- read.delim("~/Documents/ILSunShadePaper/ITAG2.3_HRD.tsv", header = T)
dim(annot2)
head(annot2)
annot <- merge(annot1, annot2, by.x = 1, by.y = 1, all.x = T, all.y = T)
names(annot)[5] <- "HRD"

soma <- merge(somm, annot, by.x = 0, by.y = 1, all.x = T, all.y = F)
dim(soma)
soma$cluster <- as.factor(soma$cluster)
soma$Row.names <- as.factor(soma$Row.names)
levels(soma$cluster)
dim(soma)

#setwd("~/Desktop/ShadePaper/clustering/hopach/Hopach_genelist_noSPE/")
setwd("~/Desktop/ShadePaper/clustering/kohonen/KohonenClusters_noSPE/")
for (i in 1:64) {
  new <- subset(soma, cluster == unique(soma$cluster)[i])
  filename <- paste0("SOM8x8_cluster", i, "_GeneList_052014.csv")
  write.csv(new, file = filename)
}


library(goseq)
library(GO.db)

  #transcripts file
  Bra_trans_len <- read.table("~/Documents/gofun/ITAG2.3_cds.mod.length", sep="\t", h=T)
  
  colnames(Bra_trans_len) <- c("Gene", "length")
  Bra_trans_lenf <- Bra_trans_len[, 1]
  Bra_trans_len2 <- Bra_trans_len[, 2]
  genes.in.annot <- Bra_trans_lenf
  
  #GO file
  go <- read.table("~/Documents/gofun/GO.table", h = T)
  colnames(go) <- c("Gene", "GO")
  go.list <- strsplit(as.character(go[, 2]),split = ",", fixed = T)
  names(go.list) <- as.character(go[, 1])
  
  #expression file


for (i in 1:64) { #46 for hopachCOR #64 for kohonen
  genes <- somm[somm$cluster == unique(somm$cluster)[i], ]
  genes$Row.names <- factor(rownames(genes))
  if (length(genes$Row.names) > 10) {
    print(paste0("cluster ", unique(somm$cluster)[i], " has ", length(genes$Row.names), " genes!"))
    
    genes.DEup <- as.numeric(genes$Row.names != 1) 
    names(genes.DEup) <- genes$Row.names
    
  genes.DEup.zero.names <- Bra_trans_lenf[!Bra_trans_lenf %in% genes$Row.names]
  genes.DEup.zero <- rep(0, length(genes.DEup.zero.names))
  names(genes.DEup.zero) <- genes.DEup.zero.names
  
  genes.DEup.all <- c(genes.DEup, genes.DEup.zero)
  
  genes.DEup.all <- genes.DEup.all[names(genes.DEup.all) %in% names(go.list)]
  
    bias <- Bra_trans_len2
    names(bias) <- genes.in.annot
    UP_bias <- bias[names(bias) %in% names(genes.DEup.all)]
    
    pwfup <- nullp(genes.DEup.all, bias.data = UP_bias) 
    go.analysis <- goseq(pwfup, gene2cat = go.list)
    
    print(length(go.analysis$category[p.adjust(go.analysis$over_represented_pvalue, method="BH") < 0.01]))
    enriched.GO_category <- go.analysis$category[go.analysis$over_represented_pvalue < 0.01]
    enriched.GO_p_value <- go.analysis$over_represented_pvalue[go.analysis$over_represented_pvalue < 0.01]
    enriched.GO <- as.data.frame(cbind(enriched.GO_category, enriched.GO_p_value))
    enriched.GO$term <- Term(as.character(enriched.GO$enriched.GO_category))
    enriched.GO$ont <- Ontology(as.character(enriched.GO$enriched.GO_category))
    
    
    if (length(enriched.GO$term) > 0) {
     # print(head(enriched.GO.up))
      filename <- paste0("kohonen/Kohonen_GOterms_noSPE/SOM8x8_cluster", i, "_GOenrich_052614.csv")
      write.csv(enriched.GO, file = filename)
    
    }
  }
}


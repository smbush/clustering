# NEXT UP!!

# find SOM clusters to zoom in on 
28, 37, 43, 48, 52, 60 #52 is the one with IL9.3.2; #48 is 4.3.2

# first calculate the values of each GO term in each IL
# alternatively for each GO term of interest, calculate their values in ILs

GO52 <- read.csv("~/Desktop/ShadePaper/clustering/kohonen/Kohonen_GOterms_noSPE/SOM8x8_cluster52_GOenrich_052614.csv", row.names = 1)
GO48 <- read.csv("~/Desktop/ShadePaper/clustering/kohonen/Kohonen_GOterms_noSPE/SOM8x8_cluster48_GOenrich_052614.csv", row.names = 1)

GO52$expt <- "SOM52"
GO48$expt <- "SOM48"

GO <- rbind(GO52, GO48)
head(GO)
names(GO) <- c("GOcat", "p.value", "term", "ont", "som")

genotype <- colnames(somm)[1:72]
genotype <- genotype[c(1:7, 25:72, 8:24)]

setwd("~/Desktop/ILSunShData/Genotype_Data/ILspecDEgenes/")
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

for (g in genotype) {
  filename <- paste0("ILSunSh_", g, "_052014_AllGenes.csv")
  genes <- read.csv(filename, row.names = 1)
  
  genes.DEup <- as.numeric(genes$Row.names[genes$logFC < 0.05]) 
  names(genes.DEup) <- genes$Row.names[genes$logFC < 0.05]
   
  genes.DEup <- genes.DEup[names(genes.DEup) %in% names(go.list)]
  
  bias <- Bra_trans_len2
  names(bias) <- genes.in.annot
  UP_bias <- bias[names(bias) %in% names(genes.DEup)]
  
  pwfup <- nullp(genes.DEup, bias.data = UP_bias) 
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



# plot IL GO based on the GO enrichments in that cluster of interest
# relate to phenos

# determine which SOM clusters the flowering or cell wall etc. genes are in
# zoom in on that/those clusters
# plot IL GO based on the GO enrichments in that cluster of interest
# relate to phenos

# find Tomato genes/probes in shade - finalized list
# 


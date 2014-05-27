## PCA
## 2014 1 15 Yasu
##  modified by Kazu 
##   check plus and minus in data (plus should be shade-induced genes) (021014)
##	 use non centered, but still scalled data for pcrcomp to preserve plus/minus of fold change
#setwd( "D:/ResearchProject(Yasu)/Projects/Others/sun-shade")
setwd("/Volumes/Data5/NGS_related/Arabidopsis_analysis/SAS_muts_time_course_RNAseq")
##library
library(ggplot2)
library(class)
library(MASS)
library(kohonen)
library(reshape)
library(plyr)
load("~/Downloads/DGEsummary.response.Rdata")
## read in data
data <- subset(DGEsummary.response[2:46])
data <- data*(-1) # plus is shade-induced, munus is shade-reduced
gene <- subset(DGEsummary.response[1])
rownames(data) <- gene[,1]
data[is.na(data)] <- 0

# raw data (expression level)
load("summary.Rdata")
head(summary[rownames(summary) %in% rownames(data),])
head(data)
## plus = shade induced, minus = shade repressed (021014)
## rearrange dataset
names(data)

Col <- subset(data[1:5])
colnames(Col) <- c("1h","4h","16h","25h","49h")
Col$genotype <- rep("Col",length(Col[,1]))

yucQ <- subset(data[6:10])
colnames(yucQ) <- c("1h","4h","16h","25h","49h")
yucQ$genotype <- rep("YuQ",length(yucQ[,1]))

coi1 <- subset(data[11:15])
colnames(coi1) <- c("1h","4h","16h","25h","49h")
coi1$genotype <- rep("Coi",length(coi1[,1]))

jar1 <- subset(data[16:20])
colnames(jar1) <- c("1h","4h","16h","25h","49h")
jar1$genotype <- rep("Jar",length(jar1[,1]))

pif3 <- subset(data[21:25])
colnames(pif3) <- c("1h","4h","16h","25h","49h")
pif3$genotype <- rep("Pi3",length(pif3[,1]))

pif45 <- subset(data[26:30])
colnames(pif45) <- c("1h","4h","16h","25h","49h")
pif45$genotype <- rep("Pi4",length(pif45[,1]))

spt <- subset(data[31:35])
colnames(spt) <- c("1h","4h","16h","25h","49h")
spt$genotype <- rep("Spt",length(spt[,1]))

phyB <- subset(data[36:40])
colnames(phyB) <- c("1h","4h","16h","25h","49h")
phyB$genotype <- rep("PhB",length(phyB[,1]))

co <- subset(data[41:45])
colnames(co) <- c("1h","4h","16h","25h","49h")
co$genotype <- rep("Co_",length(co[,1]))

data2 <- rbind(Col,yucQ,coi1,jar1,pif3,pif45,spt,phyB,co)
data2.s <- subset(data2[1:5])
rownames(data2.s)<-paste(rownames(data2),data2$genotype,sep=".")

## remove genes with sd=0 for PCA
data2.s$sd <- apply(data2.s, 1, function(d) sd(d)) 
data3 <- subset(data2.s, sd>0)
data3 <- subset(data3[1:5]) #8 genes were removed 

#> rownames(data2.s[which(data2.s$sd==0),])
#[1] "AT3G30720.11.YuQ" "AT3G46370.11.YuQ" "AT5G42825.11.YuQ" "AT1G61120.12.Coi" "AT1G53480.15.Pi4" "AT5G24860.18.Co_"
#[7] "AT1G69120.18.Co_" "AT3G54340.18.Co_"

## Create a matrix of the data to perform a PCA on and scale it
sc.sub.data1 <- t(scale(t(data3),center=T)) # scaling with centering

sc.sub.data2 <- t(scale(t(data3),center=F)) # scaling , no centering
sc.sub.data3<-t(scale(t(data3), center = FALSE, scale = apply(t(data3), 2, sd, na.rm = TRUE))) # from ?scale "To scale by the standard deviations without centering, use"
# what this scale does?
head(apply(data3,1,var)) # origianl data
# AT1G01040.2.Col AT1G01060.1.Col AT1G01070.1.Col AT1G01120.1.Col AT1G01170.1.Col AT1G01190.1.Col 
     # 0.01345585      0.28097441      0.06656659      0.10483116      0.02799180      0.08142117 

head(apply(sc.sub.data1,1,var)) # scaling with centering
# AT1G01040.2.Col AT1G01060.1.Col AT1G01070.1.Col AT1G01120.1.Col AT1G01170.1.Col AT1G01190.1.Col 
              # 1               1               1               1               1               1 

head(apply(sc.sub.data2,1,var)) # scaling without centering
# AT1G01040.2.Col AT1G01060.1.Col AT1G01070.1.Col AT1G01120.1.Col AT1G01170.1.Col AT1G01190.1.Col 
      # 0.9994313       0.6043057       0.5712850       0.8864846       0.3430782       0.1565278 

head(apply(sc.sub.data3,1,var)) # scaling with SD without centering
AT1G01040.2.Col AT1G01060.1.Col AT1G01070.1.Col AT1G01120.1.Col AT1G01170.1.Col AT1G01190.1.Col 
              1               1               1               1               1               1 
# OK!
# how about mean? centered? yes (see ?scale)
head(apply(sc.sub.data,1,mean))
# AT1G01040.2.Col AT1G01060.1.Col AT1G01070.1.Col AT1G01120.1.Col AT1G01170.1.Col AT1G01190.1.Col 
  # -1.112391e-17   -3.890117e-17    0.000000e+00    1.665335e-17   -6.664862e-17    3.890117e-17 
head(apply(sc.sub.data2,1,mean))
# AT1G01040.2.Col AT1G01060.1.Col AT1G01070.1.Col AT1G01120.1.Col AT1G01170.1.Col AT1G01190.1.Col 
    # -0.02133044      0.56263260      0.58563809      0.30135082     -0.72493961      0.82144858 
head(apply(sc.sub.data3,1,mean))
# AT1G01040.2.Col AT1G01060.1.Col AT1G01070.1.Col AT1G01120.1.Col AT1G01170.1.Col AT1G01190.1.Col 
    # -0.02133651      0.72376329      0.77482369      0.32006397     -1.23767116      2.07627395 

head(data3)
                         # 1h          4h         16h         25h         49h
# AT1G01040.2.Col -0.13185000  0.11881699 -0.12124149  0.07150712  0.05039227
# AT1G01060.1.Col -0.19060579  0.81892926 -0.15399736  0.51374247  0.93015816
# AT1G01070.1.Col  0.43862547 -0.17005639  0.03935891  0.39084943  0.30076460
# AT1G01120.1.Col -0.16066576  0.18229686  0.63566130 -0.08024666 -0.05889994
# AT1G01170.1.Col -0.07797153 -0.03288463 -0.24436071 -0.45879575 -0.22134575
# AT1G01190.1.Col  0.54284142  0.28448248  0.38633116  0.75808673  0.99051932

head(sc.sub.data1) #ce
                        # 1h         4h        16h        25h         49h
# AT1G01040.2.Col -1.1153076  1.0456264 -1.0238544  0.6377806  0.45575494
# AT1G01060.1.Col -1.0833492  0.8211816 -1.0142859  0.2454338  1.03101970
# AT1G01070.1.Col  0.9252419 -1.4339441 -0.6222727  0.7400671  0.39090786
# AT1G01120.1.Col -0.8162884  0.2429692  1.6432085 -0.5679099 -0.50197944
# AT1G01170.1.Col  0.7716339  1.0411191 -0.2228773 -1.5045594 -0.08531638
# AT1G01190.1.Col -0.1738632 -1.0792930 -0.7223600  0.5804732  1.39504304

head(sc.sub.data2)
                        # 1h         4h        16h        25h        49h
# AT1G01040.2.Col -1.1363208  1.0239986 -1.0448936  0.6162687  0.4342949
# AT1G01060.1.Col -0.2795317  1.2009954 -0.2258438  0.7534257  1.3641174
# AT1G01070.1.Col  1.2849674 -0.4981856  0.1153032  1.1450060  0.8810995
# AT1G01120.1.Col -0.4672117  0.5301144  1.8484860 -0.2333552 -0.1712794
# AT1G01170.1.Col -0.2729714 -0.1151262 -0.8554852 -1.6062033 -0.7749119
# AT1G01190.1.Col  0.7526621  0.3944415  0.5356570  1.0511046  1.3733777

head(sc.sub.data3)
                        # 1h         4h        16h        25h        49h
# AT1G01040.2.Col -1.1366441  1.0242899 -1.0451909  0.6164441  0.4344184
# AT1G01060.1.Col -0.3595859  1.5449449 -0.2905226  0.9691970  1.7547830
# AT1G01070.1.Col  1.7000656 -0.6591204  0.1525510  1.5148908  1.1657316
# AT1G01120.1.Col -0.4962244  0.5630332  1.9632725 -0.2478459 -0.1819155
# AT1G01170.1.Col -0.4660372 -0.1965520 -1.4605484 -2.7422305 -1.3229875
# AT1G01190.1.Col  1.9024108  0.9969809  1.3539139  2.6567471  3.4713170


tisdata1 <- as.matrix(sc.sub.data1, dimnames=list(rownames(data2.s)) ) # original version
tisdata2 <- as.matrix(sc.sub.data2, dimnames=list(rownames(data2.s)) )
tisdata3 <- as.matrix(sc.sub.data3, dimnames=list(rownames(data2.s)) ) # 

## Perform the PCA
tispca1.1 <- prcomp(tisdata1, scale=TRUE) # scaling...
tispca1.2 <- prcomp(tisdata1, scale=FALSE) # probalyb tispca and tispca2 is the same
tispca1.3 <- prcomp(tisdata1, scale=FALSE,center=FALSE) # no scaling

summary(tispca1.1) # original
#Importance of components:
#  PC1    PC2    PC3    PC4      PC5
#Standard deviation     1.1701 1.1634 1.0842 1.0497 1.94e-16
#Proportion of Variance 0.2738 0.2707 0.2351 0.2204 0.00e+00
#Cumulative Proportion  0.2738 0.5445 0.7796 1.0000 1.00e+00

summary(tispca1.2)
# Importance of components:
                          # PC1    PC2    PC3    PC4       PC5
# Standard deviation     1.1021 1.0395 0.9282 0.9097 3.729e-15
# Proportion of Variance 0.3048 0.2712 0.2162 0.2077 0.000e+00
# Cumulative Proportion  0.3048 0.5760 0.7923 1.0000 1.000e+00
 
summary(tispca1.3)
# Importance of components:
                          # PC1    PC2    PC3   PC4       PC5
# Standard deviation     1.1073 1.0409 0.9288 0.910 3.556e-15
# Proportion of Variance 0.3065 0.2708 0.2157 0.207 0.000e+00
# Cumulative Proportion  0.3065 0.5774 0.7930 1.000 1.000e+00

tispca3.1 <- prcomp(tisdata3, scale=TRUE) # scaling...
tispca3.2 <- prcomp(tisdata3, scale=FALSE) # probalyb tispca and tispca2 is the same
tispca3.3 <- prcomp(tisdata3, scale=FALSE,center=FALSE) # no scaling, no centering

summary(tispca3.1)
summary(tispca3.3) # use this


## Retrieve PCA scores
tis.pca.scores <- data.frame(tispca3.3$x)
names(tis.pca.scores)
rownames(tis.pca.scores)
head(tis.pca.scores)

## Write out master data files with original data, scaled data, and PCA results
#data.val <- cbind(data3,sc.sub.data1,tis.pca.scores) # original
data.valz <- cbind(data3,sc.sub.data3,tis.pca.scores)

# write.table(data.val, file="pca.scores.transcripts.all.txt")
# write.table(tispca$rotation, "loadings.transcripts.all.txt")

tiff("plain.transcripts.all.PC1PC2.tif", width=8, height=8, unit="in",compression="lzw",res=100)
t <- ggplot(data.val, aes(PC1, PC2))
t + geom_point(alpha=0.10) + theme_bw() 
dev.off()

tiff("plain.transcripts.all.PC3PC4.tif", width=8, height=8, unit="in",compression="lzw",res=100)
t <- ggplot(data.val, aes(PC3, PC4))
t + geom_point(alpha=0.10) + theme_bw() 
dev.off()

## SOM
set.seed(2) # Set a random seed so that SOM results are reproducible
ssom <- som(sc.sub.data3, somgrid(3,3,"hexagonal")) # try (2,5), (3,4) and look plot(ssom, type = "codes") below to see if som can detect stage specific clustering (not too much details) (by Yasu, 021014)
plot(ssom, type = "codes") # 

ssom2.5 <- som(sc.sub.data3, somgrid(2,5,"hexagonal")) # try (2,5), (3,4) and look plot(ssom, type = "codes") below to see if som can detect stage specific clustering (not too much details) (by Yasu, 021014)
quartz();plot(ssom2.5) 

ssom3.4 <- som(sc.sub.data3, somgrid(3,4,"hexagonal")) # try (2,5), (3,4) and look plot(ssom, type = "codes") below to see if som can detect stage specific clustering (not too much details) (by Yasu, 021014)
quartz();plot(ssom3.4) # I like five peaks and five trough and one big (no peak and no trough), one small (no chagnes) 

summary(ssom3.4)
#som map of size 3x3 with a hexagonal topology.
#Training data included; dimension is 21070 by 5
#Mean distance to the closest unit in the map: 1.370107

# ch <- ssom$changes
# ch <- as.data.frame(ch)
# ch$x <- c(1:100)
# tiff("ssom.changes.transcripts.all.3x3.tif", width=8, height=8, unit="in",compression="lzw",res=100)
# ggplot(ch,aes(x,V1))+geom_line(size=1.2,colour="grey30")+theme_bw()
# dev.off()

ch <- ssom3.4$changes
ch <- as.data.frame(ch)
ch$x <- c(1:100)
tiff("ssom3.4.changes.transcripts.all.3x3.tif", width=8, height=8, unit="in",compression="lzw",res=100)
ggplot(ch,aes(x,V1))+geom_line(size=1.2,colour="grey30")+theme_bw()
dev.off()

# plot(ssom, type ="changes")
# plot(ssom, type = "codes") # 
# plot(ssom, type = "counts")
# plot(ssom, type = "quality")

plot(ssom3.4, type ="changes")
plot(ssom3.4, type = "codes") # 
plot(ssom3.4, type = "counts")
plot(ssom3.4, type = "quality")

## Create and write-out master SOM file
#data.val2 <- cbind(data.val,ssom$unit.classif,ssom$distances)
data.val3.4 <- cbind(data.valz,ssom3.4$unit.classif,ssom3.4$distances)

head(data.val3.4)
write.csv(data.val3.4, file="supersom.data.transcripts.all.3x4.csv")

## Codes for the SOM nodes
codes <- ssom3.4$codes
head(codes)
write.table(codes, file="codes.transcripts.all.3x4.txt")

# ## PC graphs
# tiff("node.transcripts.all.PC1PC2.3x3.tif", width=10, height=8, unit="in",compression="lzw",res=100)
# t <- ggplot(data.val2, aes(PC1, PC2))
# t + geom_point(alpha=0.3, size=2.5,aes(colour=factor(ssom$unit.classif))) + theme_bw() + 
  # scale_colour_manual(values=c("tomato","springgreen3","blue3","orange1","magenta3","yellow4","turquoise2","deeppink1","black"))
# dev.off()

# tiff("node.transcripts.all.PC3PC4.3x3.tif", width=10, height=8, unit="in",compression="lzw",res=100)
# t <- ggplot(data.val2, aes(PC3, PC4))
# t + geom_point(alpha=0.3,size=2.5, aes(colour=factor(ssom$unit.classif))) + theme_bw() + 
  # scale_colour_manual(values=c("tomato","springgreen3","blue3","orange1","magenta3","yellow4","turquoise2","deeppink1","black","palegreen4","yellow","green"))
# dev.off()

t<-ggplot(data.val3.4, aes(PC1, PC2))
t<- t + geom_point(alpha=0.3,size=1.5, aes(colour=factor(ssom3.4$unit.classif))) + theme_bw() + scale_x_continuous(limits=c(-5,5)) +
 scale_colour_manual(values=c("tomato","springgreen3","blue3","orange1","magenta3","yellow4","turquoise2","deeppink1","black","palegreen4","yellow","green"))
t
ggsave("node.transcripts.all.PC1PC2.3x4.pdf",plot=t,width=11,height=8)
# if you start from data.val2
data.val3.4<-read.csv("supersom.data.transcripts.all.3x4.csv",row.names=1) # this is supposed to be (original data.val2 in Yasu's script)

## Boxplots for each node (Kazu added genotype part)
gt <- 
  unique(gsub("(AT[[:alnum:]]+)(\\.)([[:digit:]]+)(\\.)([[:print:]]+)", "\\5", rownames(data.val3.4)))
gt # [1] "Col" "YuQ" "Coi" "Jar" "Pi3" "Pi4" "Spt" "PhB" "Co_"
# option 1
for(genotype in gt) {
for(i in c(1:9)){
  cluster=i
  sub.group <- subset(data.val2, ssom.unit.classif==cluster)
  sub.group.gt<-sub.group[grep(genotype,rownames(sub.group),]
  expression <- sub.group.gt[,6:10]
  expression <- cbind(expression,rownames(expression))
  m.expression <- melt(expression, id="rownames(expression)") # multiple column into one column
  print(paste("genotype is ",genotype))
    p <- ggplot(m.expression, aes(x=variable, y=value))
	p <- p + geom_point(position="jitter",size=1.5,alpha=0.6) + geom_boxplot(outlier.size=0, alpha=0.8) + theme_bw() +labs(title=paste("boxplot.node",i,"in",genotype))
  ggsave(file=paste("boxplot.node",i,"in",genotype,".tiff"), width=4, height=8, unit="in",compression="lzw",dpi=100)  
}
}
# option2
data.val3.4$gt<-as.factor(gsub("(AT[[:alnum:]]+)(\\.)([[:digit:]]+)(\\.)([[:print:]]+)","\\5",rownames(data.val3.4)))
data.val3.4$AGI<-as.factor(gsub("(AT[[:alnum:]]+)(\\.)([[:digit:]]+)(\\.)([[:print:]]+)","\\1",rownames(data.val3.4)))
names(data.val3.4)[6:10]<-c("1","4","16","25","49")
summary(data.val3.4)

data.val3.4melt<-melt(data.val3.4[,c(6:10,16,18:19)],id.var=c("gt","AGI","ssom3.4.unit.classif"))
names(data.val3.4melt)[3]<-"ssom3.4unit.classif"
data.val3.4melt$ssom3.4unit.classif<-factor(data.val3.4melt$ssom3.4unit.classif,levels=c(10,7,6,11,3,9,5,1,4,8,12,2)) # use this factor for all heatmap
p<-ggplot(data.val3.4melt,aes(x=variable,y=value,color=value))
p<-p + geom_point(position="jitter",size=1.5,alpha=0.3) + geom_boxplot(outlier.size=0,aes(fill=factor(ssom3.4unit.classif),alpha=0.8)) + theme_bw() +xlab("sun/shade treatment (hr)") + scale_colour_gradient2(low="green1",mid="gray",high="magenta")
p<-p + facet_grid(ssom3.4unit.classif~gt,scale="free",space="free") #+ facet_wrap(ssom3.4unit.classif~.)
#p<- p + theme(strip.background = element_rect(colour =factor(data.val3.4melt$ssom3.4unit.classif))) # does not work as I wanted (022014)
#p
ggsave("boxplot.node.genotype.pdf",width=11,height=20,unit="in")
ggsave("boxplot.node.genotype.png",width=11,height=20,unit="in",dpi=600)

###
#
# GO analysis (I did not use this)
#
##
library(org.At.tair.db)

# Upendra's function
library(org.At.tair.db)

annotate <- function(result,sample) { #results is a toptags object
  results.annotated <- result$table #need to extract the gene table from results to get started
  results.annotated$AGI <- substr(rownames(results.annotated),1,9) #get rid of the gene model number
  
  #add At GO
  Atgo <- toTable(org.At.tairGO)
  head(Atgo)
  
  BP <- TRUE #only keep BP go TERMS
  if (BP) Atgo <- Atgo[Atgo$Ontology=="BP",]
  #convert to list
  Atgo.list <- tapply(Atgo$go_id,Atgo$gene_id,c)
  #collapse list
  Atgo.terms <- data.frame(AtGO=unlist(lapply(Atgo.list,paste,collapse=";")),AGI=names(Atgo.list))
  results.annotated <- merge(results.annotated,Atgo.terms,by="AGI",all.x=T,sort=F)
  
  
  #add At gene symbol
  AtSymbol <- toTable(org.At.tairSYMBOL)
  AtSymbol <- AtSymbol[!duplicated(AtSymbol$gene_id),]
  results.annotated <- merge(results.annotated,AtSymbol,by.x="AGI",by.y="gene_id",all.x=T,sort=F)
  
  #Add At description
  AtNames <- toTable(org.At.tairGENENAME)
  AtNames <- AtNames[!duplicated(AtNames$gene_id),]
  anyDuplicated(AtNames$gene_id)
  results.annotated <- merge(results.annotated,AtNames,by.x="AGI",by.y="gene_id",all.x=T,sort=F)
  head(results.annotated)
  write.csv(results.annotated,sample)
}
############




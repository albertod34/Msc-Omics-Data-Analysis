
#loading data

library(SummarizedExperiment)
load("gbm.Rdata")

pheno <- colData(rse)
names(pheno)
summary(pheno$Cluster)



# let's select only the 4 subtypes  ...

mask <- !pheno$Cluster%in%c("G-CIMP", "NA") & !is.na(pheno$Cluster)

#checking that the mask works

pheno.chnged <- pheno[mask,]
names(pheno.chnged)
summary(pheno.chnged$Cluster)

# final object (subset using the mask) with samples of interest
rse.s <- rse[,mask]
rse.s

# get counts
names(assays(rse.s))
cc <- assays(rse.s)$"raw_counts"
dim(cc)

# get group
group <- colData(rse.s)$Cluster
table(group)


#Normalization with TMM

library(tweeDEseq)
help("filterCounts")
counts.f <- filterCounts(cc)
dim(cc)
dim(counts.f)
counts.tmm <- normalizeCounts(counts.f, method="TMM")
dim(counts.tmm)

#MA plot to check normalization

#before normalization
library(edgeR)
maPlot(counts.f[,3], counts.f[,4], pch=19, cex=.5, ylim=c(-8,8),
       allCol="darkgray", lowess=TRUE,  
       xlab=expression( A == log[2] (sqrt(Sample1 %.% Sample2))  ),
       ylab=expression(M == log[2](Sample1)-log[2](Sample2)))
grid(col="black")

#after normalization
maPlot(counts.tmm[,3], counts.tmm[,4], pch=19, cex=.5, ylim=c(-8,8),
       allCol="darkgray", lowess=TRUE,  
       xlab=expression( A == log[2] (sqrt(Sample1 %.% Sample2))  ),
       ylab=expression(M == log[2](Sample1)-log[2](Sample2)))
grid(col="black")

#DEGs with edgeR

names(pheno.chnged)
r<- pheno.chnged$Cluster
d<-DGEList(counts=counts.tmm, group= as.factor(r))
names(d)
head(pheno.chnged$Cluster)
d
d<-estimateCommonDisp(d)
names(d)
d<-estimateTagwiseDisp(d)   

  ##Classical vs Neural

Class_Neural_EdgeR.tagwise <- exactTest(d, pair=c("Classical", "Neural"), dispersion="tagwise")

topTags(Class_Neural_EdgeR.tagwise)

nrow(topTags(Class_Neural_EdgeR.tagwise , n=Inf, p=0.05))

  ##Classical vs Mesenchymal

Class_Mesenchymal_EdgeR.tagwise <- exactTest(d, pair=c("Classical", "Mesenchymal"), dispersion="tagwise")

topTags(Class_Mesenchymal_EdgeR.tagwise)

nrow(topTags(Class_Mesenchymal_EdgeR.tagwise , n=Inf, p=0.05))

  ##Classical vs Proneural

Class_Proneural_EdgeR.tagwise <- exactTest(d, pair=c("Classical", "Proneural"), dispersion="tagwise")

topTags(Class_Proneural_EdgeR.tagwise)

nrow(topTags(Class_Proneural_EdgeR.tagwise , n=Inf, p=0.05))


##Enrichment Analysis

library(GOstats)
library(org.Hs.eg.db)


 #Neural

tt.cn<-topTags(Class_Neural_EdgeR.tagwise, n=Inf)

mask.n <- tt.cn$table$FDR < 0.01 & abs(tt.cn$table$logFC) > log2(2)
deGenes.n <- rownames(tt.cn$table[mask.n, ])
head(deGenes.n)
length(deGenes.n)

  #get the gene IDs
deGenes.n.s<-unlist(strsplit(deGenes.n,"|",fixed = TRUE))
head(deGenes.n.s)
deGenes.n.s<-as.integer(deGenes.n.s)
head(deGenes.n.s)
mask.ns<-!is.na(deGenes.n.s)
deGenes.n.s<- deGenes.n.s[mask.ns]
head(deGenes.n.s)

  ##get the IDs for gene Universe
geneUniverse.n <- rownames(tt.cn$table)
length(geneUniverse.n)
geneUniverse.n.s<-unlist(strsplit(geneUniverse.n,"|",fixed = TRUE))
head(geneUniverse.n.s)
geneUniverse.n.s<-as.integer(geneUniverse.n.s)
head(geneUniverse.n.s)
mask.guns<-!is.na(geneUniverse.n.s)
geneUniverse.n.s<- geneUniverse.n.s[mask.guns]
head(geneUniverse.n.s)
length(geneUniverse.n.s)

params.n <- new("GOHyperGParams", geneIds=deGenes.n.s, universeGeneIds=geneUniverse.n.s,annotation="org.Hs.eg.db", ontology="BP",
          pvalueCutoff=0.05, conditional=TRUE, testDirection="over")

hgOver.n <- hyperGTest(params.n)
hgOver.n

htmlReport(hgOver.n, file="gbm_neural.html")


  #Mesenchymal


tt.cm<-topTags(Class_Mesenchymal_EdgeR.tagwise, n=Inf)

mask.m <- tt.cm$table$FDR < 0.01 & abs(tt.cm$table$logFC) > log2(2)
deGenes.m <- rownames(tt.cm$table[mask.m, ])
head(deGenes.m)
length(deGenes.m)

#get the gene IDs
deGenes.m.s<-unlist(strsplit(deGenes.m,"|",fixed = TRUE))
head(deGenes.m.s)
deGenes.m.s<-as.integer(deGenes.m.s)
head(deGenes.m.s)
mask.ms<-!is.na(deGenes.m.s)
deGenes.m.s<- deGenes.m.s[mask.ms]
head(deGenes.m.s)


##get the IDs for gene Universe
geneUniverse.m <- rownames(tt.cm$table)
length(geneUniverse.m)
geneUniverse.m.s<-unlist(strsplit(geneUniverse.m,"|",fixed = TRUE))
head(geneUniverse.m.s)
geneUniverse.m.s<-as.integer(geneUniverse.m.s)
head(geneUniverse.m.s)
mask.gums<-!is.na(geneUniverse.m.s)
geneUniverse.m.s<- geneUniverse.m.s[mask.gums]
head(geneUniverse.m.s)
length(geneUniverse.m.s)

params.m <- new("GOHyperGParams", geneIds=deGenes.m.s, universeGeneIds=geneUniverse.m.s,annotation="org.Hs.eg.db", ontology="BP",
                pvalueCutoff=0.05, conditional=TRUE, testDirection="over")

hgOver.m <- hyperGTest(params.m)
hgOver.m

htmlReport(hgOver.m, file="gbm_mesenchymal.html")

  #Proneural

tt.cp<-topTags(Class_Proneural_EdgeR.tagwise, n= Inf)


mask.p <- tt.cp$table$FDR < 0.01 & abs(tt.cp$table$logFC) > log2(2)
deGenes.p <- rownames(tt.cp$table[mask.p, ])
head(deGenes.p)
length(deGenes.p)

#get the gene IDs
deGenes.p.s<-unlist(strsplit(deGenes.p,"|",fixed = TRUE))
head(deGenes.p.s)
deGenes.p.s<-as.integer(deGenes.p.s)
head(deGenes.p.s)
mask.ps<-!is.na(deGenes.p.s)
deGenes.p.s<- deGenes.p.s[mask.ps]
head(deGenes.p.s)


##get the IDs for gene Universe
geneUniverse.p <- rownames(tt.cp$table)
length(geneUniverse.p)
geneUniverse.p.s<-unlist(strsplit(geneUniverse.p,"|",fixed = TRUE))
head(geneUniverse.p.s)
geneUniverse.p.s<-as.integer(geneUniverse.p.s)
head(geneUniverse.p.s)
mask.gups<-!is.na(geneUniverse.p.s)
geneUniverse.p.s<- geneUniverse.p.s[mask.gups]
head(geneUniverse.p.s)
length(geneUniverse.p.s)
params.p <- new("GOHyperGParams", geneIds=deGenes.p.s, universeGeneIds=geneUniverse.p.s,annotation="org.Hs.eg.db", ontology="BP",
                pvalueCutoff=0.05, conditional=TRUE, testDirection="over")

hgOver.p <- hyperGTest(params.p)
hgOver.p

htmlReport(hgOver.p, file="gbm_proneural.html")

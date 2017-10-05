
##1## Quality assessment for Rat data

#loading the required libraries

library(GEOquery) 
library(oligo)

#obtaining data from GEO and loading it to R session

gseid_Rat<-"GSE71083"
gse <- getGEO(gseid_Rat,GSEMatrix=TRUE) 

GEO<-getGEOSuppFiles(gseid_Rat,makeDirectory=TRUE)

untar(file.path(getwd(),gseid_Rat,"GSE71083_RAW.tar"), exdir = file.path(getwd(),gseid_Rat))

setwd(file.path(getwd(),gseid_Rat))
celFiles = list.files( pattern = "CEL.gz")

GEOFS <- read.celfiles(celFiles)

# exploring data

GEOFS 
class(GEOFS)
slotNames(GEOFS)

sampleNames(GEOFS) 
sampleNames(GEOFS)<-gsub(".CEL.gz","",sampleNames(GEOFS))

# image analysis

image(GEOFS)

#we don't see any problem in the images so we continue with the analysis

# quality assessment 
colos<- c("blue","red","green","yellow")
hist(GEOFS,col=colos)
legend("topright",sampleNames(GEOFS),fill=colos)
boxplot(GEOFS,col=colos)

#MA plots. One vs others and one by one:
MAplot(GEOFS)
MAplot(GEOFS,pairs=TRUE)
#loess of fit is what we expected

GEOFS.fit<-fitProbeLevelModel(GEOFS)
RLE(GEOFS.fit,col=colos)
# Medians are around 0 so it's correct.
NUSE(GEOFS.fit,col=colos)
# Medians are around 1 so it's correct.

##1## Quality assessment for Mouse data

#obtaining data from GEO and loading it to R session

gseid_Mouse<-"GSE71082"
gse_m <- getGEO(gseid_Mouse,GSEMatrix=TRUE) 

GEO_m<-getGEOSuppFiles(gseid_Mouse,makeDirectory=TRUE)

untar(file.path(getwd(),gseid_Mouse,"GSE71082_RAW.tar"), exdir = file.path(getwd(),gseid_Mouse))

setwd(file.path(getwd(),gseid_Mouse))
celFiles2 = list.files( pattern = "CEL.gz")

GEOFS_M <- read.celfiles(celFiles2)

# exploring data

GEOFS_M 
class(GEOFS_M)
slotNames(GEOFS_M)

sampleNames(GEOFS_M) 
sampleNames(GEOFS_M)<-gsub(".CEL.gz","",sampleNames(GEOFS_M))

# image analysis 

image(GEOFS_M)

#we don't see any problem in the images so we continue with the analysis

# quality assessment 

hist(GEOFS_M,col=colos)
legend("topright",sampleNames(GEOFS_M),fill=colos)
boxplot(GEOFS_M,col=colos)

#MA plots. One vs others and one by one:
MAplot(GEOFS_M)
MAplot(GEOFS_M,pairs=TRUE)
#loess of fit is what we expected

# RLE and NUSE plots:

GEOFS.fit_M<-fitProbeLevelModel(GEOFS_M)
RLE(GEOFS.fit_M,col=colos)
# Medians are around 0 so it's correct.

NUSE(GEOFS.fit_M,col=colos)
# Medians are around 1 so it's correct.

#General observation of step 1: In the GEO Series GSE71084 we found data splitted in 2 different GEO Series: 
#GSE71082 for mouse and GSE71083 for rat. 

##2## Normalization for Rat data

#normalization using rma method

GEOFS.rma<-rma(GEOFS,target= "core")
boxplot(GEOFS.rma,col=colos)


##2## Normalization for Mouse data

GEOFS_M.rma<-rma(GEOFS_M,target= "core")
boxplot(GEOFS_M.rma,col=colos)

##data aggregation for Rat

use.cor="pairwise.complete.obs"
x<-exprs(GEOFS.rma)
clust.cor.ward<- hclust(as.dist(1-cor(x,use=use.cor)),method="ward.D2")
plot(clust.cor.ward, main="hierarchical clustering", hang=-1,cex=0.6)

##data aggregation for Mouse


y<-exprs(GEOFS_M.rma)
clust.cor.ward_M<- hclust(as.dist(1-cor(y,use=use.cor)),method="ward.D2")
plot(clust.cor.ward_M, main="hierarchical clustering", hang=-1,cex=0.6)

## in both cases no need to correct for batch effects. The clustering looks quite good.

##3## DEG for rat data

#library of linear models for microarray analysis (limma) has to be loaded in order to perform the DEG detection
library(limma)

#We generate our three needed objects, ExpressionSet, the design matrix based on conditions and the contrasts:
cond1<-pData(GEOFS)$index[1:2]<-c(1,1)
cond2<-pData(GEOFS)$index[3:4]<-c(2,2)
cond<- as.factor(c(cond1,cond2))
design<-model.matrix(~0+cond)
rownames(design)<-sampleNames(GEOFS)

fit<-lmFit(GEOFS.rma,design) 
contrast.matrix<-makeContrasts(cond2-cond1,levels=design)
fit2<-contrasts.fit(fit,contrast.matrix)
fite<-eBayes(fit2)

top.table<-topTable(fite,coef=1,number=Inf,adjust="BH") # adjusted for multiple comparisons with Benjamini&Hochberg approach
results<-decideTests(fite)
table(results)
table(decideTests(fite,adjust.method="none"))

#once we have obtained our results we adjust the p-value < 0.05

results.adj.p0.05<-top.table[top.table$adj.P.Val<0.05,]
dim(results.adj.p0.05)

# We obtain 2156 DEGs.

range(results.adj.p0.05$logFC)

##3## DEG for mouse data

cond1m<-pData(GEOFS_M)$index[1:2]<-c(1,1)
cond2m<-pData(GEOFS_M)$index[3:4]<-c(2,2)
condm<- as.factor(c(cond1m,cond2m))
designm<-model.matrix(~0+condm)
rownames(designm)<-sampleNames(GEOFS_M)

fitm<-lmFit(GEOFS_M.rma,design)
contrast.matrixm<-makeContrasts(cond2-cond1,levels=design)
#contrasts matrix produce a strange result. we do the same as above but not the same is obtained
#> contrast.matrixm
#Contrasts
#Levels  cond2m - cond1m
#cond1               1
#cond2               1
#we change the conditions and put cond2 and cond1 in order to get good results.
fit2m<-contrasts.fit(fitm,contrast.matrixm)
fitem<-eBayes(fit2m)
top.tablem<-topTable(fitem,coef=1,number=Inf,adjust="BH")

resultsm<-decideTests(fitem)
table(resultsm)
table(decideTests(fitem,adjust.method="none"))

results.adj.p0.05m<-top.tablem[top.tablem$adj.P.Val<0.05,]
dim(results.adj.p0.05m)

#we don't find any results with the adjusted p-value
#So we can try to find some DEGs using just the p-value with no correction made.

results.p0.05m<-top.tablem[top.tablem$P.Val<0.05,]
dim(results.p0.05m)

#we now obtain 2283 DEGs.

## 4. Generate a Volcano plot and a heat map 
#volcano plot Rat
volcanoplot(fite,coef=1,highlight=10,names=fite$genes$NAME,main="Volcano plot")

#heat map rat
library(gplots)
data.clus<-exprs(GEOFS.rma[rownames(results.p0.05),])
clust.rows <- hclust(as.dist(1-cor(t(data.clus))),method="ward.D2")
clust.cols <- hclust(as.dist(1-cor(data.clus)),method="ward.D2")  
heatcol<-colorRampPalette(c("green", "Black","red"), space = "rgb")

heatm<-heatmap.2(as.matrix(data.clus), col = heatcol(256),
                 dendrogram="column", Colv=as.dendrogram(clust.cols),
                 Rowv=as.dendrogram(clust.rows),
                 scale="row",cexRow=0.1, cexCol=0.5, 
                 main="",key=TRUE,keysize=1,density.info="none",trace="none")



#volcano plot for mouse
volcanoplot(fitem,coef=1,highlight=10,names=fitem$genes$NAME,main="Volcano plot")

#heat map for mouse
library(gplots)
data.clus.m<-exprs(GEOFS_M.rma[rownames(results.p0.05m),])
clust.rows.m <- hclust(as.dist(1-cor(t(data.clus.m))),method="ward.D2")
clust.cols.m <- hclust(as.dist(1-cor(data.clus.m)),method="ward.D2")  
heatcol.m<-colorRampPalette(c("green", "Black","red"), space = "rgb")

heatm.m<-heatmap.2(as.matrix(data.clus.m), col = heatcol.m(256),
                 dendrogram="column", Colv=as.dendrogram(clust.cols.m),
                 Rowv=as.dendrogram(clust.rows.m),
                 scale="row",cexRow=0.1, cexCol=0.5, 
                 main="",key=TRUE,keysize=1,density.info="none",trace="none")

## 5. Annotation of results

library(annotate)

#specific package for mouse annotation

source("https://bioconductor.org/biocLite.R")
biocLite("mogene10sttranscriptcluster.db")
biocLite("ragene10sttranscriptcluster.db")
library(mogene10sttranscriptcluster.db)
library(ragene10sttranscriptcluster.db)

#Getting all the values of mouse data with wich we want to build the annotation:

datm <- exprs(GEOFS_M.rma)[rownames(results.p0.05m),]
logFCm <- results.p0.05m$logFC
pvalm <- results.p0.05m$P.Value
adj.pvalm<-results.p0.05m$adj.P.Val

sym.m<-unlist(mget(rownames(results.p0.05m), env=mogene10sttranscriptclusterSYMBOL))
name.m<-unlist(mget(rownames(results.p0.05m), env=mogene10sttranscriptclusterGENENAME))
chr.m<-unlist(mget(rownames(results.p0.05m), env=mogene10sttranscriptclusterCHR))

chr.mouse.1 <- match (rownames (results.p0.05m), names(chr.m))
chr.mouse.2 <- chr.m[chr.mouse.1]

affyids.m<-rownames(results.p0.05m)
genelist.m <- list(affyids.m)
filename.m <- "Results Mouse.html"
title.m <- "Differentially expressed fibrinogen vs control"
othernames.m <- list(sym.m,name.m,chr.mouse.2,round(logFCm, 1), round(pvalm, 4), round(adj.pvalm, 4), round(datm, 2)) 
head.m <- c("Probe ID", "Symbol", "Gene Name", "Chr","logFC", "p-value","adj.p-value",sampleNames(GEOFS_M.rma))
repository.m <- list("affy")

# we generate then the html page in order to have results well presented

htmlpage(genelist.m, filename.m, title.m, othernames.m, head.m, repository = repository.m)

#Getting all the values of rat data with wich we want to build the annotation:

dat <- exprs(GEOFS.rma)[rownames(results.adj.p0.05),]
logFC <- results.adj.p0.05$logFC
pval <- results.adj.p0.05$P.Value
adj.pval<-results.adj.p0.05$adj.P.Val

sym<-unlist(mget(rownames(results.adj.p0.05), env=ragene10sttranscriptclusterSYMBOL))
name<-unlist(mget(rownames(results.adj.p0.05), env=ragene10sttranscriptclusterGENENAME))
chr<-unlist(mget(rownames(results.adj.p0.05), env=ragene10sttranscriptclusterCHR))


affyids<-rownames(results.adj.p0.05)
genelist <- list(affyids)
filename <- "Results Rat.html"
title <- "Differentially expressed fibrin vs control"
othernames <- list(sym,name,chr,round(logFC, 1), round(pval, 4), round(adj.pval, 4), round(dat, 2)) 
head <- c("Probe ID", "Symbol", "Gene Name", "Chr","logFC", "p-value","adj.p-value",sampleNames(GEOFS.rma))
repository<- list("affy")

# we generate then the html page in order to have results well presented

htmlpage(genelist, filename, title, othernames, head, repository = repository)

## 6. Functional analysis using GOstats

library(GOstats)
library(GSEABase)

#mouse data:

library(org.Mm.eg.db)

frame = toTable(org.Mm.egGO)
goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
head(goframeData)

goFrame=GOFrame(goframeData,organism="Mus musculus")
goAllFrame=GOAllFrame(goFrame)

gsc.GO <- GeneSetCollection(goAllFrame, setType = GOCollection())

#we will use entrez ids
entrez<-unlist(mget(featureNames(GEOFS_M.rma), env=mogene10sttranscriptclusterSYMBOL))
length(entrez)

#universe
length(unique(entrez))
entrezuniverse<-unique(entrez)

#results 
results.entrez<-unique(unlist(mget(rownames(results.p0.05m), env=mogene10sttranscriptclusterSYMBOL)))

cutoff=0.001
GO.params.bp <- GSEAGOHyperGParams(name="GOstats",  geneSetCollection=gsc.GO, geneIds = results.entrez, universeGeneIds = entrezuniverse,  ontology = "BP", pvalueCutoff = cutoff,  conditional = FALSE, testDirection = "over")
GO.results.bp<-hyperGTest(GO.params.bp)
head(summary(GO.results.bp))
GO.results.bp
htmlReport(GO.results.bp, "GO.results.BiologicalProcess.html")










library(minfi)
setwd("~/AGomez_practica")


idat.folder <-getwd()

targets <- read.450k.sheet(base=idat.folder)



#Loading data

rgset <- read.450k.exp(targets = targets)

#Quality control 

MSet <- preprocessRaw(rgset) 

RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)

gc()
GRset <- mapToGenome(RSet)
gc()

qc <- getQC(MSet)

plotQC(qc)

phenoData <- pData(rgset)

densityPlot(MSet, sampGroups = phenoData$Status)

#Sex prediction

predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
head(predictedSex)

pdf("sex.pdf")
plotSex(getSex(GRset, cutoff = -2))
dev.off()

#Detection of p-values

detP <- detectionP(rgset)

pdf("pvals.detect.pdf")
barplot(colMeans(detP),col=factor(targets$Status),las=2,cex.names=0.8,
        main="Mean detection p-values")
dev.off()

  #cutoff 0.05

nkeep <- colMeans(detP) > 0.05
nkeep
table(nkeep)

#Normalization with SQN

gc()    
gRatioSet.quantile <- preprocessQuantile(rgset)
gc()    


#remove probes in chr X and Y
ann450k<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
keep <- !(featureNames(gRatioSet.quantile) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)

gRatioSet.quantile <- gRatioSet.quantile[keep,]

# remove probes with SNPs at CpG 
gRatioSet.quantile<- dropLociWithSnps(gRatioSet.quantile)

betas<-getBeta(gRatioSet.quantile)

colnames(betas)<-targets$Status


ann450k<-as.data.frame(ann450k)
gc()
dades <- merge(betas,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
rownames(dades)<-dades$Row.names
dades<-dades[,-1]


##DMP finding


dmp <- dmpFinder(betas, pheno = targets$Status  , type = "categorical")
dmp<-subset(dmp,dmp$qval<.01)

gc()

Adult <- grep("^Ad",colnames(dades))
Fetal <- grep("^Fe",colnames(dades))
dades$Mean_Adult <- rowMeans(dades[,Adult],na.rm=T)
dades$Mean_Fetal <- rowMeans(dades[,Fetal],na.rm=T)

final.dades<- dades[row.names(dades)%in%row.names(dmp),]


res<-final.dades[abs(final.dades$Mean_Adult-final.dades$Mean_Fetal)>.2,]

dim(res)

#DM in promoters

proms<-res[grep("TSS1500|TSS200|5'UTR|1stExon",res$UCSC_RefGene_Group),]
dim(proms)

#DM in promoters and islands 

proms.island<-proms[grep("Island",proms$Relation_to_Island),]
dim(proms.island)

# CpGs hypermethylated and hypomethylated

hpr_fetal <- proms.island[proms.island$Mean_Fetal > 0.66,]
dim(hpr_fetal)
hyp_adult <- hpr_fetal[hpr_fetal$Mean_Adult < 0.33,]
dim(hyp_adult)


#Differentially expressed genes
setwd("~/")
d_expr <- read.csv("diff.expressed.AdultvsFetal.csv")
head(d_expr)
dim(d_expr)
d_expr <- d_expr[,-c(1)]
dim(d_expr)
d_expr_o2 <- d_expr[d_expr$logFC>2,]

d_expr_symb <- unique(unlist(strsplit(as.character(d_expr_o2$Gene.symbol[d_expr_o2$Gene.symbol != ""]),'///')))
pi_symb <- unique(unlist(strsplit(proms.island$UCSC_RefGene_Name,';')))


d_expr_final <- match(d_expr_symb,pi_symb)

d_expr_final <- d_expr_final[!is.na(d_expr_final)]

DEGS <- proms.island[d_expr_final,]
dim(DEGS)

##GO enrichment
#####GO enrichment analysis
gc()

library(org.Hs.eg.db)
library(GOstats)
library(GO.db)
library(annotate)


ids<- unlist(mget(as.character(pi_symb),ifnotfound=NA, revmap(org.Hs.egSYMBOL)))

univ <- Lkeys(org.Hs.egGO)

paramBP <- new("GOHyperGParams", geneIds=ids, universeGeneIds=univ, annotation="org.Hs.eg.db", 
               ontology="BP",pvalueCutoff= 0.01, conditional=FALSE,testDirection="over")
gc()
hypBP <- hyperGTest(paramBP)

# Get the p-values of the test (BP)
gGhyp.pv <- pvalues(hypBP)
gGhyp.odds<-oddsRatios(hypBP)
gGhyp.counts<-geneCounts(hypBP)

sigGO.ID <- names(gGhyp.pv[gGhyp.pv < 0.01])


# Test the number of counts

gGhyp.counts<-as.data.frame(gGhyp.counts)
gGhyp.counts$GOterms<-rownames(gGhyp.counts)
gGhyp.counts<-gGhyp.counts[rownames(gGhyp.counts) %in% sigGO.ID,]

#S ignificant GO terms of BP (Molecular Function)

sigGO.Term <- getGOTerm(sigGO.ID)[["BP"]]

results_GO<-cbind(as.data.frame(gGhyp.pv[gGhyp.pv < 0.01]), as.data.frame(sigGO.Term), gGhyp.counts)

write.csv(results_GO, "GO_enrichment_BP_Adult_vs_Fetal.csv")

# HeatMap

library(gplots)
l <- as.matrix(proms.island[,c(Adult,Fetal)])
spcol <- c(rep("blue1",length(Adult)),rep("red",length(Fetal)))
heatmap.2(l,
          main="Diffmeth CpG's ranked",
          labRow=NA,
          trace="none",
          na.rm=T,
          col=greenred,
          ColSideColors=spcol,
          distfun=function(x) dist(x,method="manhattan"),
          dendrogram = "column")


####################################
###Integrative genomics practical###
####################################

setwd("~/")

### 1. From the TCGA project choose a specific cancer type and report how many cases for each technology were processed.


# The selected cancer type is: Esophageal carcinoma (ESCA)

# Reported number of cases by technology:

# CopyNumber Analysis : 184.
# Clinical Analysis: 185.
# Methylation : 185.
# miRSeq: 184.
# mRNASeq: 184.
# LowPass DNASeq CopyNum: 51.
# Reverse Phase Protein Array: 126.
# Mutation Annotation File: 185. 


### 2. Using package RTCGAToolbox download from GDAC Firehose three different sorts of data
### (for instance, RNASeq, CNV and miRNA). (If RTCGAToolbox would not work try TCGAbiolinks).



library(RTCGAToolbox)

getFirehoseDatasets() # Looking for Esophageal carcinoma (ESCA)

getFirehoseRunningDates() # Deciding the date

?getFirehoseData


rdata<- getFirehoseData(dataset = 'ESCA', runDate = "20150821", Clinic = TRUE, RNAseq_Gene = TRUE,
                        miRNASeq_Gene = TRUE, CNV_SNP = TRUE )


#save(rdata, file = "ESCATCGA.rdata") 
load('ESCATCGA.rdata')

#Retrieving clinical data
rclinic<-getData(rdata,"Clinical")
dim(rclinic)
#185 20
rclinic[1:3,]
names(rclinic)


#Retrieving RNAseq data
rRNAseq<-getData(rdata,"RNASeqGene")
dim(rRNAseq)
#26120   198
rRNAseq[1:3,1:3]


#Retrieving miRNAseq data
rmiRNAseq<-getData(rdata,"miRNASeqGene")
dim(rmiRNAseq)
#1046  198
rmiRNAseq[1:3,1:3]

#Retrieving CNVSNP data
rCNV<-getData(rdata,"CNVSNP")
dim(rCNV)
#60803     6
rCNV[1:4,]


### 3. Apply MFA to the downloaded data sets using all cases or just a subset of them in case there are too many.
### Selection of cases, if applied, must be based on specific characteristics related to clinical variables (for
### instance pathology state, gender, etc.). Keep in mind that you will require several steps for preparing data by
### matching selected cases in the databases, selecting necessary variables, transposing, etc. All these
### steps need to be justified in the script and/or final report.

## Exploring clinical information for possible subsetting
dim(rclinic)
class(rclinic)
names(rclinic)
head(rownames(rclinic))
length(unique(rownames(rclinic)))

table(rclinic$gender) # selection criteria by gender
# female   male 
# 27    158 

# we keep female gender

clinic_fem<- subset(rclinic,rclinic$gender== 'female')

# obtaining specific id by participant for each tcga barcode

namessamples<-sapply(rownames(clinic_fem),function(x){unlist(strsplit(x,split = '.',fixed = T))[3]}) #we keep the third element of the barcode cause it's the participant ID
names(namessamples)<-NULL
rownames(clinic_fem)<- namessamples

## Managing RNAseq data

dim(rRNAseq)
class(rRNAseq)
prepRNAdat<- as.data.frame(t(rRNAseq)) #transpose and convert to a data.frame in order to have samples on rownames
dim(prepRNAdat)

sampsRNAdat<-sapply(rownames(prepRNAdat),function(x){tolower(unlist(strsplit(x,split = '-',fixed = T)))[3]}) #we keep the third element of the barcode cause it's the participant ID


whichdup<-duplicated(sampsRNAdat) #find duplicated samples

#check lenghts
length(whichdup)
#198

length(sampsRNAdat)
#198

length(rownames(prepRNAdat))
#198

RNAnodups<- prepRNAdat[!whichdup,] #keep only one type of sample per patient ID
dim(RNAnodups)
# 184 26120
rownames(RNAnodups)<-unique(sampsRNAdat) # change rownames to be able to do subsetting by gender

RNAseqfem<-RNAnodups[rownames(clinic_fem),] # select females RNAseq data only
dim(RNAseqfem) #check there are 27 samples
#27 26120

#anotation to reduce variables
library(biomaRt)

listMarts()

mart<-useMart(biomart = 'ENSEMBL_MART_ENSEMBL', host = 'www.ensembl.org')

listDatasets(mart = mart)

hsmart<-useMart(biomart = 'ENSEMBL_MART_ENSEMBL', host = 'www.ensembl.org',dataset = 'hsapiens_gene_ensembl')

#looking for genes in a exact chromosome
head(listAttributes(hsmart))
grep("chromosome",listAttributes(hsmart)[,1], value = T)
grep("gene_name",listAttributes(hsmart)[,1], value = T)
grep("hgnc",listFilters(hsmart)[,1], value = T)
RNASeqanot <- getBM(c("external_gene_name", "chromosome_name"), values = colnames(RNAseqfem),
                    filters = "hgnc_symbol",mart = hsmart)

class(RNASeqanot)
# "data.frame"
dim(RNASeqanot)
# 22241     2
head(RNASeqanot)

# Selection of chromosome 13, due to literature: http://www.ncbi.nlm.nih.gov/pubmed/24510239
RNA13<- subset(RNASeqanot, RNASeqanot$chromosome_name=='13')
dim(RNA13)
# 349   2
head(RNA13)

rnasf<- RNAseqfem[, RNA13$external_gene_name] #final object with RNAseq data for females in 13th chromosome
dim(rnasf)
# 27 349

###
### Managing miRNAseq data

dim(rmiRNAseq)
class(rmiRNAseq)
prepmiRNAdat<- as.data.frame(t(rmiRNAseq)) #transpose and convert to a data.frame in order to have samples on rownames
dim(prepmiRNAdat)
# 198 1046
sampsmiRNAdat<-sapply(rownames(prepmiRNAdat),function(x){tolower(unlist(strsplit(x,split = '-',fixed = T)))[3]}) #we keep the third element of the barcode cause it's the participant ID


whichdupmi<-duplicated(sampsmiRNAdat) #find duplicated samples

#check lenghts
length(whichdupmi)
#198

length(sampsmiRNAdat)
#198

length(rownames(prepmiRNAdat))
#198

miRNAnodups<- prepmiRNAdat[!whichdupmi,] #keep only one type of sample per patient ID
dim(miRNAnodups)
# 184 1046
rownames(miRNAnodups)<-unique(sampsmiRNAdat) # change rownames to be able to do subsetting by gender

miRNAseqfem<-miRNAnodups[rownames(clinic_fem),] # select females miRNAseq data only 
dim(miRNAseqfem)#check there are 27 samples
# 27 1046

# anotation
grep("mir",listFilters(hsmart)[,1], value = T)


miRNASeqanot <- getBM(c("mirbase_id", "chromosome_name"), values = colnames(miRNAseqfem),
                      filters = "mirbase_id",mart = hsmart)

class(miRNASeqanot)

dim(miRNASeqanot)
#833   2
head(miRNASeqanot)

miRNA13<- subset(miRNASeqanot, miRNASeqanot$chromosome_name=="13")
dim(miRNA13)
# 20   2
head(miRNA13)
mirnasf<- miRNAseqfem[, miRNA13$mirbase_id]#final object with miRNAseq data for females in 13th chromosome
dim(mirnasf)
# 27 20


###
### Managing CNVSNP data

dim(rCNV)
#60803 6
class(rCNV)
head(rCNV)

# Selection of chromosome 13, due to literature: http://www.ncbi.nlm.nih.gov/pubmed/24510239
CNVchrom13<- subset(rCNV, rCNV$Chromosome == 13)
dim(CNVchrom13)
#1918    6
head(CNVchrom13)

which.max(table(CNVchrom13$Start))# finding most repeated start
which.max(table(CNVchrom13$End)) #finding most repeated end

CNVchrom13start<- subset(CNVchrom13, CNVchrom13$Start ==19450806) #subsetting by position
CNVchrom13send<- subset(CNVchrom13start, CNVchrom13start$End ==114987458) # subsetting by position
rownames(CNVchrom13send)<- CNVchrom13send$Sample
nwrownames<-sapply(rownames(CNVchrom13send),function(x){tolower(unlist(strsplit(x,split = '-',fixed = T)))[3]}) #getting the IDs
whichdupCNV<-duplicated(nwrownames) #finding duplicates

length(whichdupCNV)
#138
length(nwrownames)
#138
length(rownames(CNVchrom13send))
#138

CNVnodups<- CNVchrom13send[!whichdupCNV,] #keep only one type of sample per patient ID
dim(CNVnodups)
# 112 6
rownames(CNVnodups)<-unique(nwrownames) # change rownames to be able to do subsetting by gender

CNVfem<-CNVnodups[rownames(clinic_fem),] # select females CNV data only 
dim(CNVfem)#check there are 27 samples
# 27  6

#apply MFA


#DATA preparation
#select samples

CNV4MFA<-CNVfem[!is.na(CNVfem[,1]),c(5,6)]
#which(apply(RNAs4MFA,2,mean) == 0) # know wich variables are 0 in all samples
RNAs4MFA<- rnasf[rownames(CNV4MFA),-c(121,211,214,223,346)]
#which(apply(miRNA4MFA,2,mean) == 0) # know which variables are 0 in all samples
miRNA4MFA<- mirnasf[rownames(CNV4MFA),c(2,3,4,5,6,7,10,11,12,15,16,17,20)]

identical(rownames(CNV4MFA),rownames(miRNA4MFA)) # check rownames

rna.l<-ncol(RNAs4MFA)
mirna.l<-ncol(miRNA4MFA)
cnv.l<-ncol(CNV4MFA)

data4Facto<-data.frame(RNAs4MFA,CNV4MFA,miRNA4MFA) # building the compact data.frame
dim(data4Facto)
#15 359

#Results
library(FactoMineR)
res <- MFA(data4Facto, group=c(rna.l,cnv.l,mirna.l), type=c(rep("c",3)), ncp=5, name.group=c("RNA","CNV","miRNA")) 


### 4. Comment on the MFA results.



# The graphics are provided in supplementary pdf, but comments displayed below using the order of appearance.

# In the last graphic, groups representation, our three kind of data don't seem to be really close. Looks like a 
# random distribution although it's true that are almost centered and not quite dispersed. An interesting
# thing is that the first dimension of the global PCA explains almost 30% of the initial inertia and the 
# second dimension explains the 20%. So in total we have 50% of inertia explained by the two first dimensions.
# In the first graphic, partial axes, we see five of the dimensions of the partial PCAs. In general we cannot
# say that any dimention fits exactly with dimension 1 or 2 except dimension 2 of CNV data that seem to inversely 
# correlate with dimension 2 of the global pca.
# The graphic of the correlation circle is completely illegible so it makes no sense to commentate due to presence
# of many variables and many dimensions.
# And finally, in the remaining two plots, representing individual factor map, we cannot denote tendencies since data are 
# scattered and ramdomly distributed.




### 5. Generate a CIRCOS plot showing the results of the analysis as well as the downloaded data, 
### justifying data filtering and selection criteria for each track of the plot.


library(OmicCircos)

# Managing downloaded rnaseq data

grep("chr",listAttributes(hsmart)[,1], value = T)[1]
grep("gene",listAttributes(hsmart)[,1], value = T)
grep("position",listAttributes(hsmart)[,1], value = T)
grep("hgnc",listFilters(hsmart)[,1], value = T)
rRNAseq[1:5,1:5]

RNAseqannotation<- getBM(c("external_gene_name", "chromosome_name", "start_position"),
                              filters = "hgnc_symbol",
                              values = rownames(rRNAseq), mart = hsmart)

head(RNAseqannotation) #look how it looks like
RNAseqannotation<- RNAseqannotation[RNAseqannotation$chromosome_name<3,]#getting only chromosomes of interest
RNAseqannotation<- RNAseqannotation[!duplicated(RNAseqannotation[,"external_gene_name"]),] #avoiding duplicates
rrna<- merge(RNAseqannotation,rRNAseq,by.x = 'external_gene_name',by.y='row.names') #merging datasets
rownames(rrna)<-rrna$external_gene_name

rrna[1:5,1:5]
class(rrna)


#Managing downloaded miRNAseq data

grep("mir",listFilters(hsmart)[,1], value = T)
rmiRNAseq[1:5,1:5]

mirnaannotation <- getBM(c("external_gene_name", "chromosome_name", "start_position","mirbase_id"),
                           filters = "mirbase_id",
                           values = rownames(rmiRNAseq), mart = hsmart)
head(mirnaannotation)#look how it looks like
mirnaannotation<- mirnaannotation[mirnaannotation$chromosome_name<3,]#getting only chromosomes of interest
mirnaannotation <- mirnaannotation[!duplicated(mirnaannotation$external_gene_name),] #avoiding duplicates
mirnaannotation <- mirnaannotation[!duplicated(mirnaannotation$mirbase_id),] #avoiding duplicates
mirna <- merge(mirnaannotation,rmiRNAseq,by.x = 'mirbase_id',by.y='row.names')#merging datasets
rownames(mirna)<-mirna$external_gene_name


mirna[1:5,1:5]
class(mirna)

#Managing downloaded CNVsnp data

grep("chrom",listFilters(hsmart)[,1], value = T)
head(rCNV)
cnvannot <- getBM(c("external_gene_name","chromosome_name", "start_position" ), filters = "chromosome_name",
                         values = rCNV$Chromosome, mart = hsmart)
head(cnvannot)#look how it looks like
dim(cnvannot)

cnvannot<- cnvannot[cnvannot$chromosome_name<3,]#getting only chromosomes of interest
cnvannot <- cnvannot[!duplicated(cnvannot$external_gene_name),]#avoiding duplicates
cnvsnp <- merge(cnvannot,rCNV,by.x = 'chromosome_name',by.y='row.names')#merging datasets
rownames(cnvsnp)<-cnvsnp$external_gene_name


cnvsnp[1:5,1:5]
class(cnvsnp)
dim(cnvsnp)



# Plot Circos:

#same name for columns
colnames(cnvsnp)[1]<-"chrom"
colnames(rrna)[2]<-"chrom"
colnames(mirna)[3]<-"chrom"
colnames(cnvsnp)[2]<-"name"
colnames(mirna)[2]<-"name"
colnames(rrna)[1]<- "name"
colnames(rrna)[3]<- "pos"
colnames(mirna)[4]<- "pos"
colnames(cnvsnp)[3]<- "pos"


# putting the same order for notation columns
mirna4circ<- mirna[,c(3,4,2,5:dim(mirna)[2])]
class(mirna4circ)
cnv4circ<- cnvsnp[,c(1,3,2,8)]
class(cnv4circ)
rrna4circ<-rrna[,c(2,3,1,4:dim(rrna)[2])]
class(rrna4circ)

#plot

pdf("CIRCOS Results.pdf", 8,8);
colors <- topo.colors(10, alpha=0.5);
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="Circos");
#Zoom to avoid sexual chromosomes
zoom <- c(1, 22, 1000000, 100000000, 0, 360);
# data

circos(R=400, type="chr", cir="hg18", col=colors, print.chr.lab=TRUE, W=4, scale=TRUE, zoom = zoom);

circos(R=325, cir="hg19", W=80, mapping=rrna4circ, col.v=4, type="s", B=F, lwd=1, col=colors[1],zoom = zoom);

circos(R=250, cir="hg19", W=80, mapping=mirna4circ,  col.v=4,  type="heatmap2",cluster=TRUE, col.bar=TRUE, lwd=1, col="blue",col.bar.po="topright",zoom=zoom);

circos (R=170, cir="hg19" , W=80, mapping=cnv4circ,col.v=4, type="ls" , B=F, col=colors[ 7 ] ,zoom = zoom, lwd=2, scale=TRUE);

dev.off()

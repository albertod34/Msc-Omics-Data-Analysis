


#EXERCISE: Researchers are interested in detecting SNPs associated with colon cancer (variable cascon of ’colon.txt’ ﬁle).
#To address this question, they performed a GWAS using DNA information about 2312 individuals.
#Genotype information (a random selection of 100,000 SNPs from the whole genome) is available in plink format (ﬁles ’colon.bed’, ’colon.bim’, ’colon.fam’ that are already available on Biocloud in the folder ’Data for exercises’).
#The phenotypic information can be found in the ﬁle ’colon.txt that includes these variables:

#id: identification number
#cascon: case-control status (0: control, 1:case)
#age: age in years
#smoke: smoking status
#bmi: body mass index
#ev3: 3rd principal component (population stratificaction)
#ev4: 4th principal component (population stratificaction)

setwd("~/Data_for_exercises")

library(snpStats)
library(SNPassoc)
library(GWASTools)
library(MASS)
library(PredictABEL)
library(biomaRt)
library(annotate)
#1. Perform a complete Genome-Wide Association Study (GWAS). The complete pipeline must
#include functions for:
 #1. Getting p-values of association between SNPs and colon cancer.

snps <- read.plink("colon")

names(snps)

geno<- snps$genotypes #genotype information

colon_ph <- read.delim("colon.txt") #colon individuals phenotypes
head(colon_ph)


identical(rownames(colon_ph),rownames(geno))#check if individuals are in the same order as in geno

rownames(colon_ph)<-colon_ph$id #putting the ids number in rownames

any(!rownames(colon_ph)%in%rownames(geno)) #check if the are the same individuals in both datasets 
#we obtain a false so there are same individuals but not at the same order. 
colon_ph <- colon_ph[rownames(geno),]# then we ordered them
identical(rownames(colon_ph), rownames(geno))
#we obtain a true so that means that all individuals are in the same exact order. 

info_snps<-col.summary(geno[colon_ph$cascon==0,])#we are interested in check the HWE in controls 
head(info_snps)

info2_snps<-row.summary(geno)
head(info2_snps)



hist(info_snps$MAF)
hist(info_snps$z.HWE) #chi-square values

plot(info2_snps) #Control Heterozygosity problems, we don't see it in the plot so we do not 
#need to filter it. We don't have negatives values as shown in the plot. 


as.test<- single.snp.tests(cascon, data=colon_ph, snp.data=geno) #association test

#quality control filtering
info_snps$pvalHWE <- 1-pnorm(info_snps$z.HWE) #calculate and keep the p-values
use <- info_snps$MAF > 0.01 & info_snps$pvalHWE > 0.001 #set the thresholds 
as.test <-as.test[use,]

#P-values corresponding to Cochran-Armitage trend test
pval <- p.value(as.test, df=1) 
plot(-log10(pval), col=ifelse(pval<0.0001, "red", "black"))


#Assessing whether population stratifcation is present or not. If so, perform a GWAS ad-
#justing for principal components (HINT: do not run a PCA, just use ev3 and ev4 variables if necessary).

chi2 <- chi.squared(as.test, df=1)
qq.chisq(chi2)#we can see that observed and expected values follow a lambda factor except a few. 
#So we can observe that there is not population stratification phenomenon 

#Creating a Manhattan-plot highlighting those SNP that are signifcant after addressing
#multiple comparisons.
annotation <- snps$map

head(annotation)


p.adj <- p.adjust(pval, method="fdr")
ans <- data.frame(SNP=names(as.test), pvalue=pval, fdr=p.adj)
ans.o <- ans[order(ans$pval),]
head(ans.o)




chromosome<- snps$map$chromosome
manhattanPlot(p.adj, chromosome[use], signif = 1e-7)

#Annotating the gene, chromosome, position and allele names for those significant SNPs that
#pass multiple comparisons.

snps.int<- ans.o[c(1,2),1]
snps.int
snps.int.c<-as.character(snps.int)

mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host="www.ensembl.org")


snpInfo <- getBM(c("refsnp_id", "chr_name", "chrom_start", "allele"),  
                 filters = c("snp_filter"), 
                 values = snps.int.c, mart = mart)
snpInfo

#Creating a Locus zoom plot for those SNPs that are significantly associated with colon
#cancer after addressing multiple comparisons problem.

candidate <- as.character(snps.int.c[1])
candidate

chr <- annotation[candidate, "chromosome"]
pos <- annotation[candidate, "position"] 
size <- 100000
mask <- annotation$chromosome == chr &
  annotation$position > pos - size &
  annotation$position < pos + size   
sum(mask)
snps.sel <- annotation[mask, "snp.name"]

head(ans.o)
info <- ans.o[ans.o$SNP%in%snps.sel, 1:2]
names(info) <- c("MarkerName", "P.value")
write.table(info, file="snpslz.txt", sep="\t",
            row.names=FALSE, quote=FALSE)

#Estimating the OR under dominant, recessive and additive model for those SNPs that pass
#multiple comparisons.


head(colon_ph)

a<-as(geno[,snps.int.c],"character")
head(a)

snpInfo


a[,1]<-gsub("A","C",a[,1])
a[,1]<-gsub("B","T",a[,1])
a[,2]<-gsub("A","A",a[,2])
a[,2]<-gsub("B","G",a[,2])
head(a)

phenosnps<-cbind(colon_ph,a)
head(phenosnps)
pheno.s<-setupSNP(phenosnps,8:ncol(phenosnps))
head(pheno.s)
ans <- WGassociation(cascon, pheno.s)
ORres<-WGstats(ans)
ORres$rs10112382
ORres$rs4733560

#Selecting the top-50 SNPs to create a genetic score and evaluating its predictive value.

head(ans.o)



geno.sel<-geno[,as.character(ans.o$SNP[1:50])]
head(geno.sel)
geno.sel

class(geno.sel[,1])



ww<-as.data.frame(geno.sel)
head(ww)
dim(ww)
ff<-function(x)
{
  xx<-as.numeric(x)-1
  xx[xx<0]<-NA
  return(xx)
}

ww.numeric<- data.frame(lapply(ww,ff))
head(ww.numeric)

attach(colon_ph)
head(colon_ph)

head(snps$fam)

dd1<-cbind(cascon,ww.numeric)
head(dd1)
dd.complete<-dd1[complete.cases(dd1),]
head(dd.complete)
ls(colon_ph)
mod<- stepAIC(glm(cascon~.,dd.complete,family = "gaussian"),method="forward")
summary(mod)

snps.score <- names(coef(mod))[-1] # Intercept is not going to be interpreted
pos <- which(names(dd.complete)%in%snps.score) # Get columns of the SNPs. 

score <- riskScore(mod, data=dd.complete, 
                   cGenPreds=pos, 
                   Type="unweighted")
table(score)

mod.lin <- glm(cascon~score, dd.complete, 
               family="gaussian")


summary(mod.lin)



coef(mod.lin)[2]



tval <- summary(mod.lin)[["coefficients"]]["score","t value"]

pval.score <- 1-pnorm(tval)

predrisk <- predRisk(mod.lin, dd.complete)

plotROC(data = dd.complete, cOutcome = 1, predrisk = predrisk)

#Performing pathway data analysis.

path.ann <- ans[,1:2]

write.table(path.ann, file="pvals.txt", sep="\t", 
            row.names=FALSE, quote=FALSE, 
            col.names=FALSE)


download.file("https://db.tt/LrZoog7j", 
              destfile="ICSNPathway_results.txt", mode="wb") 
file.show("ICSNPathway_results.txt")


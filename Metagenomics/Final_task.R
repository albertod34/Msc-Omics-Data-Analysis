setwd("~/")

#libraries

library(phyloseq)
packageVersion("phyloseq")
library(vegan)
packageVersion("vegan")
library(ggplot2)
packageVersion("ggplot2")
library(DESeq2)
packageVersion("DESeq2")
library(DESeq)

require(doBy)
packageVersion("doBy")
require(RColorBrewer)
packageVersion("RColorBrewer")
#######   gdata packages allows for microsofft office filetype read-in
library(gdata)


#### Since we are using 16s amplicon design it will be possible to 
### taxonomically  classify down to the genus level

AvailableRanks<-c("Kingdom","Phylum","Class","Order","Family","Genus")

### Define the minimum number of OTU counts for a sample to be further analyzed
minSampleCountsB<-1000
### Read the metadata
metadataB<-read.delim("Metadata_Evaluation.txt")


dim(metadataB)
metadataB[1:15,]
summary(metadataB)


#some changes to manage data easily 

head(rownames(metadataB))
rownames(metadataB)<-metadataB$SampleID
head(rownames(metadataB))

colnames(metadataB)

#Import OTU table

xB<-import_mothur(mothur_constaxonomy_file="Taxonomy_evaluation.txt",
                  mothur_shared_file="OTUTable_Evaluation.shared",
                  mothur_tree_file="phyloTree_Evaluation.tree",
                  cutoff=0.03)

xB #250 samples

# Attach the metadata
sample_data(xB)<-metadataB

# Define available ranks in taxonomy file
colnames(tax_table(xB))<-AvailableRanks

# Filter-off samples that are below 1000 counts (minSampleCountsB)
barplot(colSums(otu_table(xB)),las=2,cex.names=0.6,main="#Counts/Sample")
xB<-prune_samples(sample_sums(xB)>minSampleCountsB,xB) 
xB #242 samples

# Coverage barplots to pdf 
pdf("CoverageBarplots.pdf",paper="A4r")
barplot(colSums(otu_table(xB)),las=2,cex.names=0.6,main="#Counts/Sample")
dev.off()


#### 1 #### Which sampling site shows higher diversity? And richness?

x.1.0 <- xB

#Diversity plot
p<-plot_richness(x.1.0,"Site",measures=c("Shannon","Simpson","InvSimpson"))
p+geom_boxplot(aes(fill=Site))
ggsave("DiversityBySampling_Site.pdf")

#Richness plot
p<-plot_richness(x.1.0,"Site",measures=c("Observed","Chao1","ACE"))
p+geom_boxplot(aes(fill=Site))
ggsave("RichnessBySampling_Site.pdf")


my.sampleData<-data.frame(sample_data(x.1.0))
my.sampleData
er.x.1.0<-estimate_richness(x.1.0)
rownames(er.x.1.0)<-rownames(my.sampleData)

er.x.1.0<-merge(my.sampleData,er.x.1.0,by="row.names",all.x=T)
rownames(er.x.1.0)<-rownames(my.sampleData)
sample_data(x.1.0)<-er.x.1.0

# Generate Diversity Stats and capture them in a txt file

capture.output(file="DivAndRichnessStats.txt",paste())
for(covar in "Site") {
  for (measure in c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson")){
    capture.output(file="DivAndRichnessStats.txt",paste("############",measure,"by",covar),append=T)
    if(length(unique(er.x.1.0[,covar]))==2){
      myttest<-t.test(er.x.1.0[,measure]~er.x.1.0[,covar])
      capture.output(file="DivAndRichnessStats.txt",myttest,append=T)
    }else{
      myanova<-aov(er.x.1.0[,measure]~er.x.1.0[,covar])
      capture.output(file="DivAndRichnessStats.txt",myanova,append=T)
      capture.output(file="DivAndRichnessStats.txt",summary(myanova),append=T)
      capture.output(file="DivAndRichnessStats.txt",TukeyHSD(myanova),append=T)
    }
  }
}


#### 6 #### Using PERMANOVA testing (adonis function on bray-curtis distances) test which of the explanatory variables has a strongest link with microbiome structure?

x.6.0<-xB
wh0=genefilter_sample(x.6.0,filterfun_sample(function(x) x>1), A=0.01*nsamples(x.6.0))
x.6.0<-prune_taxa(wh0,x.6.0)
x.6.0.genus<-tax_glom(x.6.0,taxrank="Genus")
x.6.0<-x.6.0.genus
x.6.0 = transform_sample_counts(x.6.0, function(x) (((1* (x))/sum(x))))

# Non-parametric analysis of variance using PERMANOVA through adonis function

x.6.0.response<-t(otu_table(x.6.0))
metadata<-data.frame(sample_data(x.6.0))

x.6.0.explanatory<-metadata[,c("PSOE","HIVStatus","Site")]
for(var in colnames(x.6.0.explanatory)){
  explanatory=data.frame(x.6.0.explanatory[,eval(var)])
  myadonis<-adonis(x.6.0.response~.,data=explanatory)
  capture.output(paste("Adonis on: ",var),file="adonis_tests.txt",append=T)
  capture.output(myadonis,file="adonis_tests.txt",append=T)
}

adonis(x.6.0.response~.,data=x.6.0.explanatory)
x.6.0.adonis<-adonis(x.6.0.response~.,data=x.6.0.explanatory)
simper(x.6.0.response,x.6.0.explanatory[,"HIVStatus"])
simper(x.6.0.response,x.6.0.explanatory[,"PSOE"])
simper(x.6.0.response,x.6.0.explanatory[,"Site"])

#### 3 #### Create barplot for each sampling site, separately, showing fylum distributions. 

x.3.0<-xB
# We apply a more stringent filter, basically we are not interested in rare OTUs but on
# general main trends of taxonomical composition

wh0=genefilter_sample(x.3.0,filterfun_sample(function(x) x>5), A=0.1*nsamples(x.3.0))
x.3.0<-prune_taxa(wh0,x.3.0)

### tax_glom function agglomerates/collapses all OTU belonging to the same taxonomical level

x.3.0.phylum<-tax_glom(x.3.0,taxrank="Phylum")

x.3.0.genus<-tax_glom(x.3.0,taxrank="Genus")

### Convert to data.frame for easier manipulation
ps.melt.x.3.0.phylum<-psmelt(x.3.0.phylum)
ps.melt.x.3.0.genus<-psmelt(x.3.0.genus)

### Let's look at the phylum data.frame
summary(ps.melt.x.3.0.phylum)


### We will calculate relative abundances as proportions for further analysis
ps.melt.x.3.0.phylum$AbundanceProportion <- ave(ps.melt.x.3.0.phylum$Abundance,list(ps.melt.x.3.0.phylum[,"SampleID"]), FUN=function(L) L/sum(L))
ps.melt.x.3.0.genus$AbundanceProportion <- ave(ps.melt.x.3.0.genus$Abundance,list(ps.melt.x.3.0.genus[,"SampleID"]), FUN=function(L) L/sum(L))

reverse.levels <- function(x) {
  if(is.factor(x)) {
    x <- factor(as.character(x), levels=rev(levels(x)), ordered=TRUE)
  } else if(is.data.frame(x)) {
    for(i in seq_along(x)) {
      if(is.factor(x[,i])) {
        x[,i] <- factor(as.character(x[,i]), levels=rev(levels(x[,i])), ordered=TRUE)
      } else {
        warning(paste0('Column ', i, ' is not a factor.'))
      }
    }
  } else {
    stop(paste0('Unsupported format: ', class(x)))
  }
  return(x)
}

### Barplots at phylum. We use decreasing phylum abundances as X order 
### and the phylum abundance of first sample as y order
### This careful use of X/Y order will already reveal patterns
### Define Level Order for X axis (SampleID)
my.levels=orderBy(~-AbundanceProportion,data=ps.melt.x.3.0.phylum[ps.melt.x.3.0.phylum$Phylum=="Firmicutes",])$SampleID
ps.melt.x.3.0.phylum$SampleID<-factor(ps.melt.x.3.0.phylum$SampleID,
                                      levels=my.levels,ordered=T)


### Define Level Order for Y axis (phylum,genus,species)
my.levels<-orderBy(~-AbundanceProportion,data=ps.melt.x.3.0.phylum[ps.melt.x.3.0.phylum$SampleID == as.character(my.levels[1]),])$Phylum
ps.melt.x.3.0.phylum$Phylum<-factor(ps.melt.x.3.0.phylum$Phylum,levels=my.levels,ordered=T)

### We also take care of colors
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorCount=length(unique(ps.melt.x.3.0.phylum$Phylum))
getPalette=colorRampPalette(brewer.pal(12,"Set3"))
p<-ggplot(ps.melt.x.3.0.phylum,aes(x=SampleID,y=AbundanceProportion,order=ps.melt.x.3.0.phylum$Phylum))
p+geom_bar(stat="identity",aes(fill=Phylum))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(title="Phylum Level")+
  scale_fill_manual(values=getPalette(colorCount))+ facet_grid(~Site, scales = "free")


### Let's save the plot
ggsave("x.3.0.phylum.barplot.ordered.pdf")



#### 4 #### Run ordination analysis using PCoA/Wunifrac distances on the dataset and map sampling site and HIV Status as color/shape.


x.4.0<-xB
#We keep OTUs that appear at least in two different samples
wh0<-genefilter_sample(x.4.0,filterfun_sample(function(x) x>0), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)


### We transform the data to proportion abundances
x.4.0 <- transform_sample_counts(x.4.0, function(x) ((x/sum(x))))


### MDS Ordination 
x.4.0.ord<-ordinate(x.4.0,"MDS",distance = "bray",trymax=200)


### Simple Plot
p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples

### Let's map some metadata with coloured and ellipses

p.4.0.samples  + geom_point(aes(color=Site,shape=HIVStatus,size=PSOE)) +ggtitle(" MDS(Bray-Curtis)")+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(legend.text=element_text(size=12,face="italic"),legend.title=element_text(size=12))+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=Site),level=0.95)+
  theme(plot.title=element_text(lineheight=1,face="bold",size=15))

## WUNIFRAC

### Let's do a similar analysis but using phylogenetic OTU-OTU relationship
### using unifrac distances
x.4.0.ord<-ordinate(x.4.0,"MDS",distance="wunifrac",trymax=200)


### Simple Plot
p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples

### Now with metadata mapping

p.4.0.samples  + geom_point(aes(color=Site,shape=HIVStatus,size=PSOE)) +ggtitle(" MDS(WUnifrac)")+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(legend.text=element_text(size=12,face="italic"),legend.title=element_text(size=12))+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=Site),level=0.95)+
  theme(plot.title=element_text(lineheight=1,face="bold",size=15))

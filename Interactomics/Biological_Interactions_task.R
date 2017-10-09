
setwd('C:/Users/alber/Alberto/M�ster Omics/Interactomics/JPlanas/ftask/')
library(igraph)
library(data.table)

#intact<- fread('intact.txt', header =TRUE, data.table = FALSE)
# store intact object

#save(intact, file = 'intactobject.Rdata')

load('intactobject.Rdata')
# know and select two pathogens
colnames(intact)
head(intact$`Taxid interactor A`)
intact[,10][1:50]

species<- unique(intact[,10])
species[1:100]

positionviruses<- grep("virus",species )
speciesofvirus<- species[positionviruses]
speciesofvirus

# species selected

#"taxid:11065(den2n)|taxid:11065(\"Dengue virus type 2 (strain Thailand/NGS-C/1944) (DENV-2)\")"

#"taxid:128952(ebozm)|taxid:128952(\"Zaire ebolavirus (strain Mayinga-76)\")" 

# looking for interactions Human - Pathogen and Pathogen - Human  of Dengue virus type 2

HDenpos<- which(intact[,10] =="taxid:9606(human)|taxid:9606(Homo sapiens)" & intact[,11] == "taxid:11065(den2n)|taxid:11065(\"Dengue virus type 2 (strain Thailand/NGS-C/1944) (DENV-2)\")" )
HDenpos_rev<- which(intact[,11] =="taxid:9606(human)|taxid:9606(Homo sapiens)" & intact[,10] =="taxid:11065(den2n)|taxid:11065(\"Dengue virus type 2 (strain Thailand/NGS-C/1944) (DENV-2)\")" )

Human_Dengue<- intact[HDenpos,]
  
Dengue_Human<- intact[HDenpos_rev,]

# getting only protein IDs

HD<- Human_Dengue[,c(1,2)]
DH<- Dengue_Human[,c(1,2)]
head(DH)
head(HD)

# And now merging in to obtain only one data.frame by pathogen type

HinterwDengue<- rbind(DH,HD)
dim(HinterwDengue)


# Repeat process with Zaire ebolavirus


Hebpos<- which(intact[,10] =="taxid:9606(human)|taxid:9606(Homo sapiens)" &intact[,11] == "taxid:128952(ebozm)|taxid:128952(\"Zaire ebolavirus (strain Mayinga-76)\")" )

Heb_rev<- which(intact[,11] =="taxid:9606(human)|taxid:9606(Homo sapiens)" &intact[,10] =="taxid:128952(ebozm)|taxid:128952(\"Zaire ebolavirus (strain Mayinga-76)\")" )

Human_zeb<- intact[Hebpos,]

zeb_Human<- intact[Heb_rev,]

# getting only protein IDs

Hzeb<- Human_zeb[,c(1,2)]
zebH<- zeb_Human[,c(1,2)]
head(Hzeb)
head(zebH)

# And now merging in to obtain only one data.frame by pathogen type

Hinterwzebola<- rbind(Hzeb,zebH)
dim(Hinterwzebola)

# take a look of both data.frames

head(Hinterwzebola)
head(HinterwDengue)

# build the network

# graphically representing Human and Dengue interactions
humden<-graph.data.frame(HinterwDengue,directed=F)
#jpeg("hd.jpg",bg="black")

plot(humden,vertex.label.color ='white',vertex.color='green',vertex.label.dist=1)

#dev.off()
# graphically representing Human and Zaire ebolavirus


humzeb<-graph.data.frame(Hinterwzebola,directed=F)

#jpeg("he.jpg", bg = "black")
plot(humzeb,vertex.label.color ='transparent',vertex.color='blue',vertex.label.dist=1)
#dev.off()

####################################
#### Studying Dengue type 2 case#### 


#n� of nodes
vcount(humden)
#Names of nodes
unique(V(humden))
#n� of edges
ecount(humden)
#n� of components
components(humden)$no
#size of components
components(humden)$csize
#degree

degree(humden)

# degree distribution of the network

hist(degree(humden))

#eigenvalue

evc1<-evcent(humden)
evc1[[1]]

#betweeness
betweenness(humden)

#closeness

closeness(humden, mode = "all")
# proteins with highest degree in human-pathogen network by organism

sortedprots<-sort(degree(humden),decreasing = T)

protshum<- names(sortedprots)[names(sortedprots)%in%unique(Human_Dengue[,1])]
protshum2<-names(sortedprots)[names(sortedprots)%in%unique(Dengue_Human[,2])]#les que s�n humanes

Human_protein_names<- unique(c(protshum,protshum2))

sortedprots[Human_protein_names] #human proteins sorted by degree


dengueprots<- names(sortedprots)[names(sortedprots)%in%Dengue_Human[,1]]#les que s�n de patogen
dengueprots2<-names(sortedprots)[names(sortedprots)%in%unique(Human_Dengue[,2])] #patogen tb

Dengue_protein_names<- unique(c(dengueprots,dengueprots2))

sortedprots[Dengue_protein_names] #pathogen proteins sorted by degree




#######################################
#### Studying Zaire Ebolavirus case#### 


#n� of nodes
vcount(humzeb)
#Names of nodes
unique(V(humzeb))
#n� of edges
ecount(humzeb)
#n� of components
components(humzeb)$no
#size of components
components(humzeb)$csize
#degree

degree(humzeb)

# degree distribution of the network

hist(degree(humzeb))

#eigenvalue

evc2<-evcent(humzeb)
evc2[[1]]

#betweeness
betweenness(humzeb)

#closeness

closeness(humzeb)
# proteins with highest degree in human-pathogen network by organism

sortedprots2<-sort(degree(humzeb),decreasing = T)

protshumz<- names(sortedprots2)[names(sortedprots2)%in%unique(Human_zeb[,1])]
protshumz2<-names(sortedprots2)[names(sortedprots2)%in%unique(zeb_Human[,2])]#les que s�n humanes

Human_protein_namesZ<- unique(c(protshumz,protshumz2))

sortedprots2[Human_protein_namesZ] #human proteins sorted by degree


sortedprots2[Human_protein_namesZ][1:10] #top ten human proteins


ebolaprots<- names(sortedprots2)[names(sortedprots2)%in%zeb_Human[,1]]#les que s�n de patogen
ebolaprots2<-names(sortedprots2)[names(sortedprots2)%in%unique(Human_zeb[,2])] #patogen tb

Ebola_protein_names<- unique(c(ebolaprots,ebolaprots2))

sortedprots2[Ebola_protein_names] #ebola proteins sorted by degree



## merge human proteins that interact with both viruses

human_proteins_interwith_pathogens<- unique(c(Human_protein_namesZ,Human_protein_names))

###############################
## build human-human network ##

hh<- which(intact[,10] =="taxid:9606(human)|taxid:9606(Homo sapiens)" & intact[,11] == "taxid:9606(human)|taxid:9606(Homo sapiens)" )
Hum_Hum<- intact[hh,] 
HH<- Hum_Hum[,c(1,2)]

huhu<-graph.data.frame(HH,directed=F)

#knowing data distribution in degree and betweenness

hist(degree(huhu))
hist(degree(huhu)[human_proteins_interwith_pathogens])
shapiro.test(degree(huhu)[runif(5000,1,length(degree(huhu)))]) #random selection due to test sample size limitation (max 5000)
shapiro.test(degree(huhu)[human_proteins_interwith_pathogens]) 

hist(betweenness(huhu)[human_proteins_interwith_pathogens])  
hist(betweenness(huhu))
shapiro.test(betweenness(huhu)[human_proteins_interwith_pathogens]) 
shapiro.test(betweenness(huhu)[runif(5000,1,length(degree(huhu)))]) 

#After realising that data do not follow a normal distribution we run the Wilcoxon test to check centrality
# measures between human proteins interacting with pathogens and human proteins which not interact

wilcox.test(degree(huhu),degree(huhu)[human_proteins_interwith_pathogens], alternative = "two.sided") # found significant differences
wilcox.test(betweenness(huhu),betweenness(huhu)[human_proteins_interwith_pathogens],alternative = "two.sided") # found significant differences

## checking data distribution in closeness
hist(closeness(huhu))
hist(closeness(huhu)[human_proteins_interwith_pathogens])
shapiro.test(closeness(huhu)) 
shapiro.test(closeness(huhu)[human_proteins_interwith_pathogens]) 

#Again, data do not follow a normal distribution so we run Wilcoxon test

wilcox.test(closeness(huhu),closeness(huhu)[human_proteins_interwith_pathogens], alternative = 'two.sided') # found significant differences



#Task 2: Bladder.6.

#Remove all objects 
rm(list=ls())

#Set the working directory
setwd("C:/Users/alber/Alberto/Màster Omics/Bioinformatics/R/Malu Calle/Task2")

#Read the data Bladder.6.txt

bladder6<- read.table("Bladder.6.txt", header = T, sep = " ")

#Look at the dimension and the first 6 rows of the data and do a summary.

head(bladder6)
dim(bladder6)
summary(bladder6)

#Using the function attach we can call the variables only by it's name

attach(bladder6)

#1.	Convert variables "y" and "gender" to factor variables.
  #Use the function factor(x) in orther to tell R to treat the variables as categorical, not continuous.
#2.	Define the labels for variables "y" and "gender".

concase<-factor(y,levels=c(0,1),labels=c("control","case")) #Specify the names of the categories control and case

gender<-factor(gender,levels=c(1,2), labels=c("male","female")) #Specify the names of the categories male and female

summary(concase)# Do a summary of concase to know the number of cases and controls
summary(gender)# Do a summary of gender to know the number of males and females

#3.	Build a frequency table for "y" that contains, both, the absolute and the relative frequencies. The same for "gender".

options(digits=4) #Set the decimal numbers we want
freq.concas <- table(concase) #table of absolute frequencys of control&case

relfrq.concas<-prop.table(freq.concas) #table of relative frequencys of control&case

freqtable_cc<- cbind(freq.concas,relfrq.concas) #complete table of frequencys of control&case
freqtable_cc

freq.gender<- table(gender) #table of absolute frequencys of gender

relfrq.gender<-prop.table(freq.gender) #table of relative frequencys of gender

freqtable_gender<-cbind(freq.gender,relfrq.gender) #complete table of frequencys of gender
freqtable_gender

#4.	Plot gene1 levels as a function of "gender".

boxplot(gene1~gender, ylab="gene1 levels") #boxplot of gene1 levels as a function of gender

#5.	Test the normality of the expression levels of the 20 genes (use function apply).
#How many genes are not normally distributed and which are these genes?

normtest<-apply(bladder6[,4:23],2,shapiro.test) #apply the shapiro test to see which gene expresion levels follow a normal distribution
normtest #a list of the 20 genes with test applied 

 #Generate a for loop in order to store all the p.values in an empty vector to ask further questions
p.values<-c(NULL)
?length
for(i in 1:length(normtest)){
  x<-normtest[[i]]$p.value
  p.values<- c(p.values, x)
  }
?sum
sum(p.values<0.05) #With sum and less than 0.05 we ask R how many genes are not normally distributed

which(p.values <0.05) #With which and less than 0.05 we ask R which are those genes that are not normally distributed

#6.	Test whether mean expression levels of gene1 and gene2 are equal.
#Note: this is a test for equality of means with paired samples.
#This test is performed by computing the difference of the two variables (gene1-gene2) 
#and testing whether the mean of the difference is equal to zero.

Wilcoxontest<-wilcox.test(gene1,gene2, paired=T) #We do the Wilcoxontest because levels of expression of gene 2 do not follow a normal distribution

Wilcoxontest$p.value #We refuse H1, we do not reject the null, but Ho is not proved.The mean expression levels of gene1 and 2 are not significantly different

#7.	Test if the mean expression levels of gene1 are equal between cases and controls
  #Test for the equality of two variances among levels of gene1 and control/case to know which test should be used later

var.test(gene1~concase) 

  # We find that variances are not equal because p.value is lower than 0.05 so we have to apply the t.test function
  # but with the argument var.equal=F.

meanexpcc<-t.test(gene1~concase,var.equal=F)
meanexpcc$p.value #With this p.value we cannot accept H1 so the means are not significantly different

#8.	Obtain a 95% bootstrap confidence interval for the 20th percentile of gene1. 
?quantile

library(boot) #load the library boot
?boot #ask fot the arguments of function boot

#Generate our statistic of interest in our case the 20th percentile
bootstrp<-function(data,i) 
{
  q<-quantile(data[i],probs = 0.2)
  return(q)
}

#apply the function boot and boot.ci to our data in order to obtain the C.I 
boot.out<-boot(gene1, bootstrp, R=1000)
boot.ci(boot.out, conf = 0.95, type = "perc")

#9.	Contrast using a permutation test whether the 20th percentile of gene1 for male and female are equal.
#Permutation test for the equality of the 20th percentil of gene1
#between male and female

# The summary statistic is Q1-Q2

# Function that computes our statistic of interest

S.function<-function(x){
  x1<-x[gender=="male",]
  x2<-x[gender=="female",]
  quant.male<-quantile(x1$gene1,probs = 0.2)
  quant.female<-quantile(x2$gene1,probs=0.2)
  S<-quant.male-quant.female
  return(S)
}

# Appply this function to our data, bladder6
Sobs<-S.function(bladder6)
Sobs

#Obtain the null distribution of the statistic qdif under the null hypothesis H0: Q1=Q2. 
#This can be achieved by permuting the gender variable in order to break any 
#relationship between gene1 and gender 

S.null<-function(x){
  sgender<- sample(gender)
  x1<-x[sgender=="male",]
  x2<-x[sgender=="female",]
  quant.male<-quantile(x1$gene1,probs = 0.2)
  quant.female<-quantile(x2$gene1,probs=0.2)
  S<-quant.male-quant.female
  return(S)
}

# Running 1000 permutations we obtain the null permutational distribution of the test statistic S
?replicate
S.nullpermutdist <- replicate(1000,S.null(bladder6))
S.nullpermutdist[1:3]

# histogram of the permutational null distribution created by our data with mean near 0. 
hist(S.nullpermutdist,xlim = c(-0.5,0.5), ylim=c(0,350))

#Compute the permutation test pval  for the one-tailed alternative H1: Q1>Q2:
# p-value=P(|S.nullpermutdist|>|Sobs|)
#look for absolute values in order to not ignore a part of  my values
p.value<-(sum(abs(S.nullpermutdist)>abs(Sobs))+1)/(length(S.nullpermutdist)+1) #count how many numbers are to the right and compare it to all numbers in the null rates
p.value #We do not reject the null, we refuse H1, the are no significant differences between booth percentiles                                                             

#10.	Perform a nonparametric test for association of gender and the risk of disease. 
#Provide the OR (change the levels of "gender" if necessary in order that the given OR is larger than 1)

#install package epitool to be able to use the function oddsratio
install.packages("epitools")
library("epitools")

table.or<-table(gender,concase) #Perform a 2x2 table of gender and control/case

?oddsratio #ask for the arguments of function oddsratio

oddsratio(table.or)

Or1<-oddsratio(table.or, rev="c")#change the columns to have cases as column of reference 
Or1$measure #we change columns and rows in order to give an OR>1

Or2<-oddsratio(table.or, rev="both")
Or2$measure #here we see that males have 1.3 more risk of disease, even though, the difference is not signifcant
            #because the odds ratio is between (0.87 - 1.83) and includes 1.

#11.	Explore for possible relationship between methylation and gene expression.

gene_exp<- gene1+gene2+gene3+gene4+gene5+gene6+gene7+gene8+gene9+gene10+gene11+gene12+gene13+gene14+gene15+gene16+gene17+gene18+gene19+gene20
 
 #Create an object that is the result of adding all the gene expression in one vector that means all 
 # the gene expression by individual with the aim of searching relation with methylation

shapiro.test(gene_exp) #Look if data follows a normal distribution to decide about further testing
                       # We see that certainly the p-value>0.05, so data follow a Normal distribution

shapiro.test(methyl) #Notice that data do not follow a normal distribution. 
                     #It is important in the next test to realise if there's any relation with gene expression


cor.test(methyl,gene_exp,method = c("spearman")) # We see that p.value<0.05, that means that the test is significant for us
                                                 # and rho is near 0 so we talk about very low correlation.      
plot(methyl,gene_exp) # We check that certainly very light correlation does exist.

#12.	Identify genes that are related to the risk of bladder cancer using a multivariate logistic regression model with stepwise variable selection. 
#Denote the selected model as "best.model". Interpret the obtained model.

model<-glm(concase~gene1+gene2+gene3+gene4+gene5+gene6+gene7+gene8+gene9+gene10+gene11+gene12+gene13+gene14+gene15+gene16+gene17+gene18+gene19+gene20, family=binomial()) #full model
step(model,direction="both")#stepwise variable selection
best.model<-glm(concase ~ gene2 + gene4 + gene6 + gene18 + gene20,family = binomial()) #best model
summary(best.model) #the interaction of these genes (2,4,6,18,20) explains better the risk of bladder cancer
                   
#13.	Analyze the classification ability of "best.model" (ROC curve and AUC) according to the following schemes:
   #a.	Apparent validation of "best.model" using the same data that was used for model building.
library(ROCR)
#classification accuracy best.model
lp<-best.model$linear.predictors #Store all the linear predictors of our best.model
lp[1:10]
pred <- prediction(lp, concase)
perf <- performance(pred, "tpr", "fpr" ) #make the prediction and perform the ROC curve

plot(perf) #We see that the curve is below 0.5 and we do the complementary

pred <- prediction(1-lp, concase)
perf <- performance(pred, "tpr", "fpr" )

plot(perf)#now ROC curve is a matter of interest
title("ROC curve for best.model:concase~gene 2,4,6,18 and 20; Bladder6")

performance(pred,"auc")  #calculate the area under the curve
  
   #b.	Crosvalidation with k = 5 for "best.model".


# breakpoints for the different cv sets
?ceiling
K <- 5
n <- nrow(bladder6)  
bp <- c(0, ceiling((1:K)*n/K))  #To get the integers of our 4 breakpoints that create 5 parts of the sample. 

# index permutation

?sample
index <- sample(rownames(bladder6))

# CV for k=5 we built 5 models and test it with 5 different parts of the sample each time
pred <- NULL

for(i in 1:K){
  indTest <- index[(bp[i]+1):bp[i+1]] 
  indTrain <- index[!(index %in% indTest)]  
  model.i <- glm(y~gene2+gene4+gene6+gene18+gene20, data=bladder6[indTrain,], family=binomial())
  pred.i <- predict(model.i, newdata=bladder6[indTest, ])
  pred <- c(pred,pred.i)
}  


# pred contains the concatenated linear predictor of the permuted individuals
# pred.y contains the ordered linear predictions

pred.y <- pred[rownames(bladder6)] 

#Once we have the linear predictors in a vector we can perform the ROC curve 

prediction <- prediction(pred.y, concase) 
perf <- performance(prediction, "tpr", "fpr" )
plot(perf)

prediction_c <- prediction(1-pred.y, concase)#do the complementary curve 
perf_c <- performance(prediction_c, "tpr", "fpr" )
plot(perf_c)
title("ROC curve for best.model.Cross-Validation with K=5")
performance(prediction_c,"auc") #calculate the area under the curve
   #c.	Explain why both the "a" and "b" overestimate the classification accuracy of "best.model".
   #Suggest a scheme for determining an unbiased estimate of the classification accuracy of "best.model".

    #Both models overestimate the classification accuracy because somehow data are used to build
    #the model and to test it. The ideal validation is on a new data set. We look for an independent sample
    #to built our model and we test it in best.model.


#14.	For each variable "genei", i=1, ., 20, define a new factor variable "levels.genei" with 
#values "high" if gene levels are positive or zero and "low" if gene levels are negative. 

levels.gene1<- factor(gene1<0, labels=c ("high", "low"))
levels.gene2<- factor(gene2<0, labels=c ("high", "low"))
levels.gene3<- factor(gene3<0, labels=c ("high", "low"))
levels.gene4<- factor(gene4<0, labels=c ("high", "low"))
levels.gene5<- factor(gene5<0, labels=c ("high", "low"))
levels.gene6<- factor(gene6<0, labels=c ("high", "low"))
levels.gene7<- factor(gene7<0, labels=c ("high", "low"))
levels.gene8<- factor(gene8<0, labels=c ("high", "low"))
levels.gene9<- factor(gene9<0, labels=c ("high", "low"))
levels.gene10<- factor(gene10<0, labels=c ("high", "low"))
levels.gene11<- factor(gene11<0, labels=c ("high", "low"))
levels.gene12<- factor(gene12<0, labels=c ("high", "low"))
levels.gene13<- factor(gene13<0, labels=c ("high", "low"))
levels.gene14<- factor(gene14<0, labels=c ("high", "low"))
levels.gene15<- factor(gene15<0, labels=c ("high", "low"))
levels.gene16<- factor(gene16<0, labels=c ("high", "low"))
levels.gene17<- factor(gene17<0, labels=c ("high", "low"))
levels.gene18<- factor(gene18<0, labels=c ("high", "low"))
levels.gene19<- factor(gene19<0, labels=c ("high", "low"))
levels.gene20<- factor(gene20<0, labels=c ("high", "low"))

#15.	Perform a nonpametric test of association between each variable "levels.genei" and the risk of disease. 
#Adjust the p-values for multiple testing according to an fdr threshold equal to 0.1. Interpret the results.

chistest_1<- chisq.test(table(levels.gene1,concase))
chistest_1$p.value
chistest_2<- chisq.test(table(levels.gene2,concase))
chistest_2$p.value
chistest_3<- chisq.test(table(levels.gene3,concase))
chistest_3$p.value
chistest_4<- chisq.test(table(levels.gene4,concase))
chistest_4$p.value
chistest_5<- chisq.test(table(levels.gene5,concase))
chistest_5$p.value
chistest_6<- chisq.test(table(levels.gene6,concase))
chistest_6$p.value
chistest_7<- chisq.test(table(levels.gene7,concase))
chistest_7$p.value
chistest_8<- chisq.test(table(levels.gene8,concase))
chistest_8$p.value
chistest_9<- chisq.test(table(levels.gene9,concase))
chistest_9$p.value
chistest_10<- chisq.test(table(levels.gene10,concase))
chistest_10$p.value
chistest_11<- chisq.test(table(levels.gene11,concase))
chistest_11$p.value
chistest_12<- chisq.test(table(levels.gene12,concase))
chistest_12$p.value
chistest_13<- chisq.test(table(levels.gene13,concase))
chistest_13$p.value
chistest_14<- chisq.test(table(levels.gene14,concase))
chistest_14$p.value
chistest_15<- chisq.test(table(levels.gene15,concase))
chistest_15$p.value
chistest_16<- chisq.test(table(levels.gene16,concase))
chistest_16$p.value
chistest_17<- chisq.test(table(levels.gene17,concase))
chistest_17$p.value
chistest_18<- chisq.test(table(levels.gene18,concase))
chistest_18$p.value
chistest_19<- chisq.test(table(levels.gene19,concase))
chistest_19$p.value
chistest_20<- chisq.test(table(levels.gene20,concase))
chistest_20$p.value

#store all the p-values in a new vector

pvals_chistest<-c(chistest_1$p.value,chistest_2$p.value,chistest_3$p.value,chistest_4$p.value,chistest_5$p.value,
                  chistest_6$p.value,chistest_7$p.value,chistest_8$p.value,chistest_9$p.value,chistest_10$p.value,
                  chistest_11$p.value,chistest_12$p.value,chistest_13$p.value,chistest_14$p.value,chistest_15$p.value,
                  chistest_16$p.value,chistest_17$p.value,chistest_18$p.value,chistest_19$p.value,chistest_20$p.value)

#we adjust p.values with p.adjust function

?p.adjust

adjusted.pvalues<-p.adjust(pvals_chistest,method="fdr")

#we set the threshold at 0.1

adjusted.pvalues<0.1

which(adjusted.pvalues<0.1) #To know which genes are below our threshold

#16.	Using the last 30 variables corresponding to gene expression levels in three different pathways, 
#perform a clustering analysis (hierarchical and k-means) and explore groups of genes that have similar expression.
#Explain the results. 

#hierarchical
plot(hclust(dist(t(bladder6[,25:54]), method="euclidian"),method="average")) #clear representation of the clusters

cutree(hclust(dist(t(bladder6[,25:54]), method="euclidian"),method="average"),3)

#k-means

kmeans2<-kmeans(t(bladder6[,25:54]),3)
kmeans2
ls(kmeans2)
kmeans2$cluster

#to see in a more acurate way the distrubution of pathways in each cluster
sort(kmeans2$cluster)



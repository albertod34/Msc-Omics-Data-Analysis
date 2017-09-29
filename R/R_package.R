

setwd("C:/Users/alber/Alberto/M�ster Omics/Bioinformatics/R/JL.Mosquera/finaltask")

dat<- read.table('microRNA_expressions.txt', header = T)
class(dat)
head(dat)
class(dat$group)
class(dat$expression)

#Function 1

welstat<- function(d)
{
  #control input conditions 
  stopifnot(is.factor(d[,1]),is.numeric(d[,2]) )
  
  levels(d[,1])<- c(0,1)
  
  firstgroup<- subset(d, d[,1] == "0")
  secgroup<- subset(d, d[,1] == "1")
  
  out<-t.test(firstgroup[,2],secgroup[,2], var.equal=FALSE, paired=FALSE)
  
  tstat<- out[["statistic"]]
  
  return(round(tstat,3))
}




#Function 2

permutedWS<- function(d)
{
  sampdf<- apply(d,2,sample)
  ndf<- data.frame(sampdf[,1],d[,2])
  out<- welstat(ndf)
  
  return(out)
}


#Function 3 

tvpval<- function(d, val = 1000)
{
  if(val<1000) warning('It is recommended a minimum of 1000 permutations')
  
  out<-welstat(d)
  
  nd<- replicate(val, permutedWS(d))
  
  pval<- (sum(nd>out) + 1)/(val+1)
  
  flist<- list(out, nd, pval)
  names(flist)<- c('Observed t-value',paste(val,'t-values vector'), 'pvalue')
  
  return(flist)
}


package.skeleton(name = "Alberto_Dominguez",
                 list = c("welstat", "permutedWS", "tvpval", "dat"))


rm.man.files <- c("Alberto_Dominguez-package.Rd", "tvpval.Rd", "welstat.Rd",
                  "permutedWS.Rd","dat.Rd")



file.remove(file.path(file.path(getwd(), "Alberto_Dominguez/man"), rm.man.files))
 

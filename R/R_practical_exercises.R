#############################################  Practical exercises #################################################################
####################################################################################################################################


#Exercise 1. Load the following function in your R workspace


generateExpressions <- function(nrow = 5, ncol.1 = 1, mean.1 = 0, sd.1 = 1, samples.2 = F, ncol.2 = NULL, mean.2 = NULL, sd.2 = NULL)
{
	# Creating the main object, a matrix with especified rows and columns filled with random numbers from a normal distribution of specified 
	# mean and sd. 
	
	out<-matrix(rnorm(nrow*ncol.1,mean.1,sd.1),nrow,ncol=ncol.1) # generating the matrix
	
	if(samples.2) #managing parameter sample.2 
	{
	  if(is.null(ncol.2)) ncol.2 <- ncol.1 # if it's null the parameter ncol.2, takes the same value of ncol.1
	  if(is.null(mean.2)) mean.2 <- mean.1 #the same as above but with the mean.2
	  if(is.null(sd.2)) sd.2 <- sd.1 # the same with sd.2
	
	    matrix.2 <- matrix(rnorm(nrow * ncol.2, mean.2, sd.2) , nrow, ncol = ncol.2) #how to create matrix.2
	    out <- cbind(out, matrix.2) # adding matrix.2 to previous created matrix by columns with cbind() 
	
	}else{ # the other situations not especified before
	  ncol.2<- 0

	}
	
	# Creating row and column names
	
	colnames(out)<-paste("sample",1:(ncol.1+ncol.2),sep=".") 
	rownames(out)<-paste("gene",1:nrow,sep=".")
	
	return(out) # what we want out of the local environment 
}
	#1. Invoke generateExpressions with different input parameters.
	
	
	##  generateExpressions(3,2,4,2,TRUE,4,3,2)
	##         sample.1    sample.2 sample.3 sample.4  sample.5  sample.6
	##  gene.1 2.446906  5.12876112 2.888788 2.849143 0.9838071 0.7113081
	##  gene.2 7.460721  0.91760685 4.093362 1.055589 3.2272832 1.2616278
	##  gene.3 1.369406 -0.04079142 3.715019 2.925887 4.0801573 1.3033153
	
	##  generateExpressions(3,2,4,2,FALSE,4,3,2)
	##           sample.1 sample.2
	##  gene.1 0.03739991 7.996053
	##  gene.2 4.16236254 5.439820
	##  gene.3 1.78248595 3.524161
	
	
	#2. Rewrite the function in order to provide a more human readable structure. That is,
	#try to introduce the indentation, spaces and carriage-return linefeeds when they are
	#required as well as to write appropriate comments for clarifying the action that is being
	#performed.
	
	#3. Briefly describe what the function does and show an example.
	
	#The function generates a matrix of random numbers from a normal distribution with mean and standard deviation. It has eight arguments, all previously
	# specified, that means that you can run the function with no need to type any argument (i.e. generateExpressions()). With the arguments you can control
	#the size of the matrix (number of rows and columns) and the mean and standard deviation of the normal distribution. Besides, the argument samples.2, 
	#gives you the option to add another matrix with different or not mean and standard deviation values from a normal distribution too. In the example, we 
	# create two matrices binded by column with 3 rows and 5 (from matrix1) + 4 (from matrix2) = 9 columns. The first matrix has mean 4 and sd 2 and the second
	# one has mean 3 and sd 2.
	
	##  generateExpressions(3,5,4,2,TRUE,4,3,2)
	##         sample.1 sample.2 sample.3 sample.4 sample.5  sample.6  sample.7
	##  gene.1 5.360190 6.446443 3.039070 4.678023 6.326475 0.5313944 0.8510898
	##  gene.2 1.575867 2.779010 5.855871 2.347179 4.224280 2.2627575 2.9504904
	##  gene.3 3.035684 4.923540 2.249277 6.543539 1.646801 3.0722606 0.7003186
	##         sample.8 sample.9
	##  gene.1 6.122686 3.870344
	##  gene.2 4.164866 2.756472
	##  gene.3 3.611586 3.226111




##Exercise 2. File gene exp.data contains a subset of gene expression values from a human Affymetrix array. 
#Each probe consists of six values associated with three different experimental conditions (control, placebo and treatment).

#1. Input this file in your R workspace.
setwd('C:/Users/alber/Alberto/M�ster Omics/Bioinformatics/R/JL.Mosquera/finaltask/DataforExers1')
gexp<-read.delim('gene_expr.data')
class(gexp)
[1] "data.frame"
#2. How many probes are in the dataset?
head(gexp)
length(table(gexp$Probe))
[1] 33619	
#3. How many genes are in the dataset?
length(table(gexp$Symbol))
[1] 15629
#4. How many probes are in the dataset per each gene?
33619/15629
[1] 2.151065
#5. Write a function that reads an arbitrary file like gene exp.data and gives you these three numbers.

PrGenFunction <- function(gexpdata)
{
	gdata<- read.delim(gexpdata) # reading data
	
	#Control table parameters
	
	if(any(colnames(gdata) == "Probe") == FALSE) stop("The probes column must have the name Probe")
	if(any(colnames(gdata) == "Symbol") == FALSE) stop("The genes column must have the name Symbol")
	
	# Compute calculations
	n_probes <- length(table(gdata$Probe))
	n_genes<- length(table(gdata$Symbol))
	ppg<- n_probes/n_genes
	
	object <- cbind(n_probes,n_genes,ppg)
	colnames(object)<- c("Probes","Genes", "Probes per gene" )
	rownames(object)<- "Information"
	return(object)

}



##Exercise 3. Given a numeric vector, write a function called descriptive, that yields a
##new vector containing the following descriptive statistics:


# 1. n: the number of measurements
# 2. missing: the number of NA's
# 3. Mean: the average
# 4. Std.Dev: the standard deviation
# 5. Min: the minimum value
# 6. Q1: the quantile 1 (i.e. percentile 25%)
# 7. Median: the median
# 8. Q3: the quantile 3 (i.e. percentile 75%)
# 9. Max: the maximum value


descriptive <- function(v)
{
	if(!is.numeric(v)) stop("The function requires a numeric vector") # Controlling the input
	
	n<- length(v)
	
	if(any(is.na(v))) #looking for missing values
	{
	  Missing<- length(ob<-split(v,is.na(v))[[2]])
    }else{                   
	  Missing<-0    
	}
	
	# Doing the rest of descriptive statistics
	
	Mean<- mean(v, na.rm = TRUE)
	Std.Dev<- sd(v, na.rm = TRUE)
	Min<- min(v, na.rm = TRUE)
	Q1<- quantile(v,  0.25, na.rm = TRUE)
	Median<-median(v, na.rm = TRUE)
	Q3<- quantile(v, 0.75, na.rm = TRUE)
	Max<- max(v, na.rm = TRUE)
	
	# creating the descriptive vector and setting the names
	
	out<- c(n, Missing, Mean, Std.Dev, Min, Q1, Median, Q3, Max) 
	
	names(out)<- c( "Number of measurements", "Number of NA's", "Mean", "Standard Deviation","Minimum", "Quantile 1", "Median","Quantile 3" ,"Maximum") 
	
	# What we want out
	return(round(out,2))
	
}




##Exercise 4. Consider the data.frame that you read in exercise 1. 
##Write a function called exprDescriptive that invokes descriptive function for describing the 
##expression values of either samples (columns) or probes (rows). 
##The result must be a matrix where on rows are the samples (or probes depending on the option
##selected by the user) and on columns are the statistics.


exprDescriptive<- function(gexpdata, colofrownames = NULL , colsnotnumeric = NULL, probes = TRUE)
{
	gdata<- read.delim2(gexpdata) # reading data
	rownames(gdata)<- gdata[,colofrownames] # changing the rownames
	gdata<- gdata[,-colsnotnumeric]# deleting the non numeric columns in order to compute the statistics
	
	# making the computation process by probe or by gene
	if(probes)
	{
	  gmat<- apply(gdata,1,descriptive)
	}else{
	  gmat<- apply(gdata,2,descriptive)
	}
	gmat<- t(gmat)
	
	return(gmat)
}




##Exercise 5. The GC-content is a characteristic of a specific fragment of DNA or RNA, or that of the whole genome of an organism. 
##It is expressed as the percentage of nitrogenous bases on a DNA molecule that are either guanine (G) or cytosine (C).
##That is, the GC-content is calculated as:

##GC-content = (G + C / A + T + G + C) * 100

##Write a function called gc that calculates the GC-content of a given DNA sequence.

gc <- function(dnaseq)
{  
  if(length(dnaseq) != 1) stop( "Enter the sequence all together") #Control input parameter for further operations
  
  dna<-tolower(dnaseq)
  dna<-strsplit(dna,NULL)
  vec<-dna[[1]]
  
  if(length(grep("a|c|g|t",vec))!= length(vec)) stop(" DNA sequence is required ") #Control the kind of sequence
  tb<-as.numeric(table(vec)) # count the bases and store it as a numeric vector with alphabetical order
  gc_cont<- round(((tb[2]+tb[3])/(tb[1]+tb[2]+tb[3]+tb[4]))*100,2) # Calculating GC content
  
  return(gc_cont)
}





##Exercise 6. Write a function called wc that invokes the function readLines for reading a text file and gives you the 
##number of lines of that file. You can use file dna_seqs.fasta for testing your function.

wc <- function(x)
{
  lines<-readLines(x)
  
  out<-length(lines)
  
  return(out)
}




##Exercise 7. Write a function called summarySingleFasta that reads a fasta file with a single sequence and yields 
##a vector with the following four elements:
  
  # 1. Description: the description provided in the first line of the fasta sequence
  # 2. Sequence: the DNA sequence
  # 3. Width: the length of the DNA sequence
  # 4. GC.content: the GC-content of the sequence

##You can use the file single_dna_seq.fasta for testing your function.

summarySingleFasta <- function(f)
{
  fil<- readLines(f) # reading the file
  
  #Creating the vector
  Description <- fil[1]
  Seq <- fil[2]
  Width <- length(strsplit(Seq,NULL)[[1]])
  GC.Content<- gc(Seq)
  
  #pasting all the results and giving its names 
  out <- c(Description, Seq , Width, GC.Content)
  names(out)<- c('Description', 'Sequence', 'The width', 'GC-content %')
  
  return(out)
}



##Exercise 8. Write a function called summaryMultipleFasta that reads a fasta file with
##multiple DNA sequences, whose description consists of two elements (a potential name of a
##gene and a potential name of a specie) separated by a single space, and yields a data.frame
##such that for each sequence (rows) provides the following information (columns):
	
	#1. Gene: the name of the potential gene
	#2. Specie: the name of the potential specie
	#3. Sequence: the DNA sequence
	#4. Width: the length of the DNA sequence
	#5. GC.content: the GC-content of the sequence

##You can use the file dna_sequence.fasta for testing your function.

summaryMultipleFasta <- function(f)
{
  fi<- readLines(f) #reading the file
  
  mask<- grep(">",fi) # separating the names from the sequences 
  seqnames<- fi[mask]
  seqs<- fi[-mask]
  as<-gsub(">", "",seqnames)
  
  seps<-strsplit(as, " ") #separating the potential name of the gene and the potential name of the specie
  mat<- t(sapply(seps,as.matrix)) #creating a the matrix with the separate names
  mat2<- cbind(mat,seqs) #merging it with the sequences
  width<- sapply(strsplit(seqs,NULL),length)# Calculating width
  GC.content<- sapply(seqs,gc) # the gc.content
  
  #Building the final matrix
  matf<- cbind(mat,seqs,width,GC.content)
  df <- as.data.frame(matf)
  
  rownames(df)<-NULL
  colnames(df)<- c('Gene','Specie', 'Sequence','Width', 'GC.Content')
  
  return(df)
}

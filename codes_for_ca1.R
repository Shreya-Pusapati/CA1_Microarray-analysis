#installing the packages
BiocManager::install("GEOquery")
BiocManager::install("oligo")
BiocManager::install("affy")
BiocManager::install("limma")
install.packages("ggplot2")
install.packages("pheatmap")
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("hgu133plus2.db") #annotation package for the Affymetrix microarrays

if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db") #Annotation of the human genome

#loading packages into R
library(GEOquery)
library(affy)
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)

#downloading the chosen GSE as a list of expression sets using getGEO
gse50737 <- getGEO("GSE50737", GSEMatrix=TRUE)

#checking for class of the downloaded GSE data
show(gse50737)
class(gse50737)
class(gse50737[[1]])
gse50737_eset <- gse50737[[1]] #to get expression set data
pData(Gse50737_eset)[,c("geo_accession","tissue:ch1","treatment:ch1")] #1 covariate modified= treatment:ch1

#Reading CEL data using the oligo package.
gse50737_r <- getGEOSuppFiles('GSE50737') #does not work; stops downloading at 23.5 MB
#or

setwd("C:/Users/P SHREYA/Desktop/CAr") #set working directory
untar("GSE50737_RAW.tar", exdir = "celdata2/") # use this command on downloaded raw.tar data
setwd("C:/Users/P SHREYA/Desktop/CAr/celdata2")
Gse50737_celdata <- oligo::read.celfiles(list.celfiles()) #cel files are read and stored in Gse50737_celdata.

#normalization by RMA
Gse50737_eset <- oligo::rma(Gse50737_celdata) # using oligo package # annotation by pd.hugene.2.0.st #fast & simple 
#or
#the following method can be used:
pd <- phenoData(Gse50737_eset) #store the phenotype data into pd
pData(Gse50737_eset) 
#in this experiment only 1 covariate is modified- treatment:ch1.

hist(Gse50737_celdata,lwd=2,xlab='log intensity', which='all',
 main="CEL file densities before quantile normalisation") #worked.

#Microarray Data processing with RMA
cel_normalised <- normalize(Gse50737_celdata,method='quantile',which='pm')
hist(cel_normalised,lwd=2,xlab='log intensity', which='all',
     main="CEL file densities after quantile normalisation") #worked with 'all', not with 'pm'

#Summarisation: from probe-level to feature-level data.
cel_summarised <- oligo::rma(cel_normalised,background=TRUE,normalize=TRUE)
hist(cel_summarised) #background and normalize given as false so as to allow normalisation and background correction to occur.

#Identifying differentially expressed lncRNA & mRNA using linear models  (for processed getGEO data)
design <- model.matrix( ~ 0 + gse50737_eset[['treatment:ch1']])
colnames(design) <- levels(as.factor(gse50737_eset[['treatment:ch1']]))
design

colnames(design)[1] <- "BE" #renamed columns 
colnames(design)[2] <- "BT"
colnames(design)[3] <- "C"

contrast_matrix1 <- makeContrasts(BE-C,BT-C,BE-BT, levels=design)
contrast_matrix1

fit <- lmFit(gse50737_eset,design)
fit2 <- contrasts.fit(fit,contrasts=contrast_matrix1)
fit2 <- eBayes(fit2)


#Extracting differentially expressed genes
topTable(fit2)
summary(decideTests(fit2,lfc=1))


#volcano plot
volcanoplot(fit2,coef=2)
GOIs <- topTable(fit2,number=Inf,p.value = 0.05,lfc=2)

volcanoplot(fit2, coef=2, main=sprintf("%d features pass our cutoffs",nrow(GOIs)))
points(GOIs[['logFC'>2]],-log10(GOIs[['P.Value']]),col='red') #works with lfc>2

#heatmaps
interested_genes <- gse50737_eset[rownames(GOIs),]
heatmap(exprs(interested_genes))

#for CEL data
#Identifying differentially expressed lncRNA & mRNA using linear models 
#make groups of the files
Groups <- c("benzene.tox","benzene.tox","benzene.tox","benzene.tox","benzene.exp","benzene.exp","benzene.exp","control","control","control")
design2 <- model.matrix(~ 0 + factor(Groups))
design2
#colnames design
colnames(design2) <- c("benzene.exp","benzene.tox","control")
design2
fitcel <- lmFit(Gse50737_eset, design2)
fitcel2 <- eBayes(fitcel)

contrast_matrix3 <- makeContrasts(benzene.tox-control,benzene.exp-control,benzene.tox-benzene.exp, levels=design2)
contrast_matrix3

fitcel3 <- contrasts.fit(fitcel2,contrasts=contrast_matrix3)
fitcel4 <- eBayes(fitcel3)
topTable(fitcel4)
summary(decideTests(fitcel4,lfc=1)) #lfc=1 means fold change >1= more expression; 1=no difference; <1= more expressed in control.

volcanoplot(fitcel4,coef=2)
GOI_eset <- topTable(fitcel4,number=Inf,p.value = 0.05,lfc=2) #specifying lfc and p value

volcanoplot(fitcel4, coef=2, main=sprintf("%d features pass our cutoffs",nrow(GOI_eset)))
points(GOI_eset[['logFC'=2]],-log10(GOI_eset[['P.Value']]),col='red') #only works with logFC >2.

interested_genes_eset <- Gse50737_eset[rownames(GOI_eset),]
heatmap(exprs(interested_genes_eset))
pheatmap(exprs(interested_genes_eset))





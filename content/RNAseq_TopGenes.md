# Script Name: TopGenes

## Overview ##

#1) Read Count Normalization
#2) Extraction of Top Differentially Expressed Genes
#3) Visualization (Phantasus)

## Require packages ##

library(dplyr) # data re-orgnization & manipulation
library(tidyr) # contains 'separate' function
library(AnnotationDbi) # gene symbols and descriptions
library(Mus.musculus) #  mouse reference 
'%!in%' <- function(x,y)!('%in%'(x,y)) #  'NOT IN' function
library(DESeq2) # read count normalization

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load and Prepare Data Set

setwd('path/to/count/data') # set working directory

countData<-read.table('counts_mapped.tsv',header=T) # import mapped, raw count data

str(countData) # check class of object and data types

# Create metadata dataframe, first column must match column headings of countData
colData<-data.frame(id=colnames(countData));colData # extract column names from countData

# expand sample info across several discrete columns
colData <-colData %>% # pipe the df to dplyr operation
  separate(id,into = c("infection","group","status","Sample"), sep = "_", remove=F) # pass string to separate function on '_' character
  
# sync column names of countData and id variable of colData
colnames(countData)<-sub('.*_', '', colnames(countData)) # removes all but the sample names
colData$id<-colData$Sample # replaces the id column of colData with the Sample column
colData<-colData[1:4] # subset colData to first four rows
colData$group<-gsub('[[:digit:]]+', '', colData$group) #remove redundant info from group columm of colData

#if controls samples will be included designate them as such for 'status' and 'infection' variables
colData[colData$status=="1"|colData$status=="2",]$status<-"CON"
colData[colData$infection=="Tc",]$infection<-"CON"


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Read count normalization (using variance stablizing transformation)

# create deseq object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=colData, 
                              design=~1,
                              tidy=F)

# variance stabilizing transformation (ignoring experimental design)
ddsVST<-varianceStabilizingTransformation(dds,blind=T)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# subset to the top N most variable genes

# calculate variance for each gene
rv <- rowVars(assay(ddsVST))

# sort by variance, and select the top 12k
select <- order(rv, decreasing = TRUE)[seq_len(min(12000, length(rv)))]

# create dataframe of top12k
top12k<-as.data.frame((assay(ddsVST)[select, ]))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# write count matrix and annotation files for import to Phantasus

# write normalized read counts to tab-delimited file called top12k
write.table(top12k,'top12k_esxh.tsv',sep='\t',quote = F,col.names = NA)

# write colData to tab-delimited file called COLannotation
write.table(colData,'COLannotation_esxh.tsv',sep='\t',quote=F,row.names=F) # write column metadata

# add gene symbols to countData dataframe
countData$symbol<-mapIds(Mus.musculus,
                          keys=as.character(row.names(countData)), 
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")

# create 'NA' string value for missing values in gene symbol variable
countData$symbol<-ifelse(is.na(countData$symbol),"NA",countData$symbol)

# write tab-delimited file of gene IDs and gene symbols caleld ROWannotation
write.table(cbind(row.names(countData),countData$symbol),'ROWannotation_esxh.tsv',sep='\t',quote=F,row.names=F)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# import files to https://artyomovlab.wustl.edu/phantasus/
# this is a free and easy to use web tool for visualizaing gene expression sample by sample
# imports can be saved for a period of time and shared by URL with collaborators that are offsite




```rscript
# Script Name: RNAseq_Differential_Expression


## OVERVIEW##

#1) Data preparation
#1) Data inspection and visualization
#2) Differential expression analysis (DESeq2)

## required packages ##

library(gprofiler2) # gene set profiling, also mapping gene IDs and descriptions
library(dplyr) # data orgnization and manipulation
library(tidyr) # contains separate function
library(AnnotationDbi) # mapping gene symbols and descriptions
library(Mus.musculus) # mouse reference 
'%!in%' <- function(x,y)!('%in%'(x,y)) #  'NOT IN' function
library(DESeq2) #  DE analysis
source('/Users/gideon/Box/Rscripts/GeneSet_function.R') # function for conducting ORA and GSEA analyses
library(ggplot2) # plotting
library(ggrepel) # plotting aesthetics
library(reshape2) # data reorganization
library(ComplexHeatmap) # plotting
library(WebGestaltR) # gene set profiling
library(gtools)
library(gridExtra)



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
# sample selection for visualization

#insert samples that will NOT be analyzed
fordeg<-c(
  "S1","S2",
  "S3","S4",
  "S37","S38","S39","S40","S41","S42","S43","S44","S45","S46","S47","S48",
  "S27","S28"
  #colData$id[colData$infection=='dH']
)

colData2<-colData[colData$id %!in% fordeg,] # removes selected samples from colData
countData2<-countData[,colnames(countData) %!in% fordeg] # removes selected samples from countData



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data Visualization (PCA)

names(colData2)[3]<-"Tcell"

# normalize read counts with DESeq object
dds <- DESeqDataSetFromMatrix(countData=countData2, 
                              colData=colData2, 
                              design=~Tcell+status,
                              tidy=F)

# normalize read counts
ddsVST<-varianceStabilizingTransformation(dds,blind=F)

# PCA of genes with high variance (top 10K in this case)
Pca<-plotPCA(ddsVST,intgroup=c("Tcell",'status','infection'),ntop=10000, returnData=T)

# plot first two components of the PCA analysis
ggplot(PCA1, aes(Pca,PC2,color=interaction(status,infection),shape=Tcell))+
  geom_point(size=3)+
  #geom_label_repel(aes(label=name))+
  xlab(paste0('PC1 (',substr(tail(unlist(attributes(PCA1),use.names=F),2)[-2],1,4),')'))+ # % variance explained by component 1
  ylab(paste0('PC2 (',substr(tail(unlist(attributes(PCA1),use.names=F),2)[-1],1,4),')'))+ # % variance explained by component 2
  ggtitle("Tcell Effect")+
  theme_bw(14)

CONCLUSION: PCA is helpful to identify outliers and factors that divide the data set



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sample Visualization (heatmap)

heat<-dist(t(assay(ddsVST)),upper=T,diag=T) # extract distance measures
heat<-as.matrix(heat) # convert to matrix # convert dataframe to matrix

#OPTINAL: next three lines can be used to arrange samples as desired (heatmap hierarchical clustering must be set to "null")
org<-as.character(mixedsort(colData2$id)) # use if sample names are mix of letters and numbers
org<-as.character(arrange(colData2,id)$id) # use if sample names are just numbers or letters
heat2<-heat[org,org]# arrange rownames and columns in same orientation as the 'org' variable

# heatmap row annotation
row_ha<-HeatmapAnnotation(Group=colData2$group,Status=colData2$status,
                          which='row')

# PLOT HEATMAP (there is a parameter for adjusting dendrogram size, which should be used to make it more readable)
Heatmap(heat,
        cluster_rows = T,cluster_columns = T,     # set both parameters to false if manually setting sample order
        color='-RdYlBu2:100',row_labels = rownames(heat),right_annotation = row_ha,
        column_title = "PLOT TITLE",
        heatmap_legend_param = list(title = 'LEGEND TITLE'))

# CONCLUSION: I find that dendrograms and heatmap are especially informative for identifying outliers



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Gene Oultlier Insepction

# run differential expression analysis on DESeq object previously created (be mindful of the model design that was given above,
# in this case the differential expression analysis accounted for presence of Tcells and cell infection 'status')
dds<-DESeq(dds)

# extract cooks distance metric, which estimates the influence of a particular data point
mcols(dds)$maxCooks <- apply(assays(dds)[["cooks"]], 1, max)

 # create dataframe of base mean gene expression and maxCook metric
df<-data.frame(baseMean=mcols(dds)$baseMean,maxCooks=mcols(dds)$maxCooks) 

# extract row names to a new variable in the dataframe
df$gene<-row.names(df)

# obtain gene symbols for the gene IDs (in my case I was starting with ENSEMBLE gene IDs)
gsymb<-gconvert(df$gene,organism='mmusculus') # if needed, use this to convert gene IDs to gene symbols

# Plotting basemean versus maxCook
temp<-left_join(df,gsymb,by=c("gene"="input")) # combine gene symbols with base mean and maxCook dataframe
df$gene<-ifelse(!is.na(temp$name),temp$name,df$gene) # if gene symbol is NA, leave the ENSEMBL ID

# Create plot (using geom_text_repel to label only the outlier points) 
plot1<-ggplot(df,aes(x=baseMean,y=maxCooks))+geom_point()+geom_text_repel(
  data = subset(df, baseMean > 80000 | maxCooks>8),aes(label = gene),size = 3, # rerun plot a few times adjusting basemean and cooks threshold as desired
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"))
plot1

# Plotting fold change versus maxCooks
# NOTE: Be mindful of the model design that was given above, in this case the differential expression analysis accounted for presence of Tcells
# and cell infection 'status'. In the next step we will extract results from the DE analysis, and without specifying any further parameters for the 
# 'result' command we will be getting fold change for "last level of the last factor / first level of the last factor" of the model design we provided
res1<-results(dds) # extract results of DE analysis
stopifnot(all.equal(rownames(dds), rownames(res1))) # OPTIONAL: checks that the order of the row names is the same for both objects
df<-data.frame(FoldChange=res1$log2FoldChange,maxCooks=mcols(dds)$maxCooks) # create dataframe of fold change and maxCook
df$gene<-row.names(df) # extract gene ID to new column in data frame
temp<-left_join(df,gsymb,by=c("gene"="input")) # convert to gene symbols
df$gene<-ifelse(!is.na(temp$name),temp$name,df$gene) # if gene symbol is NA, leave the ENSEMBL ID
plot2<-ggplot(df,aes(x=FoldChange,y=maxCooks))+geom_point()+geom_text_repel(
  data = subset(df, FoldChange > 6 | maxCooks>5),aes(label = gene),size = 3,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"))
plot2

# Combine plots 1 and 2
grid.arrange(plot1,plot2, nrow=1,top="Gene Outliers")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Extreme gene boxplots with sample info indicated 

cooks<-(assays(dds)[["cooks"]]) # extract cooks metric
row.names(cooks)<-ifelse(!is.na(temp$name),temp$name,row.names(cooks)) # label rows with gene symbols (can use the 'temp' df created previously)
cooks<-t(cooks) # transpose the data frame
trim<-c("Gpr179","Ankle1","Srp54b","Nos2","Ccl17","Ankmy1","Zfp960") # create vector of outlier genes that were previously identified
library(car) # load 'car' package
Boxplot(cooks[,colnames(cooks) %in% trim],id.method="y", 
        ylab="cooks distance",xlab="gene symbol",main="gene outliers by sample")
        
CONCLUSION: do one or two samples account for all the outlier data points, or are they distributed across a large number of samples



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# For a particular sample, we can graph the genes that are most different from the base mean expression across all sampkles

# Looking at the top genes that make sample 'S11' different (see heatmap, it is an outlier)
S11<-data.frame(cnt=counts(dds)[,11],basemean=df$baseMean,gene=df$gene)
S11$diff<-S11$cnt-S11$basemean
S11$absDiff<-abs(S11$diff)
S11<-arrange(S11,-absDiff)

ggplot(S11,aes(x=cnt,y=basemean))+geom_point()+geom_text_repel(
  data = subset(S11, basemean > 10000 & cnt>40000),aes(label = gene),size = 3, # plot a few times and adjust thresholds as desired
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sample selection for DE analysis

#insert samples that will NOT be analyzed
fordeg<-c(
  "S1","S2",
  "S3","S4",
  "S37","S38","S39","S40","S41","S42","S43","S44","S45","S46","S47","S48",
  "S27","S28"
  #colData$id[colData$infection=='dH']
)

colData2<-colData[colData$id %!in% fordeg,] # removes selected samples from colData
countData2<-countData[,colnames(countData) %!in% fordeg] # removes selected samples from countData


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# construct DESeq object and perform DE analysis

# assign factor levels as desired
colData2$infection<-factor(colData2$infection,levels=c("Rv","dH")) # for Rv intercept
colData2$status<-factor(colData2$status,levels=c("iM","bM")) # for iM intercept
#colData2$group<-factor(colData2$group,levels=c("NT","T")) # for NT intercept

# set-up S4 object for differential expression analysis
dds2 <- DESeqDataSetFromMatrix(countData=countData2, 
                               colData=colData2, 
                               design=~infection+status+infection:status,  # 2 main effects with an interaction term
                               tidy=F)

dds2@design # confirm model formula (should be same as indicated above)
dds2 <- DESeq(dds2) # perform DE analysis
resultsNames(dds2) # print coefficient terms

### NOTE: Rv, iM, and NT are the condition intercepts


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# extract results from specific contrasts

# res1: difference between dH and Rv w/out Tcells for infected macs only
res1<-results(dds2,contrast = c('infection','dH','Rv'));summary(res1) # prints an analysis summary
res1df<-results(dds2,contrast = c('infection','dH','Rv'),tidy=T) # creates a data frame of results (infection main effect)

# subset results based on adjusted p-value and foldchange metrics
res1df <- res1df %>%
  arrange(padj) %>%
  subset(padj<0.1 & log2FoldChange>0)
  
# res2: difference between dH and Rv w/out Tcells for bystander macs only
resultsNames(dds2) # coefficiebnt terms -- Rv, iM, and NT are intercepts

res2<-results(dds2,contrast=list(c("infection_dH_vs_Rv","infectiondH.statusbM")));summary(res2)
res2df<-results(dds2,contrast = list(c("infection_dH_vs_Rv","infectiondH.statusbM")),tidy=T)

res2df <- res2df %>%
  arrange(padj) %>%
  subset(padj<0.1 & log2FoldChange>0)
  
# res3: is the main effect dH vs Rv different for bystander and infected cells (interaction term)

res3<-results(dds2,name = "infectiondH.statusbM");summary(res3)
res3df<-results(dds2,name = "infectiondH.statusbM"),tidy=T)


```



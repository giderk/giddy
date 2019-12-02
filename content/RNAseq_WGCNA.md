|Co-Expressed Hub Genes Sample Figure|
|---|
|![example](https://gideonerkenswick.files.wordpress.com/2019/12/cytoscapehubgenes_example.jpg)|
```rscript
# Script Name: RNAseq WGCNA

## Overview ##

#1) Weighted gene co-expression analysis (WGCNA) of RNAseq data. Groups genes into mutually exclusive 'modules' based on similarity in gene expression
#2) correlating modules with sample attributes

## required packages ##
library(WGCNA)
library(stringr) # character manipulation
library(AnnotationDbi) # gene annotation
library(Mus.musculus) # mouse reference
library(dplyr) # data manipulation
library(tidyr) # data manipulation
'%!in%' <- function(x,y)!('%in%'(x,y)) # 'NOT IN' function
library(DESeq2) # differential expression and read count normalization
library(nnet) #
library(logistf) #
library(ggplot2) #
library(GGally) #
library(clipr) # copy objects to clipboard


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

#vif controls samples will be included designate them as such for 'status' and 'infection' variables
colData[colData$status=="1"|colData$status=="2",]$status<-"CON"
colData[colData$infection=="Tc",]$infection<-"CON"



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# sample selection

# add samples IDs that will NOT be analyzed
fordeg<-c(
  "S1","S2",
  "S3","S4",
  "S37","S38","S39","S40","S41","S42","S43","S44","S45","S46","S47","S48",
  "S27","S28"
  ,colData$id[colData$infection=='dH']
)

colData2<-colData[colData$id %!in% fordeg,] # removes samples that are in 'fordeg' dataframe
countData2<-countData[,colnames(countData) %!in% fordeg] # removes samples that are in 'fordeg' dataframe



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read count normalization (using variance stabilizing transformation from DESeq2)

#set up DESeq matrix object
dds <- DESeqDataSetFromMatrix(countData=countData2, 
                              colData=colData2, 
                              design=~group+status+group:status,
                              tidy=F)

ddsVST<-varianceStabilizingTransformation(dds,blind=F)

normdata<-t(assay(ddsVST))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  setup WGCNA 

# check for genes that are data sufficient for WGCNA, genes with very low variance or too many missing values are flagged
gsg<-goodSamplesGenes(normdata,verbose = 3)

gsg$allOK # if this is TRUE, skip the next seubsetting section

# remove deficient genes and/or samples
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",paste(names(normdata)[!gsg$goodGenes],collapse=", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(normdata)[!gsg$goodSamples],collapse=", ")));
  readydata = normdata[gsg$goodSamples,gsg$goodGenes]
}

# look at simple sample tree cluster
sampleTree = hclust(dist(readydata),method="average") 

# plot the sample tree
plot(sampleTree,main="Sample clustering to detect outliers",sub="",xlab="",cex.lab=1.5,
     cex.axis=1.5,cex.main=2)

# this is where we could employ a tree clipping techniques to further remove outlier samples
# in this case, although S11 looks a little odd from Tcell samples, it still clusters to the 
# appropriate group

#create metadata
dim(colData2);names(colData2)
samples<-rownames(readydata) # extract sample names
traitRows<-match(samples,colData2$id) # matched metadata and expression data
datTraits<-colData2[traitRows,-1] # remove ID column from metadata
rownames(datTraits)<-colData2[traitRows,1] # add sample IDs as row names
collectGarbage()

#look at sample tree with metadata included
sampleTree2<-hclust(dist(readydata),method="average") 
traitColors = labels2colors(datTraits) # assign colors to metadata labels
plotDendroAndColors(sampleTree2,traitColors,groupLabels=names(datTraits),
                    main="Sample dendrogram and trait heatmap")

# we need to select appropriate powers for the WGCNA function (explore WGCNA package explanation for further details)
powers = c(c(1:10), seq(from = 12, to=20, by=2)) # provide a range of powers values to consider
sft = pickSoftThreshold(readydata, powerVector = powers, verbose = 5) # soft pick power function automates optimization

# plots power model fit
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red");

abline(h=0.90,col="red") # draw line at where R^2 level is .9, choose value at this line

# Mean connectivity plot (optional)
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=.9,col="red")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# run WGCNA

# this is necessary to make sure that WGCNA uses the correct 'cor' function (not the base stats 'cor' function)
cor <- WGCNA::cor # replace r stat cor function with WGCNA cor function

# WGCNA function (notice the power parameters, which was calculated in last code block)
net = blockwiseModules(readydata, power = 9,TOMType = "unsigned", 
                       minModuleSize = 30, maxBlockSize = 17000,
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,saveTOMs = TRUE,
                       saveTOMFileBase = "tbseq",verbose = 3)

# this is necessary to make sure that WGCNA uses the correct 'cor' function (not the base stats 'cor' function)
cor<-stats::cor # replace r stat 'cor' function


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualize co-expression network

table(net$colors) # look at the number of clusters formed
net$dendrograms[[1]] # dendrogram object stored here

# converts numbered module labels to colored module labels
mergedColors = labels2colors(net$colors); table(moduleColors) 

#plot the network
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# associate factors with gene sets in each module

nGenes = ncol(readydata) # total number of genes
nSamples = nrow(readydata) # total number of samples
MEs0 = moduleEigengenes(readydata, moduleColors)$eigengenes # eigen value/gene for each module
MEs0 = orderMEs(MEs0) # arrange eigen genes in the same order as samples
MEs<-net$MEs # 


# convert metadata to bivariate characters (not ideal when a factor has multiple non-ordinal levels)
# we do this so that we can comput a correlation between the modules and the sample variables
dat2<-data.frame(strain=ifelse(datTraits$infection=="Rv",0,1),
                 group=ifelse(datTraits$group=="T",1,0)
                 ,status=ifelse(datTraits$status=="iM",1,0))


#datTraits$infection<-as.factor(datTraits$infection) # IGNORE
#dat2 <- data.frame(lapply(datTraits, function(x) as.numeric(as.factor(x)))) # IGNORE

#creates correlations between each variable and each module
moduleTraitCor = cor(MEs0, dat2, use = 'p')

#adjust pvalues for multple comparisons
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#create lables of the cells of the heatmap
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor) # compute dimensions of heatmap

#plot the heatmap and inspect correlations and significant levels
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(dat2),
               yLabels = names(MEs0),ySymbols = names(MEs),
               colorLabels = FALSE,colors = greenWhiteRed(50),
               textMatrix = textMatrix,setStdMargins = FALSE,
               cex.text = 0.5,zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# OPTIONAL: use this for any categorical variables with more than two, nonordinal levels.

# if we have a variable with more than 2 levels, we can use multinomial logistic regresssion

List=list()
for (i in names(MEs)) {
  lg<-multinom(datTraits$status~MEs[[i]],data = MEs);
  z<-summary(lg)$coefficients/summary(lg)$standard.errors;
  p <- (1 - pnorm(abs(z), 0, 1)) * 2;
  List[[length(List)+1]] = p
}  
names(List)<-names(MEs)
List


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# for the binomial variables, a more appropriate measure of association with modules is to use a 
# logistic regresssion, but there may be problems with 'perfect separation' so we use 'logistf' package, which
# appropriately corrects for this. IT APPEARS THIS METHOD IS CONSISTENT WITH DOING A SIMPLE CORRELATION
# AS CONDUCTED ABOVE IN THE MAKING OF THE HEATMAP.

# run a logicstic regression for each of the modules against the binary factor that is indicated.

# strain variable
List=list()
for (i in names(MEs)) {
  lg<-logistf(dat2$strain~MEs[[i]],data = MEs,family="binomial");
  lg$prob
  print(lg)
  List[[length(List)+1]] = lg$prob
}
names(List)<-names(MEs)
Istrain<-t(as.data.frame(List))

# group variable
List=list()
for (i in names(MEs)) {
  lg<-logistf(dat2$group~MEs[[i]],data = MEs,family="binomial");
  lg$prob
  print(lg) 
  List[[length(List)+1]] = lg$prob
}  
names(List)<-names(MEs)
Tgrp<-t(as.data.frame(List))

#status variable
List=list()
for (i in names(MEs)) {
  lg<-logistf(dat2$status~MEs[[i]],data = MEs,family="binomial");
  lg$prob
  print(lg) 
  List[[length(List)+1]] = lg$prob
}  
names(List)<-names(MEs)
Istatus<-t(as.data.frame(List))

# merge Pvalues for each factor for plotting
modtraits<-data.frame(Istrain[,2],
                      Tgrp[,2],
                      Istatus[,2],
                      row.names = row.names(Tgrp) # !!!!need to correct this line!!!!
                      )
# rename the factor variables to be clear
names(modtraits)<-c("strain",
                    "Tgroup",
                    "cellStatus")

ggpairs(modtraits) # plots all pairwise combinations of the factors based on module Pvalues

sigMods<-subset(modtraits,#modtraits$strain <= 0.05 | 
                modtraits$Tgroup<=0.05 
                | modtraits$cellStatus <= 0.05)
sigMods # these are the gene modules with significant logistic regression Pvalues

# there are no modules that are both significant for Tgroup and Istatus, so we will extract
# the genesets for the significant modules and see if there are any genes that are shared

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# extract and write module gene sets that are significant for sample attributes

T1<-names(net$colors[net$colors==1])
T2<-names(net$colors[net$colors==2])
S0<-names(net$colors[net$colors==0]) # this grey (# zero) module, is reserved for
# unassigned genes,hence not a propoer module. Nonetheless, it is significantly associated with cell status factor.

write.csv(T1, "WGCNA/T1_modulegenes.csv",row.names = F) # T1 gene module set
write.csv(T2, "WGCNA/T2_modulegenes.csv",row.names = F) # T2 gene module set
write.csv(S0, "WGCNA/S0_modulegenes.csv",row.names = F) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# additional thoughts...

# we can play with the deepsplit option to allow for finer or courser module detection
#   - this didn't add anything
# we can combine the dH and Rv datasets, which increases sample size, and rerun the analysis
#   - combining dH and Rv results in more module splitting for significant Tcell groups only. we had no significant groups for other factors
# instead of combining dH and Rv, we can run a consensus module analysis between the two sets
#   - not likely to be informative



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hub genes

hubs<-chooseTopHubInEachModule(readydata,net$colors,type="unsigned")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# exporting network and expression data to cytoscape

TOM = TOMsimilarityFromExpr(readydata, power = 9);
modules=c(2)
probes<-colnames(readydata)
inModule = is.finite(match(net$colors, modules));
modProbes = probes[inModule];
modGenes = mapIds(Mus.musculus,
                  keys = modProbes,
                  keytype = "ENSEMBL",
                  column = "SYMBOL",
                  multiVals = 'first');
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,threshold = 0.02,
                               nodeNames = modProbes,altNodeNames = modGenes,
                               nodeAttr = net$colors[inModule])

chooseOneHubInEachModule(readydata,net$colors)

# so we extract the top Hub gene, and then we subset the edges file based on genes that interact with
# the hub. Then we subset to a manageable number of genes based on the weigh values, and we can consider 
# transforming the weight values to make more sense visually.

range(cyt$edgeData$weight)
edges<-read.table("Box/RNAseq/esxh/Ranalysis_noT.v.yesT/CytoscapeInput-edges-2.txt",header=T,sep="\t")
edges2<subset(edges,edges$weight>.8)
edges<-arrange(cyt$edgeData,weight)
edges<-edges
chooseOneHubInEachModule(readydata,net$colors)
temp.adj<-adjacency(readydata,selectCols = probes,power=9)

```





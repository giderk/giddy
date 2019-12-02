```rscript
# Script Name: DE-Gene Ontology Interactions

## NOTE ##
# this is a pipeline that can be peformed entirely with the ClusterProfiler package with a few commands. However,
# when more customization of the terms and dispaly is desired, this script will get the process started

## OVERVIEW ##

#1) load GO gene lists
#2) perform ORA of DE genes
#3) plot connections

## Required Packages ##

library(clusterProfiler)
library(AnnotationDbi)
library(Mus.musculus)
library(ggplot2)
library(igraph)
library(ggraph)
library(ggridges)

setwd("~")  # set working directory

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# import data

#import GO terms and associated genes, add term description variable, add goID variable 
go1<-read.table("GOsets/GO-0034341_resp.INFg.txt",header=F,sep="\t");go1$goID<-"GO-0034341";go1$name<-"resp.INFg"
go2<-read.table("GOsets/GO-0007259_STATcascade.txt",header=F,sep="\t");go2$goID<-"GO-0007259";go2$name<-"STATcascade"
go3<-read.table("GOsets/GO-0002237_resp.mol.bact.origin.txt",header=F,sep="\t");go3$goID<-"GO-0002237";go3$name<-"resp.mol.bact.origin"
go4<-read.table("GOsets/GO-0034612_resp.TNF.txt",header=F,sep="\t");go4$goID<-"GO-0034612";go4$name<-"resp.TNF"
go5<-read.table("GOsets/GO-0038127_ERBBsignaling.txt",header=F,sep="\t");go5$goID<-"GO-0038127";go5$name<-"ERBBsignaling"
go6<-read.table("GOsets/GO-0070555_response_IL1.txt",header=F,sep="\t");go6$goID<-"GO-0070555";go6$name<-"response_IL1"
go7<-read.table("GOsets/ITIzone_receptors.txt",header=F,sep="\t");go7$goID<-"GO-0000000";go7$name<-"receptors";names(go7)[1]<-"V3"


DE<-read.csv("path/to/DE/genelist.csv",header=T) # import DE gene list

# map symbols to gene IDs
DElist<-mapIds(Mus.musculus,
           keys = as.character(DE$geneid),
           keytype = "ENSEMBL",
           column = "SYMBOL",
           multiVals = 'first')

#import background gene list for ORA analysis
background<-read.table("Box/RNAseq/esxh/mouseref.txt",sep="\t",header=F)

# map symbols to gene IDs
bg<-mapIds(Mus.musculus,
           keys = as.character(background$V1),
           keytype = "ENSEMBL",
           column = "SYMBOL",
           multiVals = 'first')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# perform ORA data

GO<-rbind(go1,go2,go3,go4,go5,go6) # combine GO genesets

# create 'term2gene' object for ORA
term2gene<-subset(GO,select=c("goID","V3"))
term2gene<-rbind(term2gene,subset(go7,select=c("goID","V3")))

# create 'term2name' object for ORA
term2name<-subset(GO,select=c("goID","name"));term2name<-unique(term2name)
term2name<-rbind(term2name,c("GO-0000000","receptors"))

# perform ORA analysis
test<-enricher(gene=DElist,pvalueCutoff = 0.05,pAdjustMethod = "BH",
         universe = bg,TERM2GENE = term2gene,TERM2NAME = term2name)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plotting the gene to term relationships

geneSets<-geneInCategory(test) # extract gene and go term information, creating a list of vectors for each go term
y <- as.data.frame(test@result) # pick which list elements to plot
geneSets <- geneSets[y$ID] # select on the chosen list elements
names(geneSets) <- y$Description # rename the geneset list to the goTerm descriptions

g<-list2df(geneSets) # create a dataframe of all pairwise associations between genes and terms

graph <- graph.data.frame(g,direct=T) # convert df to object for ggraph

ggraph(graph,layout='fr') + 
  geom_edge_link(aes(colour = factor(g$categoryID))) + 
  geom_node_point(size=1)+geom_node_text(aes(label=name))
  
  ```

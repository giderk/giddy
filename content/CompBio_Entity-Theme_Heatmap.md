```rscript
# Script Name: CompBio_Entity-Theme_heatmap

## OVERVIEW ##

#1) Reconfigure compBio data export into functional data frame object
#2) Plot a theme by gene heatmap for selected genes

## NOTE ##
# CompBio is a gene/protein set analysis tool developed by GTAC of WashU

## Required Packages ##

library(stringr) # character string manipulation
library(reshape2) # data rearrangement
library(plyr) # data manipulation
library(dplyr) # data manipulation
library(ComplexHeatmap) # heat map functions
library(circlize) # heatmap customization

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# load compBio exported data
setwd("~")
df<-read.csv("Desktop/ITIzone_upregulated.csv",header=T,stringsAsFactors = F)
colnames(df)<-str_remove(colnames(df),'Name.') # clean up theme name
df2<-lapply(df,as.data.frame) # convert each column to list element
df3<-ldply(df2) # convert list elements back to long dataframe with a theme variable
names(df3)<-c("Theme","element") # rename variables

df3$type<-str_extract_all(df3$element,'^.*=') # create data type column
df3$type<-str_remove_all(df3$type,'=') # cleanup 'type' column
df3$element<-str_remove_all(df3$element,'^.*=');head(df3) # cealnup 'element' column
ldf3<-split(df3,df3$type) # split dataframe on 'type' variable
df4<-ldf3[[5]] # extract entity info into new dataframe
df4$value<-str_extract(df4$element,'\\((.*?)\\)') # extract entity score to 'value'
df4$value<-str_remove_all(df4$value,'(\\()|(\\))');df4$value<-as.numeric(df4$value) # remove parenthesis
df4$gene<-str_remove_all(df4$element,'\\((.*?)\\)') # extract gene symbol to new column

df5<-dcast(df4,formula=Theme~gene,value.var = 'value') # reshape data into matrix theme x gene
df5[is.na(df5)]<-0 # replace NAs with 0s
mat<-as.matrix(df5[-1]);rownames(mat)<-df5$Theme # convert to matrix, themes as row names

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# import genes of interest
receptors<-read.csv("path/to/genes/of/interest.txt",header=F,stringsAsFactors = F)
receptors<-receptors$V1; receptors<-tolower(receptors) # convert symbols to lower-case
mat<-mat[,colnames(mat) %in% receptors] # subset matrix to genes of interest
mat<-mat[match(ldf3$Score$Theme,rownames(mat)),] # order matrix by descending theme score
rownames(mat)<-paste(ldf3$Score$Theme,ldf3$Score$element,sep = "_") # add theme score to theme names

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot heatmap
Heatmap(mat[,c(1:140)],row_order = NULL,column_dend_height = unit(30, "mm"),na_col = "grey",
        heatmap_legend_param = list(title = "Score"),column_title = "ITIzone_compBIO.themes_up-regulated.receptors",
        col=colorRamp2(c(0, 1, 100,1000,2000), c("black", "purple", "white","orange",'red')))

```

# Meta-analysis: how group composition impacts tamarin breeding success

|Sample Forest Plot|
|---|
|![example](https://gideonerkenswick.files.wordpress.com/2019/12/meta.corr_sample.jpg)|

```rscript

setwd("~")

### Required Packages ###
library(xlsx) # no longer supported, import with read.table or read.csv
library(plyr)
library(lmtest)
library(car)
library(reshape2)
library(metafor)

data<-read.xlsx("Demographics6.xlsx",1) # import data

data<-subset(data,data$Complete=="y");data<-droplevels(data) # drop incomplete records

# calculate proportions for each age class of each group
data<-ddply(data,.(Species),mutate,pM=AM/grp.size,pF=AF/grp.size,pSM=SM/grp.size,
            pSF=SF/grp.size, pTJ=TJ/grp.size)

data[is.na(data)]<-0 # replace NAs with zeros

# create short codes (id) to separate studies
data$study<-as.factor(transform(data$Reference, study= as.numeric(factor(data$Reference)))[,2])

# quantify sample size for each study
samplesize<-ddply(data,.(Reference),nrow);samplesize[samplesize$V1>4,]

# calculate correlation between proportion of sex-age groups and dependents, by study

effect1<-ddply(data,.(Reference),summarize,AMc=cor(AM,Dependants, method="spearman"),
      AFc=cor(AF,Dependants, method="spearman"),
      SMc=cor(SM,Dependants, method="spearman"),
      SFc=cor(SF,Dependants, method="spearman"),
      GRPc=cor(grp.size,Dependants,method="spearman"),
      AF.GRP=cor(grp.size,AF,method="spearman"))

effect2<-merge(effect1,samplesize, by=c("Reference")) # merge output with samples sizes
```

### Only use this code if moderator variables are desired

```rscript
# standardize the species names
levels(effect2$Species)

effect2$Species[effect2$Species=="SFUS" |
            effect2$Species=="S. fuscicollis illigeri" |
              effect2$Species=="S. weddelli" |
              effect2$Species=="S. fuscicolis"]<-"S. fuscicollis"

effect2$Species[effect2$Species=="SIMP"]<-"S. imperator"

effect2<-droplevels(effect2);levels(effect2$Species)

# make the study moderator a number variable
effect2$mod<-transform(effect2$Species, study= as.numeric(factor(effect2$Species)))[,2]
```
### Perform meta-analysis using correlations as the effect size:
* use ZCOR to transform raw correlations to the z-scale (standardization)
* data points weighted by sample size
* random-effect model: HS (hunter-schmidt) method to estimate heterogeneity

```rscript
# removing all data points with N < 5 (model assumption)
data2<-subset(effect2,effect2$V1>4);data2

# adult females
dat <- escalc(measure="ZCOR", ri=AFc, ni=V1,
              data=data2, vtype="HO");

res_AF <- rma(yi, vi, weights=sqrt(V1), data=dat, 
              method="HS", slab=Reference); res_AF

forest(res_AF, transf=transf.ztor.int, order = "obs",
       showweights = T, cex=1.2,
       ilab=cbind(data2$V1),
       ilab.xpos=(-2.5))
op <- par(cex=1.2, font=2)
text(c(-3.9,-2.4,3.6,4.9), 17.5, c("Study", "n", "Weight", "effect[LCL, UCL]"),pos=2)
par(op)

# save as PDF 8x9.5in

# adult males

dat <- escalc(measure="ZCOR", ri=AMc, ni=V1,
              data=data2, vtype="HO");

res_AM <- rma(yi, vi, weights=sqrt(V1), data=dat, method="HS",
              slab=Reference);res_AM

forest(res_AM,transf=transf.ztor.int,order = "obs",
       showweights = TRUE,cex=.75,
       ilab=cbind(data2$V1),
       ilab.xpos=(-2.5))
op <- par(cex=.75, font=2)
text(c(-3.9,-2.4,3.5,4.8), 18.5, c("Study", "n", "Weight", "effect [LCL, UCL]"),pos=2)
par(op)

# sub-adult females

dat <- escalc(measure="ZCOR", ri=SFc, ni=V1,
              data=data2, vtype="HO");
res_SF <- rma(yi, vi, weights=sqrt(V1), data=dat, method="HS", slab=Reference);res_SF
forest(res_SF)

# sub-adult males

dat <- escalc(measure="ZCOR", ri=SMc, ni=V1,
              data=data2, vtype="HO");
res_SM <- rma(yi, vi, weights=sqrt(V1), data=dat, method="HS", slab=Reference);res_SM
forest(res_SM)


#group size
dat <- escalc(measure="ZCOR", ri=GRPc, ni=V1,
              data=data2, vtype="HO");
res_GRP <- rma(yi, vi, weights=sqrt(V1), data=dat, method="HS",
               slab=Reference);res_GRP
forest(res_GRP,transf=transf.ztor.int, order = "obs",
       showweights = TRUE,cex=1.2,
       ilab=cbind(data2$V1),
       ilab.xpos=(-2.5))
op <- par(cex=1.2, font=2)
text(c(-3.9,-2.4,3.6,4.8), 18.5, c("Study", "n", "Weight", "effect [LCL, UCL]"),pos=2)
par(op)


#group size correlated with adult females

dat <- escalc(measure="ZCOR", ri=AF.GRP, ni=V1,
              data=data2, vtype="HO");
res_AF.GRP <- rma(yi, vi, weights=V1, data=dat, method="HS", slab=Reference);res_AF.GRP
forest(res_AF.GRP, transf=transf.ztor.int)
funnel(res_AF.GRP)
```

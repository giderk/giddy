# Modeling primate breeding status from reproductive morphology
|Before and After BR Model|
|---|
|![example](https://gideonerkenswick.files.wordpress.com/2019/12/brmodel_outcomes_sample.jpg)|
Used with data from wild S. imperator and L. weddelli

### Overview:
* Data preparation
* Dimension reduction (PCA)
* Breeding status assignment (lDFA)
* Model performance

```rscript

setwd("~/") # set working directory

## Require Packages ##
require(xlsx)
require(plyr)
require(reshape2)
require(ggplot2)
library(GGally)
library(MASS)
require(gridExtra)
require(FactoMineR)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# import data
morph<-read.xlsx("AnimalExams3.9.16.xlsx",1) # load morphology data

# remove unwanted variables
names(morph)
morph<-subset(morph, select=c(1:8,11,14:15,21:24,26:29))

# convert relevent variables to factors
morph$AnID<-as.factor(morph$AnID) #makes AnID a factor
morph$ExamID<-as.factor(morph$ExamID) # makes ExamID a factor

# Pull the year into a separate variable
morph$year<-as.factor(format(morph$Date, "%Y"))

#remove all records from 2009
morph<-subset(morph, morph$year!='2009')

#check for 0s and NAs in weight variable, we will be removing
# an Sf-GBR row that has no data, and keeping an Sf-OPG row
# that is only missing weight (we will give her the group average)
morph[morph$Weight<1,]
morph[is.na(morph$Weight),]

# remove individuals with weight = 0
morph<-subset(morph, morph$Weight!=0 | morph$ExamID==14)
morph<-droplevels(morph)

# insert 0s for NAs in nipple length, ifelse statements don't work with NAs
morph$NippleLength_L[is.na(morph$NippleLength_L)]<-0
morph$NippleLength_R[is.na(morph$NippleLength_R)]<-0

# Split the data by the group and sex, since lumping these data together would
# give odd results (emps are bigger than fuscis, and females and males don't have same
# genitalia)
l<-split(morph,list(morph$Sex, morph$Species))

# While the grouped data is in a list, we can add variables and do calculations as needed

# determine animal weight by substracting bag weight
l <- lapply(l, transform, nWeight = Weight-Bag.weight) 
# create elipsoid measurement of tesicular volumn following Watsa 2012
l <- lapply(l, transform, vTestes = round(((pi*TestesLength*(TestesWidth/2)^2)/6),2))
# create vulva index
l <- lapply(l, transform, iVulva = round((VulvaLength+VulvaWidth),2))
# calculatate SP gland area
l <- lapply(l, transform, aSP = round(SP_Length*SP_Width,2))
# combine nipple measurements
l<- lapply(l,transform,aNip = ifelse(NippleLength_L>0 & NippleLength_R>0,
                                     ((NippleLength_L+NippleLength_R)/2),
                                     ifelse(NippleLength_L>0 & NippleLength_R==0,
                                            NippleLength_L,
                                            ifelse(NippleLength_L==0 & NippleLength_R>0,
                                                   NippleLength_R,0))))

# using unsplit function return all the data to one data frame
morph2<-unsplit(l,list(morph$Sex, morph$Species))

# remove the subpopulations of data, and drop levels that are no longer present
fSfus<-as.data.frame(l[[1]]); fSfus<-droplevels(fSfus)
mSfus<-as.data.frame(l[[2]]); mSfus<-droplevels(mSfus)
fSimp<-as.data.frame(l[[3]]); fSimp<-droplevels(fSimp)
mSimp<-as.data.frame(l[[4]]); mSimp<-droplevels(mSimp)


# remove the variables that don't apply to the male and female groups
names(fSfus)
fSfus<-subset(fSfus, select=-c(14:15,22))
fSimp<-subset(fSimp, select=-c(14:15,22))

names(mSfus)
mSfus<-subset(mSfus, select=-c(12:13,16:17,23,25))
mSimp<-subset(mSimp, select=-c(12:13,16:17,23,25))

# replace erroneous nipple measurments (aNip) for  records...
fSfus[fSfus$ExamID==344,]$aNip<-3.38 #exam ID 344 (sf-yrc from 2015)
fSimp[fSimp$ExamID==63,]$aNip<-5.3 #exam ID 63 (si-GRC from 2015)
fSimp[fSimp$ExamID==333,]$aNip<-5.3 #exam ID 333 (src/sps from 2015)
fSimp[fSimp$ExamID==94,]$aNip<-5.3 # exam ID 94 (si OPL/ORC from 2015)
fSimp[fSimp$ExamID==347,]$aNip<-5.0 # exam ID 347 (si LPL-LRC from 2015)

# replace erroneous aSP measurements...
fSimp[fSimp$ExamID==347,]$aSP<-NA # exam ID 347 (si LPL-LRC from 2015)

# replace erroneous vTestes measurments for records...
mSfus[mSfus$ExamID==268,]$vTestes<-NA #examID 268 (sf OBR 2010)
mSfus[mSfus$ExamID==265,]$vTestes<-NA #examID 265 (sf P6T1 2010)
mSimp[mSimp$ExamID==336,]$vTestes<-NA #examID 336 (si YBS2 2015)


# check measurments are changed (11 corrections in total, of 331 capture
# instances)
fSfus[fSfus$ExamID==344,]
fSimp[fSimp$ExamID==63,]
fSimp[fSimp$ExamID==333,]
fSimp[fSimp$ExamID==319,]
fSimp[fSimp$ExamID==94,]
fSimp[fSimp$ExamID==92,]
fSimp[fSimp$ExamID==89,]
fSimp[fSimp$ExamID==223,]
mSfus[mSfus$ExamID==268,]
mSfus[mSfus$ExamID==265,]
mSimp[mSimp$ExamID==336,]


# for each age-sex cohort, and by species, replace missing morphs with mean value for their cohorts, respectively.

fSfus[!complete.cases(fSfus),] # shows not many incomplete records

# creates a list of female fuscis, separated by age group
ff<-split(fSfus,list(fSfus$Age_class))

# fills NAs with mean for each variable for that age-group
ff <- lapply(ff, transform, iVulva = 
               ifelse(is.na(iVulva), 
                      mean(iVulva, na.rm=TRUE), iVulva))
ff <- lapply(ff, transform, aSP = 
               ifelse(is.na(aSP), 
                      mean(aSP, na.rm=TRUE), aSP))
ff <- lapply(ff, transform, aNip = 
               ifelse(is.na(aNip), 
                      mean(aNip, na.rm=TRUE), aNip))
ff <- lapply(ff, transform, nWeight = 
               ifelse(is.na(nWeight), 
                      mean(nWeight, na.rm=TRUE), nWeight))

# re-combines the females fusci data, NAs have been fixed
fSfus<-unsplit(ff,list(fSfus$Age_class))

fSimp[!complete.cases(fSimp),] # shows not many incomplete records

# repeat for female imps
fSimp[!complete.cases(fSimp),]
fs<-split(fSimp,list(fSimp$Age_class))
fs <- lapply(fs, transform, iVulva = 
               ifelse(is.na(iVulva), 
                      mean(iVulva, na.rm=TRUE), iVulva))
fs <- lapply(fs, transform, aSP = 
               ifelse(is.na(aSP), 
                      mean(aSP, na.rm=TRUE), aSP))
fs <- lapply(fs, transform, aNip = 
               ifelse(is.na(aNip), 
                      mean(aNip, na.rm=TRUE), aNip))
fs <- lapply(fs, transform, nWeight = 
               ifelse(is.na(nWeight), 
                      mean(nWeight, na.rm=TRUE), nWeight))
fSimp<-unsplit(fs,list(fSimp$Age_class))

# repeat for male fuscis
mSfus[!complete.cases(mSfus),] 
mf<-split(mSfus,list(mSfus$Age_class))
mf <- lapply(mf, transform, vTestes = 
               ifelse(is.na(vTestes), 
                      mean(vTestes, na.rm=TRUE), vTestes))
mf <- lapply(mf, transform, aSP = 
               ifelse(is.na(aSP), 
                      mean(aSP, na.rm=TRUE), aSP))
mSfus<-unsplit(mf,list(mSfus$Age_class))

# repeat for male imps
mSimp[!complete.cases(mSimp),]
mi<-split(mSimp,list(mSimp$Age_class))
mi <- lapply(mi, transform, vTestes = 
               ifelse(is.na(vTestes), 
                      mean(vTestes, na.rm=TRUE), vTestes))
mi <- lapply(mi, transform, aSP = # code changed to assign 0 for NAs
               ifelse(is.na(aSP), 
                      0, aSP))
mSimp<-unsplit(mi,list(mSimp$Age_class))
```

> we now have out working datasets for each sex of each species, and can perform dimension reduction using PCA. Coordinates values will be used in assigning breed status.

```rscript
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# create new versions of each working dataset
fSfus2<-fSfus
fSimp2<-fSimp
mSfus2<-mSfus[mSfus$Br_status!="None",] # remove male juveniles
mSimp2<-mSimp[mSimp$Br_status!="None",] # remove male juveniles

# subset data to the morphological variables of interest
names(fSfus2)
ffpca<-subset(fSfus2,select=c(19:22))
fipca<-subset(fSimp2,select=c(19:22))

names(mSfus2)
mfpca<-subset(mSfus2,select=c(17:19))
mipca<-subset(mSimp2, select=c(17:19))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Principal component analysis

FFpca = PCA(na.omit(ffpca), scale.unit=TRUE, ncp=5, graph=F) # Female weddellis PCA
FIpca = PCA(na.omit(fipca), scale.unit=TRUE, ncp=5, graph=F) # Female imps PCA
MFpca = PCA(na.omit(mfpca), scale.unit=TRUE, ncp=5, graph=T) # Male weddellis PCA
MIpca = PCA(na.omit(mipca), scale.unit=TRUE, ncp=5, graph=F) # Male imps PCA

# extract loadings
a1<-sweep(FFpca$var$coord,2,sqrt(FFpca$eig[1:ncol(FFpca$var$coord),1]),FUN="/")
a2<-sweep(FIpca$var$coord,2,sqrt(FIpca$eig[1:ncol(FIpca$var$coord),1]),FUN="/")
a3<-sweep(MFpca$var$coord,2,sqrt(MFpca$eig[1:ncol(MFpca$var$coord),1]),FUN="/")
a4<-sweep(MIpca$var$coord,2,sqrt(MIpca$eig[1:ncol(MIpca$var$coord),1]),FUN="/")

# combine PCA loadings into single DF, add variable to distinguish each animal cohort
loading<-rbind.fill(data.frame(a1),data.frame(a2),data.frame(a3),data.frame(a4))
Loading<-cbind(loading,PCA=c(rep(c("FF","FI"), each=4),rep(c("MF","MI"),each=3)))

# write data frame to file
write.xlsx(Loading,"Loadings.xlsx")

#extract eigen values
eigen<-cbind(
  rbind.fill(FFpca$eig,FIpca$eig,MFpca$eig,MIpca$eig),PCA=c(
  rep(c("Sw female","Si female"),each=4),
  rep(c("Sw male","Si male"), each=3)))

# write eigen values to file
write.xlsx(eigen,"eigen.xlsx")


# extracting correlations of each variable with dimensions
var.cor<-cbind(rbind.fill(data.frame(FFpca$var$cor),
           data.frame(FIpca$var$cor),
           data.frame(MFpca$var$cor),
           data.frame(MIpca$var$cor)),PCA=c(
             rep(c("Sw female","Si female"),each=4),
             rep(c("Sw male","Si male"), each=3)))

# write correlations to file
write.xlsx(var.cor,"var.cor.xlsx")


# extract component coordinate values for cohort and add them to the working datasets (see above)
a<-as.data.frame(FFpca$ind$coord); ff<-cbind(fSfus2,a)
b<-as.data.frame(FIpca$ind$coord); fi<-cbind(fSimp2,b)
c<-as.data.frame(MFpca$ind$coord); mf<-cbind(mSfus2,c)
d<-as.data.frame(MIpca$ind$coord); mi<-cbind(mSimp2,d)

# Average scores for individuals that have repeat measures
fSfus3<-ddply(subset(ff,ff$Br_status!="Unknown"), .(AnID,Sex,Br_status),summarize,
              Dim1=mean(Dim.1),Dim2=mean(Dim.2));fSfus3<-droplevels(fSfus3)

fSimp3<-ddply(subset(fi,fi$Br_status!="Unknown"), .(AnID,Sex,Br_status),summarize,
              Dim1=mean(Dim.1),Dim2=mean(Dim.2));fSimp3<-droplevels(fSimp3)

mSfus3<-ddply(subset(mf,mf$Br_status!="Unknown"), .(AnID,Sex,Br_status),summarize,
              Dim1=mean(Dim.1),Dim2=mean(Dim.2));mSfus3<-droplevels(mSfus3)

mSimp3<-ddply(subset(mi,mi$Br_status!="Unknown"), .(AnID,Sex,Br_status),summarize,
              Dim1=mean(Dim.1),Dim2=mean(Dim.2));mSimp3<-droplevels(mSimp3)

# test that variances between groups are not statistically different
# which is an assumption of linear discriminant function analysis

bartlett.test(Dim1~Br_status,data=fSfus3); ggplot(fSfus3,aes(y=Dim1,x=Br_status))+geom_boxplot()
bartlett.test(Dim2~Br_status,data=fSfus3); ggplot(fSfus3,aes(y=Dim2,x=Br_status))+geom_boxplot()

bartlett.test(Dim1~Br_status,data=fSimp3);ggplot(fSimp3,aes(y=Dim1,x=Br_status))+geom_boxplot()
bartlett.test(Dim2~Br_status,data=fSimp3);ggplot(fSimp3,aes(y=Dim2,x=Br_status))+geom_boxplot()

bartlett.test(Dim1~Br_status,data=mSfus3); ggplot(mSfus3,aes(y=Dim1,x=Br_status))+geom_boxplot()
bartlett.test(Dim2~Br_status,data=mSfus3); ggplot(mSfus3,aes(y=Dim2,x=Br_status))+geom_boxplot()

bartlett.test(Dim1~Br_status,data=mSimp3); ggplot(mSimp3,aes(y=Dim1,x=Br_status))+geom_boxplot()
bartlett.test(Dim2~Br_status,data=mSimp3); ggplot(mSimp3,aes(y=Dim2,x=Br_status))+geom_boxplot()


# males violate homogeneity of variance when juveniles are included, they can be removed,
# as there is no uncertaintly regarind their breeding status.

# females do not violate homogeneity (but this may change as more samples are incorporated)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# discriminant function analysis

######## Female SFUS
# train model with subset of individuals of known breeding status
fSfus3.df<-lda(Br_status ~ Dim1 + Dim2,
              data = fSfus3);fSfus3.df

######## for jacknife predictions, use the following....
# fSfus3.df<-lda(Br_status ~ nWeightM + aSPM + iVulvaM+aNipM,
#               data = fSfus3,CV=T);fSfus3.df
# percent correct for each category of G
# ct <- table(fSfus3$Br_status, fSfus3.df$class)
# diag(prop.table(ct, 1))
# total percent correct
# sum(diag(prop.table(ct)))
###### jacknife section end ####

fSfus4<-subset(ff,select=c(23:24)) # subset data frame to required variables
names(fSfus4)<-names(fSfus3[,c(4,5)]) # correct the names of the variables

# predict breeeding status for entire cohort, including unknown individuals
fSfus4.df<-predict(fSfus3.df,newdata =fSfus4);fSfus4.df
fSfus5<-as.data.frame(fSfus4.df) # convert to a dataframe
cbind(ff=row.names(ff),fSfus4=row.names(fSfus5)) # make sure row names are identical

fSfus6<-cbind(ff,br.dfa=fSfus5$class)

# check that breeding status groups are distinct with MANOVA
names(fSfus6)
summary(manova(as.matrix(fSfus6[,19:22])~fSfus6$br.dfa),"Wilks")


############# Females Simp 
# Repeating what was done above for female simps

fSimp3.df<-lda(Br_status ~ Dim1 + Dim2,data = fSimp3);fSimp3.df # run DFA 

######## for jacknife predictions, use the following....
# fSimp3.df<-lda(Br_status ~ nWeightM + aSPM + iVulvaM + aNipM,
#               data = fSimp3,CV=T);fSimp3.df
# percent correct for each category of G
# ct <- table(fSimp3$Br_status, fSimp3.df$class)
# diag(prop.table(ct, 1))
# total percent correct
# sum(diag(prop.table(ct)))
###### jacknife section end ####


fSimp4<-subset(fi,select=c(23:24))
names(fSimp4)<-names(fSimp3[,c(4,5)])

fSimp4.df<-predict(fSimp3.df,newdata =fSimp4);fSimp4.df
fSimp5<-as.data.frame(fSimp4.df)
cbind(fi=row.names(fi),fSimp4=row.names(fSimp5))
fSimp6<-cbind(fi,br.dfa=fSimp5$class)

names(fSimp6)
summary(manova(as.matrix(fSimp6[,19:22])~fSimp6$br.dfa),"Hotelling-Lawley")


########### Male Sfus #############
# Repeat for Male fuscis

names(mSfus3)
mSfus3.df<-lda(Br_status ~ Dim1+Dim2, data = mSfus3);mSfus3.df # run DFA 

######## for jacknife predictions, use the following....
# mSfus3.df<-lda(Br_status ~ nWeightM + aSPM + vTestesM,
#               data = mSfus3,CV=T);mSfus3.df
# percent correct for each category of G
# ct <- table(mSfus3$Br_status, mSfus3.df$class)
# diag(prop.table(ct, 1))
# total percent correct
# sum(diag(prop.table(ct)))
###### jacknife section end ####

mSfus4<-subset(mf,select=c(20:21))
names(mSfus4)<-names(mSfus3[,c(4,5)])

# scale all the variables
#mSfus4$nWeightM<-scale(mSfus4$nWeightM)
#mSfus4$aSPM<-scale(mSfus4$aSPM)
#mSfus4$vTestesM<-scale(mSfus4$vTestesM)

mSfus4.df<-predict(mSfus3.df,newdata =mSfus4);mSfus4.df
mSfus5<-as.data.frame(mSfus4.df)
cbind(mf=row.names(mf),mSfus4=row.names(mSfus5))
mSfus6<-cbind(mf,br.dfa=mSfus5$class)

names(mSfus6)
summary(manova(as.matrix(mSfus6[,17:19])~mSfus6$br.dfa),"Hotelling-Lawley")


########### Male Simp #############
# Repeating for male SIMP

mSimp3.df<-lda(Br_status ~ Dim1+Dim2, data = mSimp3);mSimp3.df

######## for jacknife predictions, use the following....
# mSimp3.df<-lda(Br_status ~ nWeightM + aSPM + vTestesM,
#               data = mSimp3,CV=T);mSimp3.df
# percent correct for each category of G
# ct <- table(mSimp3$Br_status, mSimp3.df$class)
#vdiag(prop.table(ct, 1))
# total percent correct
#vsum(diag(prop.table(ct)))
###### jacknife section end ####

mSimp4<-subset(mi,select=c(20:21))
names(mSimp4)<-names(mSimp3[,c(4,5)])

# scale all the variables
#mSimp4$nWeightM<-scale(mSimp4$nWeightM)
#mSimp4$aSPM<-scale(mSimp4$aSPM)
#mSimp4$vTestesM<-scale(mSimp4$vTestesM)

mSimp4.df<-predict(mSimp3.df,newdata =mSimp4);mSimp4.df
mSimp5<-as.data.frame(mSimp4.df)
cbind(mi=row.names(mSimp2),mi=row.names(mSimp5))
mSimp6<-cbind(mi,br.dfa=mSimp5$class)

names(mSimp6)
summary(manova(as.matrix(mSimp6[,17:19])~mSimp6$br.dfa),"Wilks")

# combine all breeding assignments together
br.dfa<-rbind.fill(fSfus6,fSimp6,mSfus6,mSimp6)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ploting the 'before' and 'after'

### graph 'after' DFAs
fSfus5.1<-cbind(fSfus2,fSfus5[,c(1,5,6)])
g1<-ggplot(fSfus5.1, aes(x=x.LD1,y=x.LD2, color=class))+geom_point()+
  theme(panel.grid.major = element_blank(),
        axis.line=element_line(color="black"))+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  scale_color_manual(values = c("red", "limegreen", "blue","black"))+
  guides(color=F)+
  labs(x="LDF 1",y="LDF 2",
       title=expression(paste("Female ",italic("S. weddelli"), ": After")))

fSimp5.1<-cbind(fSimp2,fSimp5[,c(1,5,6)])
g2<-ggplot(fSimp5.1, aes(x=x.LD1,y=x.LD2, color=class))+geom_point()+
  theme(panel.grid.major = element_blank(),
        axis.line=element_line(color="black"))+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  scale_color_manual(values = c("red", "limegreen", "blue","black"))+
  guides(color=F)+
  labs(x="LDF 1",y="LDF 2",
       title=expression(paste("Female ",italic("S. imperator"), ": After")))

mSfus5.1<-cbind(mSfus2,mSfus5[,c(1,4)])
g3<-ggplot(mSfus5.1, aes(x=LD1,y=class, color=class))+geom_jitter()+
  theme(panel.grid.major = element_blank(),
        axis.line=element_line(color="black"))+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  scale_color_manual(values = c("limegreen", "blue","black"))+
  guides(color=F)+
  labs(x="LDF 1",y="Status",
       title=expression(paste("Male ",italic("S. weddelli"), ": After")))

mSimp5.1<-cbind(mSimp2,mSimp5[,c(1,4)])
g4<-ggplot(mSimp5.1, aes(x=LD1,y=class, color=class))+geom_jitter()+
  theme(panel.grid.major = element_blank(),
        axis.line=element_line(color="black"))+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  scale_color_manual(values = c("limegreen", "blue","black"))+
  guides(color=F)+
  labs(x="LDF 1",y="Status",
       title=expression(paste("Male ",italic("S. imperator"), ": After")))

grid.arrange(g1,g2,g3,g4,nrow=2)

### graph 'BEFORE' DFAs
fSfus5.1<-cbind(fSfus2,fSfus5[,c(1,5,6)])
b1<-ggplot(fSfus5.1, aes(x=x.LD1,y=x.LD2, color=Br_status))+geom_point()+
  theme(panel.grid.major = element_blank(),
        axis.line=element_line(color="black"))+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  scale_color_manual(values = c("red", "limegreen", "blue","black"))+
  guides(color=F)+
  labs(x="LDF 1",y="LDF 2",
       title=expression(paste("Female ",italic("S. weddelli"), ": Before")))

fSimp5.1<-cbind(fSimp2,fSimp5[,c(1,5,6)])
b2<-ggplot(fSimp5.1, aes(x=x.LD1,y=x.LD2, color=Br_status))+geom_point()+
  theme(panel.grid.major = element_blank(),
        axis.line=element_line(color="black"))+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  scale_color_manual(values = c("red", "limegreen", "blue","black"))+
  guides(color=F)+
  labs(x="LDF 1",y="LDF 2",
       title=expression(paste("Female ",italic("S. imperator"), ": Before")))

mSfus5.1<-cbind(mSfus2,mSfus5[,c(1,4)])
b3<-ggplot(mSfus5.1, aes(x=LD1,y=class, color=Br_status))+geom_jitter()+
  theme(panel.grid.major = element_blank(),
        axis.line=element_line(color="black"))+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  scale_color_manual(values = c("limegreen", "blue","black"))+
  guides(color=F)+
  labs(x="LDF 1",y="Status",
       title=expression(paste("Male ",italic("S. weddelli"), ": Before")))

mSimp5.1<-cbind(mSimp2,mSimp5[,c(1,4)])
b4<-ggplot(mSimp5.1, aes(x=LD1,y=class, color=Br_status))+geom_jitter()+
  theme(panel.grid.major = element_blank(),
        axis.line=element_line(color="black"))+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  scale_color_manual(values = c("limegreen", "blue","black"))+
  guides(color=F)+
  labs(x="LDF 1",y="Status",
       title=expression(paste("Male ",italic("S. imperator"), ": Before")))

grid.arrange(b1,g1,b2,g2,b3,g3,b4,g4,ncol=2)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# model performance

# Sfus female
fSfus7<-subset(fSfus6,fSfus6$Br_status!="Unknown")
ct <- table(fSfus7$Br_status, fSfus7$br.dfa)
# percent correct per category
SwF<-diag(prop.table(ct, 1))
# total percent correct
t1<-sum(diag(prop.table(ct)))

#Simp Female
# now let's look at how many were predicted correctly by category
fSimp7<-subset(fSimp6,fSimp6$Br_status!="Unknown")
ct <- table(fSimp7$Br_status, fSimp7$br.dfa)
SiF<-diag(prop.table(ct, 1))
# total percent correct
t2<-sum(diag(prop.table(ct)))

#Sfus Male
# now let's look at how many were predicted correctly by category
mSfus7<-subset(mSfus6,mSfus6$Br_status!="Unknown");mSfus7<-droplevels(mSfus7)
ct <- table(mSfus7$Br_status, mSfus7$br.dfa)
SwM<-diag(prop.table(ct, 1));SwM
# total percent correct
t3<-sum(diag(prop.table(ct)));t3

#Simp Male
# now let's look at how many were predicted correctly by category
mSimp7<-subset(mSimp6,mSimp6$Br_status!="Unknown");mSimp7<-droplevels(mSimp7)
ct <- table(mSimp7$Br_status, mSimp7$br.dfa)
SiM<-diag(prop.table(ct, 1));SiM
# total percent correct
t4<-sum(diag(prop.table(ct)));t4

data.frame(rbind(SwF,SiF),Cum.avg=c(t1,t2))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# write new animal file with predicted breeding status for each individual

# add non-breeding males back into the working dataset
head(br.dfa)
BreedingStatus<-rbind.fill(br.dfa,mSfus[mSfus$Br_status=="None",],mSimp[mSimp$Br_status=="None",])
BreedingStatus$br.dfa[is.na(BreedingStatus$br.dfa)]<-"None"

#write to an excel file
write.xlsx(BreedingStatus,"BreedingStatus.xlsx")
```
> This output can now be used to study how group composition (considering sex, age, and breeding status) influences group reproductive success

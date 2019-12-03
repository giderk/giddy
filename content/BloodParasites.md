# Primate Blood Parasite Dynamics
PCR and blood smear infection data
|Infections by Age Class|Ind. Changes Infection Status|
|---|---|
|![sample](https://gideonerkenswick.files.wordpress.com/2019/12/bloodparasite_ind.changes.jpg)|![sample](https://gideonerkenswick.files.wordpress.com/2019/12/bloodparasite_ageclass.jpg)|

#### Overview
* Import and merge PCR and blood smear parasite data (one record per individual per parasite per year)
* Create unified individual infection status variable
* Rearrange data to show individual infection history by parasite

```rscript
setwd("~/") # set the working director

## Required Packages ##
require(xlsx) # use read.csv or read.table
require(plyr)
require(reshape2)
require(ggplot2)
```
### Import datasets
```rscript
#load datasets that have been exported from AnimalBase (SQL database)
Smear<-read.xlsx("Smear_12.9.15.xlsx",1)
PCR<-read.xlsx("PCR_2.27.16.xlsx",1)

#Convert variables to factors
PCR$AnID<-as.factor(PCR$AnID)
PCR$FTA<-as.factor(PCR$FTA)
PCR$PCR_ID<-as.factor(PCR$PCR_ID)

Smear$AnID<-as.factor(Smear$AnID)
Smear$ParasiteID<-as.factor(Smear$ParasiteID)
Smear$ExamID<-as.factor(Smear$ExamID)

#extract year to new column
PCR$year<-as.factor(format(PCR$CaptureDate, "%Y"))
Smear$year<-as.factor(format(Smear$CaptureDate, "%Y"))

#remove PCR records that are not informative
levels(PCR$Protocol)
PCR2<-subset(PCR,PCR$Exclude!=T)
PCR2<-subset(PCR2,PCR2$Protocol!="PrimateSexing" & PCR2$Protocol!="Tryp1")
PCR2<-droplevels(PCR2)

#check for duplicate entries per animal per year in the smear dataset
Smear<-arrange(Smear,AnID)
dups<-subset(Smear, select=c(year, AnID, Genus))
dups[duplicated(dups),]

#check that no negative control results are present (these don't have an AnID)
PCR2[is.na(PCR2$AnID==F),]
```
### Organizing PCR and smear datasets, then join them
```rscript

#Normalize names for Tryp and Filarid screening
levels(PCR2$Protocol)
PCR2$Infection<-ifelse(PCR2$Protocol=="Fil.ITS1"|PCR2$Protocol=="Fil.COIint1"|
                         PCR2$Protocol=="Fil.CO1int2","filarid","tryp")

#recast the data into wide format comparing pos, neg, and unclear events
#for each individual by year by parasites
PCRsum<-dcast(PCR2,formula = Species+year+AnID+Sex+Group.Name+Name+FTA+Infection~Outcome)

#separate data into different data frames for tryp and filarid infections
temp<-split(PCRsum,PCRsum$Infection) #puts separate data sets in a list
PCRfil<-as.data.frame(temp[[1]]) # removes from list
PCRtryp<-as.data.frame(temp[[2]]) #remove from list

#assign a positive, negative or unclear status based on PCR results for filarids.
# at this point we are conservatively saying that anything that was not clearly positive
# based on one PCR screening, is negative
PCRfil$PCRfil<-ifelse(PCRfil$Pos>1,"Positive",
                      ifelse(PCRfil$Pos==1 & PCRfil$Neg==0 & PCRfil$Unclear==0,"Positive",
                             ifelse(PCRfil$Pos==1 & PCRfil$Neg==0 & PCRfil$Unclear>0,"Positive",
                                    ifelse(PCRfil$Pos==1 & PCRfil$Neg>0 & PCRfil$Unclear==0,"Unclear1",
                                           ifelse(PCRfil$Pos==1 & PCRfil$Neg>0 & PCRfil$Unclear>0,"Unclear1",
                                                  ifelse(PCRfil$Pos==0 & PCRfil$Neg>0 & PCRfil$Unclear>0,"Unclear2","Negative"))))))

# assign a positive, negative or unclear status based on PCR results for Tryp.
# We only count something as positive if we had two positive tests, otherwise it is unclear
# or negative (no positives)
PCRtryp$PCRtryp<-ifelse(PCRtryp$Pos>=2,"Positive",
                        ifelse(PCRtryp$Neg>=1&PCRtryp$Pos==1,"Unclear","Negative"))

# subset to variables of interest only, then merge PCR-based filarid and tryp infection data back together
PCRfil2<-subset(PCRfil,select=c("Species","year","AnID","FTA","PCRfil"))

# combine the data frames of different PCR infections by the following ID variables
PCRcomb<-merge(PCRtryp,PCRfil2, by=c("Species","year","AnID","FTA"))

#Remove variables from merged data set that are no needed
PCRcomb<-subset(PCRcomb,select=-c(8:11))

# fix group membership column name
colnames(PCRcomb)[6]<-"Group"


#### Organize Smear Data before merge

#recast smear data into wide format
Smear2<-dcast(Smear, formula=Species+Sex+year+AnID+Group+Name+
                ExamID+Diff+Fields400x+TotWBC+Neut+Lymph+Mono+Eosin+Baso~
                Genus, value.var="Count")


#### Merge PCR infection with smear data by common ID variables

Blood<-merge(Smear2,PCRcomb, by=c("Species","Sex","year","AnID","Group","Name"))
````

>Now that we have smear and PCR infection data in one data frame, we can use both types  of information to give each animal a single relaible infection status. The following infection status codes will apply...

* score of 1 means clear infection based on PCR or smear
* score of 0 means failure to show infection by PCR and smear
* score of 2 means PCR can't confirm the infection, but smear shows infection
* the one condition that we have to check manually is when PCR is positive for a filarid
* but the smear is negative.  Because our marker pics up both mansonella and dipetalonema
* we can't automatically determine which of the 2 filarid parasites are present.
* for simps we know it has to be Dip, but for fuscis it could be either or both

```rscript
#ifelse arguements will not work unless NAs are removed from all infection variables
Blood$Mansonella[is.na(Blood$Mansonella)]<-0
Blood$Dipetalonema[is.na(Blood$Dipetalonema)]<-0
Blood$Trypanosoma[is.na(Blood$Trypanosoma)]<-0



Blood$M.mariae<-ifelse(Blood$PCRfil=="Positive"&Blood$Mansonella>0,1,
                       ifelse(Blood$PCRfil=="Negative"&Blood$Mansonella>0,2,
                              ifelse(Blood$PCRfil=="Unclear1"&Blood$Mansonella>0,2,
                                     ifelse(Blood$PCRfil=="Unclear2"&Blood$Mansonella>0,2,0))))

Blood$D.gra<-ifelse(Blood$PCRfil=="Positive"&Blood$Dipetalonema>0,1,
                    ifelse(Blood$PCRfil=="Negative"&Blood$Dipetalonema>0,2,
                           ifelse(Blood$PCRfil=="Unclear1"&Blood$Dipetalonema>0,2,
                                  ifelse(Blood$PCRfil=="Unclear2"&Blood$Dipetalonema>0,2,0))))

Blood$T.min<-ifelse(Blood$PCRtryp=="Positive",1,
                    ifelse(Blood$PCRtryp=="Negative"&Blood$Trypanosoma>=1,2,
                           ifelse(Blood$PCRtryp=="Unclear"&Blood$Trypanosoma>=1,2,0)))

## There are two individuals without smear infections but are positive for Dipetalonema
# Identify those individuals
Blood[Blood$PCRfil=="Positive" & Blood$Dipetalonema==0 & Blood$Mansonella==0,]

# correct their Dipe infection by querrying on FTA, a number 3 indicates a PCR positive
# smear negative individual for a given parasite, whereas a 2 indicates a smear positive
# PCR negative, whereas a 1 and 0 indicate that smear and pcr detection are in agreement

Blood[Blood$FTA==119,]$D.gra<-3
Blood[Blood$FTA==131,]$D.gra<-3

# export the working data frame
write.xlsx(Blood,"Blood.xlsx")
```


#### Combine infection data with individual breeding status

```rscript
# load the most current breeding status data
Morph<-read.xlsx("BreedingStatus2.19.16.xlsx",1)

# Subset only 2012-2014 from the morph dataset, since we only currently have infection
# data from 2012 through 2014 (three years)
Morph<-subset(Morph, Morph$year!="2015" & Morph$year!="2009" &
                Morph$year!="2010" & Morph$year!="2011")
Morph<-droplevels(Morph) # this removes levels that are no longer present in data

#Select only the variables we need from the Morph dataframe
Morph1<-subset(Morph,select=c("Species","Sex","year","AnID",
                              "Group","Name","Age_class","Date",
                              "nWeight","br.dfa","ExamID"))

# There are different numbers of records in the blood and morph datasets, this is
# because we did not collect blood from absolutely everyone trapped. Merge the
# the two data frames but do not include individuasl that don't appear in both,
# which is called an inner join. So in the code below we need to specify
# all.x = True
data<-merge(Blood,Morph1, by=c("Species","Sex","year","AnID","Group","Name"), all.x=T)

# check for blood infection records that did not match with an individual from
# the morph dataset, by querrying on ExamID.x.  If a record did match it would have NAs
# for all the blood parasite exams
data[is.na(data$ExamID.x),]
data<-arrange(data,ExamID.x, AnID)

# Check for duplicates in the data data frame, there should be none
data<-arrange(data,AnID)
dups<-subset(data, select=c(year, AnID))
dups[duplicated(dups),]

########## create working data set ############
parasite<-data

# write a file so we can check how many mismatches between PCR and smear detection there are
write.xlsx(data,"parasitexx.xlsx")

# the presence absence infections currently have 0,1,2,3 (see above for detail)
# need to convert all the 2s and 3s into 1s before we can do data analysis, but don't include
# this step if wanting descriptive stats on how many conflicting results
parasite$M.mariae[parasite$M.mariae==2]<-1
parasite$D.gra[parasite$D.gra==2]<-1
parasite$T.min[parasite$T.min==2]<-1
parasite$M.mariae[parasite$M.mariae==3]<-1
parasite$D.gra[parasite$D.gra==3]<-1
parasite$T.min[parasite$T.min==3]<-1

# Calculate parasite species richness as a new variable
parasite<-ddply(parasite, .(Species), transform, PSR= M.mariae+D.gra+T.min)

# save the file
write.xlsx(parasite,"parasite.xlsx") # does not include prior infection history
```

### Preparing data to analyze prior infection history
```rscript
# to look at history of infection we need to re-organize the data so that
# we only have records from individuals that span multiple seasons, namely
# 2012-2013 and 2013-2014
names(parasite)
# take all out but species,year, name, AnID, sex, and the three parasite presence/
# absence categories
para<-subset(parasite,select=-c(5,7:22,26:32))
names(para)

# instead of 0s and 1s, let's work with categorical names
para$M.mariae<-ifelse(para$M.mariae==1,"Pos","Neg")
para$D.gra<-ifelse(para$D.gra==1,"Pos","Neg")
para$T.min<-ifelse(para$T.min==1,"Pos","Neg")

# melt the data so that parasite infections are in long format
para1<-melt(para, id.vars=c("Species","Sex","AnID","Name","year"),na.rm=F)
names(para1)

# cast the data back into wide format by year*parasite, this will give NAs for individuals
# that are not sampled in a particular year
para2<-dcast(para1,formula = Species+Sex+AnID+Name~variable+year, value.var="value")

# with the original parasite data frame, subset by year 2013 and 2014
para2013<-subset(parasite,parasite$year=="2013");para2013<-droplevels(para2013)
para2014<-subset(parasite,parasite$year=="2014");para2013<-droplevels(para2013)

names(para2)
# with the newly recast para2 dataframe, select all the identity variables and
# infections from a single year, everything else needs to dissapear
par2012<-subset(para2, select=c(1:4,5,8,11)) ; par2012<-droplevels(par2012)

# remove all the rows that didn't have any data in 2012, this means that the animal
# was not sampled that year
par2012<-par2012[complete.cases(par2012),]

# merge the 2012 infections data with the 2013 parasite data, and
# now we have 2013 infections with prior infection history, since only matching
# individuals will be retained
P12.13<-merge(para2013,par2012, by=c("Species","Sex","AnID","Name"))

names(para2)
# repeat the same steps to get the 2013-2014 section of data
par2013<-subset(para2, select=c(1:4,6,9,12)) ; par2013<-droplevels(par2013)
par2013<-par2013[complete.cases(par2013),]
P13.14<-merge(para2014,par2013, by=c("Species","Sex","AnID","Name"))
names(P13.14)

# now we will change the names of the prior infection variables to match in
# in the p13.14 and p12.13 data frames, so we can have one large data set in which
# all animals have known prior infections status.

# Instead of having to write all the variable by scratch, we can take the ones that
# aren't changed and simply add the new, altered names.
newnames<-c(names(P12.13)[1:31],c("M.mar_prior","D.gra_prior","T.min_prior"))
names(P12.13)<-newnames ; names(P13.14)<-newnames

# put these two data frames together with rbind
para.priors<-rbind(P12.13,P13.14)

## now write the data frame to an excel file for safe keeping and this is ready
# for analysis
write.xlsx(para.priors,"para.priors.xlsx")
```
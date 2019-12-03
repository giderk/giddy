# Using group composition to model tamarin  reproductive output

## Overview
* Convert tamarin individual data into summarized group data
* Manually corrections (account for individuals that were present but not trapped)
* Data Inspection
* Statistical Modeling

```rscript

setwd("~/") # set workind dir

## Required Packages ##
require(xlsx) # use read.csv or read.table instead
require(plyr)
require(reshape2)
require(ggplot2)
library(lmtest)
library(car)
library(GGally)
library(lme4)
library(MASS)

#load data

##########These have been done already in the Breeders excel file, but confirm always#
# add one to FC group size 2011, because GRC was not trapped (add her as P_female also, etc.)
# add one to MI4 2014 because LRC was not trapped (add her as  P_Female also, etc)
# add one primary male,YBG, to Queen Bees 2015, etc
# add one F_secondary,RPBL, to JI6 2015 for , etc
# add one F_primary to JI6 2012 for RPG/RRC, etc
# add one F_primary to JF 2012 RPO, etc.
# add another to JI 2013 for RPBL, who is a subadult, secondary status
# add anther to JF 2011 for RBY, who is an adult, primary (trapped in 2010 and 2012)

status<-read.xlsx("BreedingStatus.xlsx",1) #get the most current br.status file

#subset status to only the columns of interest
data<-subset(status, select=c("AnID","Date","Species","Sex","Group","Name","Age_class",
                              "year","br.dfa"))

# change all disperser categories to secondary
data$F.young<-ifelse(data$Sex=="F" & data$Age_class!="Juvenile" &
                       data$br.dfa=="None","Fy","NA")
data$M.young<-ifelse(data$Sex=="M" & data$Age_class!="Juvenile" &
                       data$br.dfa=="None","My","NA")

# recast data, summarizing various aspects of group composition
data1<-dcast(data, Species+Group+year~Age_class) # age class
data2<-dcast(data, Species+Group+year~Sex+br.dfa) # breeding status by sex
data3<-dcast(data, Species+Group+year~br.dfa) # breeding status, sexes combined
data4<-ddply(data,.(Species,Group,year), nrow) # group size
data5<-dcast(data, Species+Group+year~Sex) # males and females
data6<-dcast(data, Species+Group+year~M.young) # number of young males
data7<-dcast(data, Species+Group+year~F.young) # number of young females
data8<-dcast(data, Species+Group+year~Age_class+Sex)

# merge summarized data
rm(comp2) # remove any prior merged file (incase this step is repeated)
comp2<-merge(data4,data1, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data5, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data2, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data3, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data6, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data7, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data8, by=c("Species","Group","year"), all=T)

colnames(comp2)[4]<-"grp.size" # correct group size column name
names(comp2)
comp2a<-subset(comp2,select=-c(19,21))


# calculating additional model parameters
comp2a<-ddply(comp2a,.(Species),mutate,rSex=round(M/F,2),grp.size.new=grp.size-Juvenile,
              Breeders=Secondary+Primary,pBreeder=round(Breeders/grp.size.new,2),
              F.new=F_Primary+F_Secondary+Fy, M.new=M_Primary+M_Secondary)

comp2a<-subset(comp2a, comp2a$Species!="CBRU") #remove the callicebus records


write.xlsx(comp2a,"grp.comp.xlsx") # write the data frame to a file 
```
## make the following manual corrections to the file that was just made
* before importing grp.comp.xlsx do the following
* confirm that all 2009 records are removed
* correct FC row 2010 to have 3 (GRC,GPG,GBR) adults, 2 P_females,
* 1 P_Male, 2 Female juvs (GPO and GPBL) the numbers are off because we retrapped the juvs many times in 2010
* also, add one primary breeding female to AR6 2013

```rscript

comp<-read.xlsx("grp.comp.xlsx",1) # import corrected file

# clean-up column headings
comp$grp.size<-comp$grp.size.new
comp$F<-comp$F.new
comp$M<-comp$M.new
names(comp)
comp<-subset(comp, select=-c(28,31,32))

# caclulate proportions since group size correlated with counts of everything
comp<-ddply(comp,.(Species), mutate,pF=F/grp.size,pM=M/grp.size,
            pFp=F_Primary/grp.size,pMp=M_Primary/grp.size,
            pFs=F_Secondary/grp.size,pMs=M_Secondary/grp.size,rSex=M/F,
            helpers=M_Primary+M_Secondary+F_Secondary-1,pHelp=helpers/grp.size,
            pAM=Adult_M/grp.size,pAF=Adult_F/grp.size,pSM=Sub.adult_M/grp.size,
            pSF=Sub.adult_F/grp.size)

comp<-subset(comp,comp$Group!="Loners") # remove loners that don't have a group

# remove group MI4 2012 because we didn't actually trap this group, we caught 2 individuals that were attempting to disperse
# Remove angels, not really a group
# Remove Orange Emps 2012, not really a group
# remove bees from 2013 becuase it was a newly formed group
comp<-subset(comp, comp$Group!="Angels" & comp$Group!="");comp<-droplevels(comp)
comp<-subset(comp,comp$Group!="OrangeEmps" | comp$year!="2012");comp<-droplevels(comp)
comp<-subset(comp,comp$Group!="MI4" | comp$year!="2012");comp<-droplevels(comp)
comp<-comp[comp$Group!="Bees"|comp$year!=2013,]

# data frame is ready for modeling
write.xlsx(comp, "comp.xlsx")
```
## Data Inspection
```rscript
comp<-read.xlsx("comp.xlsx",1) #load data

###### Data exploration briefly

ggpairs(comp[,c("Juvenile","grp.size","pFp","pMp","pFs","pMs","M","F","pHelp")]) # pairwise correlation matrix

ggpairs(comp[,c("Juvenile","M","F","M_Primary","F_Primary","M_Secondary", # pairwise correlation matrix
                "F_Secondary","pMp","pFp")])

# map a categorical variables onto the plot matrix
scatterplotMatrix(~Juvenile+grp.size+pMp+F_Primary+pMp+pFs+pF+pM|Species,data=comp,
                  main="Social Groups")

# use variations of coplots 
coplot(Juvenile ~ pPF.M|grp.size, data = comp, panel = panel.smooth,
       number = 7, xlab = c("pPF", "grp"),row=1)

# inspect data by years, groups, sex, species
lapply(comp[, c("Juvenile", "Species", "year","Group")], table)

# looking at crosstables of groups and number of juveniles (large screen space required)
ftable(xtabs(~ Species + Juvenile + Group, data = comp))

# check for overdispersion (variance to mean ratio)
rbind(colMeans(comp[5:28]),apply(comp[5:28], 2, sd)^2)

# the marginal category that groups all the data together shows a clear age effect, and you can also
# see that no group appears to be biasing the sample too much.
ggplot(comp,aes(x=Juvenile,y=grp.size))+
  geom_boxplot(size=.75)+
  geom_jitter(alpha=.5)+
  facet_grid(Species~Group, margins=T)

# distribution of response variable
hist(comp[comp$Species=="SFUS",]$Juvenile,breaks=5)
hist(comp[comp$Species=="SIMP",]$Juvenile,breaks=5)
```
## Checking for random structure that needs to be

```rscript

# subset the species
compe<-subset(comp, comp$Species=="SIMP");compe<-droplevels(compe)
compf<-subset(comp, comp$Species=="SFUS");compf<-droplevels(compf)

# checking for random structure
a<-glmer(formula=Juvenile~grp.size+pFp+pMp+pFs+pMs+pHelp+(1|Group)+(1|year), data=comp,
      family="poisson");summary(a)

a<-glmer(formula=Juvenile~grp.size+pFp+pMp+pFs+pMs+pHelp+(1|Group), data=comp,
         family="poisson");summary(a)

a<-glmer(formula=Juvenile~grp.size+pFp+pMp+pFs+pMs+pHelp+(1|Stability), data=comp,
         family="poisson");summary(a)

# group and year appear to have 0 affect
```
## Constructing Statistical Models
```rscript
# both species together using glm with a poisson distribution
fitHi<-glm(formula=Juvenile~pFp+pMp+pMs+pFs+grp.size+Species,
          data=comp,family=poisson);fitHi

summary(fitHi)
Anova(fitHi)

# finding the best fit mopdel by step-wise term deletion
fittest <- stepAIC(fitHi,
                   scope=list(upper=~pFp+pMp+pFs+pMs+grp.size+Species,
                              lower=~1),trace=T, direction="backward")
fittest$anova

# best fit model
fitb<-glm(formula=Juvenile~pFp+pMp+grp.size,
          data=comp,family=poisson);summary(fitb)
Anova(fitb)

fitc<-glm(formula=Juvenile~pFp+grp.size,
          data=comp,family=poisson);summary(fitc)

lrtest(fitb,fitc)
Anova(fitc)
drop1(fitb,test="Chisq")

# inspecting residuals and verifying model fit
plot(fitc)

qplot(fitb$fitted.values,comp$Juvenile)+stat_smooth(method="lm",se=T)
qplot(comp$Juvenile,residuals(fitb,"pearson"))+stat_smooth(method="lm",se=T)
qplot(fitted.values(fitb),residuals(fitb))+stat_smooth(method="lm",se=T)

qplot(comp$grp.size,residuals(fitb))+stat_smooth(method="lm",se=T)
qplot(comp$pMp,residuals(fitb))+stat_smooth(method="lm",se=T)
qplot(comp$pMp,residuals(fitb))+stat_smooth(method="lm",se=T)
qplot(comp$pFp,residuals(fitb))+stat_smooth(method="lm",se=T)
qplot(comp$pFs,residuals(fitb))+stat_smooth(method="lm",se=T)
qplot(comp$pF,residuals(fitb))+stat_smooth(method="lm",se=T)
qplot(comp$pM,residuals(fitb))+stat_smooth(method="lm",se=T)
qplot(comp$Species,residuals(fitb))+stat_smooth(method="lm",se=T)

###### GLM of just Fuscis

fitHi<-glm(formula=Juvenile~grp.size+pFp+pM+pFs,
           data=compf,family=poisson);fitHi
summary(fitHi)

fittest <- stepAIC(fitHi,
                   scope=list(upper=~grp.size+pFp+pM+pFs,
                              lower=~1),trace=T, direction="backward")
fittest$anova

# best fit model
fitf<-glm(formula=Juvenile ~ grp.size + pFp,
          data=compf,family=poisson);summary(fitf)
plot(fitf)

#### GLM of just emps

fitHi<-glm(formula=Juvenile~grp.size+pFp+pM+pFs,
           data=compe,family=poisson);fitHi
summary(fitHi)

fittest <- stepAIC(fitHi,
                   scope=list(upper=~grp.size+pFp+pM+pFs,
                              lower=~1),trace=T, direction="backward")
fittest$anova

# best fit model
fite<-glm(formula=Juvenile ~ pFp,
          data=compe,family=poisson);summary(fite)
plot(fite)
```
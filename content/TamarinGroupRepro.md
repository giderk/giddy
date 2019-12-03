# Using group composition to model tamarin reproductive output

|Regression Probability Curves|Pairwise Correlation Matrix|
|---|---|
|![example](https://gideonerkenswick.files.wordpress.com/2019/12/mulinomial_regression_curves.jpg)|![sample](https://gideonerkenswick.files.wordpress.com/2019/12/pairwise_correlations.jpg)|

## Overview
* Convert tamarin individual data into summarized group data
* Manual correction (account for individuals that were present but not trapped)
* Data Inspection
* Statistical Modeling
* Plotting predicted outcomes

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
library(grid)
library(coin)
require(nnet)

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
### Make the following manual corrections to the file that was just made
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
### Data Inspection
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
### Checking for random structure before modeling data

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
#### GLMs
```rscript
### both species together using glm with a poisson distribution
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


###### GLM of just emps

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
### Logistic regression

> the main problem is that we have to boil the response variable down to just 0 and 1 (two categories) when in fact there are 4 ordered categories (0 juvs, 1 juv, 2 juvs, 3 juvs) if we go with this option, we should divide up the response variable to first look at the difference between groups that have 0 juvs and those that have any juvs, and then we should analysize a subset of data to see what explains differences betwen groups with 1 juv and 2 or more juvs
```rscript

# comparing zero juvs to 1 juv
comp3<-subset(comp, comp$Juvenile<2)
comp3$Juv<-ifelse(comp3$Juvenile>0,1,0)

L1 <- glm(Juv~pFp+pMp+pFs+pMs+grp.size+Species, 
                    data = comp3, family = binomial)
summary(L1) # take out pMs

L2 <- glm(Juv~pFp+pMp+pFs+grp.size+Species, 
          data = comp3, family = binomial)
summary(L2) # take out pFs

L3 <- glm(Juv~pFp+pMp+grp.size+Species, 
          data = comp3, family = binomial)
summary(L3)  # converged on best model, sample size not large enough,
# just do a single species

L4 <- glm(Juv~pFp+pMp+grp.size, 
          data = comp3, family = binomial)
summary(L4) 
# unclear what to do at this point

L5 <- glm(Juv~pMp+pFp, 
          data = comp3, family = binomial)
summary(L5) 

L6 <- glm(Juv~pFp, 
          data = comp3, family = binomial)
summary(L6) 


AIC(L1,L2,L3,L4) #L3 is the winner by AIC score
lrtest(L2,L3)

drop1(L4,test="Chisq")

# comparing 1 juv to 2 or more juvs sample size is too small, but by looking at just SFUS pFp and group size approach significance

comp2<-subset(comp, comp$Juvenile>0)
comp2$Juv3<-ifelse(comp2$Juvenile>1,1,0)


G1 <- glm(Juv3~F_Primary+F_Secondary+M_Secondary+M_Primary+grp.size+Species, 
          data = comp2, family = binomial)
summary(G1) # take out pFs

G2 <- glm(Juv3~F_Primary+F_Secondary+M_Secondary+M_Primary+grp.size, 
          data = comp2, family = binomial)
summary(G2) # take out grp.size

G3 <- glm(Juv3~F_Primary+M_Secondary+M_Primary+grp.size, 
          data = comp2, family = binomial)
summary(G3) # pFp and group size approach significance but sample size
# too small

G4 <- glm(Juv3~F_Primary+M_Primary+grp.size, 
          data = comp2, family = binomial)
summary(G4) #

G5 <- glm(Juv3~F_Primary+M_Primary, 
          data = comp2, family = binomial)
summary(G5) #

G6 <- glm(Juv3~F_Primary, 
          data = comp2, family = binomial)
summary(G6) #


lrtest(G1,G4)
AIC(G1,G2,G3,G4)

drop1(G4,test="Chisq")

# Looking at zero juvs to any juvs
# emps don't have the sample size to test separately
# fuscis barley have the sample size
# pFp and group size are very close to signicant, pFp alone is sig.

t1<-ggplot(comp,aes(y=pMp,x=as.factor(Juv)))+geom_boxplot()+
  theme_bw(base_size = 12)+
  scale_y_continuous(limits=c(0:1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Proportion Primary Males")+xlab("Infants")

t2<-ggplot(comp,aes(y=pFp,x=as.factor(Juv)))+geom_boxplot()+
  theme_bw(base_size = 12)+
  scale_y_continuous(limits=c(0:1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Proportion Primary Females")+xlab("Infants")

t3<-ggplot(comp,aes(y=grp.size,x=as.factor(Juv)))+geom_boxplot()+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Group size")+xlab("Infants")

grid.arrange(t1,t2,t3,nrow=1)

wilcox_test(pFp~as.factor(Juv),data=comp, distribution="exact")
wilcox_test(pMp~as.factor(Juv),data=comp, distribution="exact")
wilcox_test(grp.size~as.factor(Juv),data=comp, distribution="exact")

comp$Juv<-ifelse(comp$Juvenile>0,1,0)

K1 <- glm(Juv~pFp+pMp+pFs+pMs+grp.size+Species, 
          data = comp, family = binomial)
summary(K1)
# take out pMs
Anova(K1)

K2 <- glm(Juv~pFp+pMp+pFs+grp.size+Species, 
          data = comp, family = binomial)
summary(K2) # take out pFs
Anova(K2)

K3 <- glm(Juv~pFp+pMp+grp.size+Species, 
          data = comp, family = binomial)
summary(K3) # take out species

K4 <- glm(Juv~pFp+pMp+grp.size, 
          data = comp, family = binomial)
summary(K4) # best fit

Anova(K4)

K5 <- glm(Juv~pFp+pMp, 
          data = comp, family = binomial)
summary(K5)

Anova(K5)

K6 <- glm(Juv~pMp, 
          data = comp, family = binomial)
summary(K6)

lrtest(K4,K5) 

AIC(K4,K5,K6)

drop1(K4,test="Chisq")
```
### Multinomial Logistic Regression

> this may be a worthwhile analysis to pursue, becaue we can use all categories of juveniles at once, however need to test for violations of the model. http://www.ats.ucla.edu/stat/r/dae/mlogit.htm

```rscript
comp<-read.xlsx("comp.xlsx",1)

comp$Juv2<-comp$Juvenile
comp$Juv2[comp$Juv2==3]<-2
comp$Juv2<-as.factor(comp$Juv2);levels(comp$Juv2)
comp$Juv2 <- relevel(comp$Juv2, ref = "0");levels(comp$Juv2)

mn1 <- multinom(Juv2 ~ pFp+pMp+pFs+pMs+grp.size+Species, data = comp);mn1
summary(mn1, digits=4)
z1 <- summary(mn1)$coefficients/summary(mn1)$standard.errors;z1
p1 <- (1 - pnorm(abs(z1), 0, 1)) * 2;p1
#aic=132; remove species

mn2 <- multinom(Juv2 ~ pFp+pMp+pFs+pMs+grp.size, data = comp)
summary(mn2, digits=4)
z2 <- summary(mn2)$coefficients/summary(mn2)$standard.errors;z2
p2 <- (1 - pnorm(abs(z2), 0, 1)) * 2;p2
#aic=128 remove pFs

mn3 <- multinom(Juv2 ~ pFp+pMp+pFs+grp.size, data = comp)
summary(mn3, digits=4)
z3 <- summary(mn3)$coefficients/summary(mn3)$standard.errors;z3
p3 <- (1 - pnorm(abs(z3), 0, 1)) * 2;p3
#aic=128 sig = all , pFs is not significant

mn4 <- multinom(Juv2 ~ pFp+pMp+grp.size, data = comp)
summary(mn4, digits=4)
z4 <- summary(mn4)$coefficients/summary(mn4)$standard.errors;z4
p4 <- (1 - pnorm(abs(z4), 0, 1)) * 2;p4
#aic=128 sig = all , pMp is not significant

AIC(mn4,mn3,mn2,mn1)

Anova(mn4,type="III")

dropterm(mn4,test="Chisq")

lrtest(mn1,mn4) # according to lr test the diff is not significant, so we did not need to include pMp
```
### plotting the predicted contribution of a given variable ...
```rscript
# create a new datafame (the new data set) with the exact same predictor variables as the model.
# for your variable of interest fill it with a string of data from the minimum to the maximum
# value.  for the other variables insert the average value for each row of data.
grp <- data.frame(grp.size=seq(2,7,length.out=63), pFp = mean(comp$pFp),pMp=mean(comp$pMp))

# create probabilities (predictions) from the model
predict(mn3, newdata = grp, "probs")

# attach the predictions to the new dataset
pp.grp<-cbind(grp,predict(mn3, newdata = grp, "probs",se=T));head(pp.grp)
pp.grp<-subset(pp.grp,select=);head(pp.grp)# The minus 7 was to remove the one instance of having a group with more than 2 juveniles, we have now collapsed it to 0,1, or 2+ juvs

#plot the variable of interest on the x-axis and the probability on the y axis, and color code
# the lines by whatever variable you want, in this case probably number of juveniles

Pgrp <- melt(pp.grp, id.vars = c("grp.size","pFp","pMp"),
             value.name = "probability");head(Pgrp)
names(Pgrp)[4]<-"Offspring"
head(Pgrp)  # view first few rows

# plot the probabilities of a given group size based on the number of offspring present
g1<-ggplot(Pgrp, aes(x = grp.size, y = probability, colour = Offspring))+
  geom_line()+ggtitle("Group size probabilities by number of offspring")+
  labs(color='Offspring')+xlab("Group Size");g1

#repeat with pFp
pFp <- data.frame(pFp=seq(0,.6,length.out=63), grp.size = mean(comp$grp.size),pMp=mean(comp$pMp))
predict(mn3, newdata = pFp, "probs")
pp.pFp<-cbind(pFp,predict(mn3, newdata = pFp, "probs",se=T));head(pp.pFp)
pp.pFp<-subset(pp.pFp,select=);head(pp.pFp)
PpFp <- melt(pp.pFp, id.vars = c("pFp","grp.size","pMp"),
             value.name = "probability");head(PpFp)
names(PpFp)[4]<-"Offspring";head(PpFp)
g2<-ggplot(PpFp, aes(x = pFp, y = probability, colour = Offspring))+
  geom_line()+ggtitle("Proportion of Primary Females by Number of Infants")+
  labs(color='Offspring')+xlab("Proportion of Female Primary Breeders");g2

#repeat with pMp
pMp <- data.frame(pMp=seq(0,.8,length.out=63), grp.size = mean(comp$grp.size),pFp=mean(comp$pFp))
predict(mn3, newdata = pMp, "probs")
pp.pMp<-cbind(pMp,predict(mn3, newdata = pMp, "probs",se=T));head(pp.pMp)
pp.pMp<-subset(pp.pMp,select=);head(pp.pMp)
PpMp <- melt(pp.pMp, id.vars = c("pMp","grp.size","pFp"),
             value.name = "probability");names(PpMp)[4]<-"Offspring"
head(PpMp)
g3<-ggplot(PpMp, aes(x = pMp, y = probability, colour = Offspring))+
  geom_line()+ggtitle("Proportion of Primary Males by Number of Infants")+
  labs(color='Offspring')+xlab("Proportion of Male Primary Breeders");g3


Pgraphs<-rbind.fill(melt(PpFp,id.vars=c("grp.size","Offspring","probability","pMp"),value.name="prop"),
                    melt(Pgrp,id.vars=c("pFp","Offspring","probability","pMp"),value.name="prop"),
                    melt(PpMp,id.vars=c("pFp","Offspring","probability","grp.size"),value.name="prop"))
head(Pgraphs)
levels(Pgraphs$variable)<-c("Prop. Primary Female","Group Size","Prop. Primary Male")

ggplot(Pgraphs, aes(x = prop, y = probability, color = Offspring))+
  geom_line(lwd = 1) +
  theme(axis.title.x=element_blank(), panel.grid.major=element_blank(),
        text=element_text(size=20))+
  facet_grid(.~variable,scale="free")+
  ggtitle("Results From Multinomial Logistic Regression")
  ```


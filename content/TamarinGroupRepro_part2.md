# Modeling group composition and breeding status
setwd("~/Documents/Work/Research/Data/Analysis/BrStatus&GroupComp")

library(xlsx)
library(plyr)
library(ggplot2)
library(lmtest)
library(car)
library(reshape2)

#load data
comp<-read.xlsx("comp.xlsx",1)

###### Data exploration briefly

#inspect the data with plots and correlation matrix
library(GGally)
ggpairs(comp[,c("Juvenile","grp.size","pFp","pMp","pFs","pMs","M","F","pHelp")])

test1<-subset(test,test$Juvenile==2)

ggpairs(comp[,c("Juvenile","M","F","M_Primary","F_Primary","M_Secondary",
                "F_Secondary","pMp","pFp")])

#this allows us to map a particular categorical variable onto the plot matrix
library(car)
scatterplotMatrix(~Juvenile+grp.size+pMp+F_Primary+pMp+pFs+pF+pM|Species,data=comp,
                  main="Social Groups")

# use variations of coplots  Coplots
coplot(Juvenile ~ pPF.M|grp.size, data = comp, panel = panel.smooth,
       number = 7, xlab = c("pPF", "grp"),row=1)

# take note of thigs that are strongly correlated because they frequently cannot be graphed together

# let's check out the data by years, groups, sex, species
lapply(comp[, c("Juvenile", "Species", "year","Group")], table)


# looking at crosstables of groups and number of juveniles (need large screen space)
ftable(xtabs(~ Species + Juvenile + Group, data = comp))


#looking at variance to mean ratios of all numeric variables to check for over dispersion
rbind(colMeans(comp[5:28]),apply(comp[5:28], 2, sd)^2)

# graphing the data (need to do lots of these, they are very informative).  You can see that
# the marginal category that groups al the data together shows a clear age effect, and you can also
# see that no group appear to be biasing the sample too much.
ggplot(comp,aes(x=Juvenile,y=grp.size))+
  geom_boxplot(size=.75)+
  geom_jitter(alpha=.5)+
  facet_grid(Species~Group, margins=T)

# distribution of response variable
hist(comp[comp$Species=="SFUS",]$Juvenile,breaks=5)
hist(comp[comp$Species=="SIMP",]$Juvenile,breaks=5)

###### generalized linear model with a poisson distribution
# First thing we need to do is find if there is any random structure that needs to be 
# incorporated in to the models, it would be ideal to graph the variables that we 
# do not want to influence the model, such as group and species

#subset the species
compe<-subset(comp, comp$Species=="SIMP");compe<-droplevels(compe)
compf<-subset(comp, comp$Species=="SFUS");compf<-droplevels(compf)

# load the glm package
library(lme4)

#checking for random structure
a<-glmer(formula=Juvenile~grp.size+pFp+pMp+pFs+pMs+pHelp+(1|Group)+(1|year), data=comp,
      family="poisson");summary(a)

a<-glmer(formula=Juvenile~grp.size+pFp+pMp+pFs+pMs+pHelp+(1|Group), data=comp,
         family="poisson");summary(a)

a<-glmer(formula=Juvenile~grp.size+pFp+pMp+pFs+pMs+pHelp+(1|Stability), data=comp,
         family="poisson");summary(a)

#group and year appear to have 0 affect


######Graphing both species together using glm witha  poisson distribution
fitHi<-glm(formula=Juvenile~pFp+pMp+pMs+pFs+grp.size+Species,
          data=comp,family=poisson);fitHi
summary(fitHi)

Anova(fitHi)

library(MASS)
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

# inspecting residuals and verify model fit
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


################ Logistic regression ###################################
# the main problem is that we have to boild the response variable down to just 0 and 1 (two,
# categories) when in fact there are 4 ordered categories (0 juvs, 1 juv, 2 juvs, 3 juvs)
# if we go with this option, we should divide up the response variable to first look at the
# difference between groups that have 0 juvs and those that have any juvs, and then
# we should analysize a subset of data to see what explains differences betwen groups with 1
# juv and 2 or more juvs


### Looking at comparing zero juvs to 1 juv

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

# looking at comparing 1 juv to 2 or more juvs
# sample size is too small but by looking at just SFUS
# pFp and group size approach significance

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

# looking at zero juvs to any juvs
# emps don't have the sample size to test separately
# fuscis barley have the sample size, but pFp and group
# size are very close to signicant, pFp alone is sig.

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


library(grid)
library(coin)


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


################# multinomial logistic regression###############
# this may be a worthwhile analysis to pursue, becaue we can use all categories of juveniles
# at once, however i'm not perfectly sure how to test for violations of the model
# http://www.ats.ucla.edu/stat/r/dae/mlogit.htm
require(nnet)
require(MASS)

comp<-read.xlsx("comp.xlsx",1)

#comp<-comp[comp$Stability=="old",]

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

# if we want to plot out the contribution of a given predictor variable we can use the following
# steps, or we can model reproduction based on a new dataset

# create a new datafame (the new data set) with the exact same predictor variables as the model.
# for your variable of interest fill it with a string of data from the minimum to the maximum
# value.  for the other variables insert the average value for each row of data.
grp <- data.frame(grp.size=seq(2,7,length.out=63), pFp = mean(comp$pFp),pMp=mean(comp$pMp))
# this data set you will use to create probabilities (predictions) from the model
predict(mn3, newdata = grp, "probs")
# attach the predictions to the new dataset
pp.grp<-cbind(grp,predict(mn3, newdata = grp, "probs",se=T));head(pp.grp)
pp.grp<-subset(pp.grp,select=);head(pp.grp)# The minus 7 was to remove the one instance of 
# having a group with more than 2 juveniles, we have now collapsed it to 0,1, or 2+ juvs

#plot the variable of interest on the x-axis and the probability on the y axis, and color code
# the lines by whatever variable you want, in this case probably number of juveniles
library(reshape2)
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


######### modeling without pFp ##############

mn1 <- multinom(Juv2 ~ grp.size+pMp+pMs+pFs+Species, data = comp);mn1
summary(mn1, digits=4)
z1 <- summary(mn1)$coefficients/summary(mn1)$standard.errors;z1
p1 <- (1 - pnorm(abs(z1), 0, 1)) * 2;p1
#aic=148; Species not significant

mn2 <- multinom(Juv2 ~ grp.size+pMs+pFs+pMp, data = comp)
summary(mn2, digits=4)
z2 <- summary(mn2)$coefficients/summary(mn2)$standard.errors;z2
p2 <- (1 - pnorm(abs(z2), 0, 1)) * 2;p2
#aic=146; remove group size

mn3 <- multinom(Juv2 ~ grp.size+pMs+pFs, data = comp)
summary(mn3, digits=4)
z3 <- summary(mn3)$coefficients/summary(mn3)$standard.errors;z3
p3 <- (1 - pnorm(abs(z3), 0, 1)) * 2;p3
#aic=144; remove pMp

mn4 <- multinom(Juv2 ~ pMs+pFs, data = comp)
summary(mn4, digits=4)
z4 <- summary(mn4)$coefficients/summary(mn4)$standard.errors;z4
p4 <- (1 - pnorm(abs(z4), 0, 1)) * 2;p4
#aic=142; remove pFs

mn5 <- multinom(Juv2 ~ pMs, data = comp)
summary(mn5, digits=4)
z5 <- summary(mn5)$coefficients/summary(mn5)$standard.errors;z5
p5 <- (1 - pnorm(abs(z5), 0, 1)) * 2;p5
#aic=141

AIC(mn5,mn4,mn3,mn2,mn1)

Anova(mn5,type="III")

dropterm(mn4,test="Chisq")

lrtest(mn2,mn5) # according to lr test the diff is not significant, so we drop all but pMs

# if we want to plot out the contribution of a given predictor variable we can use the following
# steps, or we can model reproduction based on a new dataset

pMp <- data.frame(pMp=seq(0,.8,length.out=64))
predict(mn5, newdata = pMp, "probs")
pp.pMp<-cbind(pMp,predict(mn5, newdata = pMp, "probs",se=T));head(pp.pMp)
pp.pMp<-subset(pp.pMp,select=);head(pp.pMp)
PpMp <- melt(pp.pMp, id.vars = c("pMp"),
             value.name = "probability");head(PpMp)
g2<-ggplot(PpMp, aes(x = pMp, y = probability, colour = variable))+
  geom_line()+ggtitle("Proportion of Secondary Males by Number of Infants")+
  labs(color='Offspring')+xlab("Proportion of Secondary Male Breeders");g2

#repeat with pFs
pFs <- data.frame(pFs=seq(0,.6,length.out=64), grp.size = mean(comp$grp.size),pFp=mean(comp$pFp))
predict(mn4, newdata = pFs, "probs")
pp.pFs<-cbind(pFs,predict(mn4, newdata = pFs, "probs",se=T));head(pp.pFs)
pp.pFs<-subset(pp.pFs,select=);head(pp.pFs)
PpFs <- melt(pp.pFs, id.vars = c("pFs","grp.size","pFp"),
             value.name = "probability");names(PpFs)[4]<-"Offspring"
head(PpFs)
g3<-ggplot(PpFs, aes(x = pFs, y = probability, colour = Offspring))+
  geom_line()+ggtitle("Proportion of Secondary Females by Number of Infants")+
  labs(color='Offspring')+xlab("Proportion of Female Secondary Breeders");g3


Pgraphs<-rbind.fill(melt(PpFp,id.vars=c("grp.size","Offspring","probability"),value.name="prop"),
                    melt(Pgrp,id.vars=c("pFp","Offspring","probability"),value.name="prop"))
head(Pgraphs)
levels(Pgraphs$variable)<-c("Prop. Primary Female","Group Size")

ggplot(Pgraphs, aes(x = prop, y = probability, color = Offspring))+
  geom_line(lwd = 1) +
  theme(axis.title.x=element_blank(), panel.grid.major=element_blank(),
        text=element_text(size=20))+
  facet_grid(.~variable,scale="free")+
  ggtitle("Results From Multinomial Logistic Regression")
















scale_linetype_manual(values = c(1,5,3))+








theme_set(theme_gray(base_size = 20))+
theme(axis.text.x = element_blank(),panel.grid.major = element_blank())+
  

# now looking at what explain group size
library(GGally)
ggpairs(comp[,c("grp.size","M","F","F_Primary","F_Secondary","Juvenile")])

ggpairs(comp[,c("Juvenile","M","F","M_Primary","F_Primary","M_Secondary",
                "F_Secondary","pMp","pFp")])

library(car)
scatterplotMatrix(~grp.size+M+F+F_Primary+F_Secondary+Juvenile|Species,data=comp,
                  main="Scatterplot by species")

sGroup<-glm(formula=grp.size~M+F+F_Primary+F_Secondary+Juvenile,
               data=comp,family=poisson);summary(sGroup)

sGroup1<-glm(formula=grp.size~M+F,
            data=comp,family=poisson);summary(sGroup1)


  
  






melt(Pgrp,id.vars=c("pFp","pMp","variable","probability"),value.name="prop.")

melt(Pgraphs,id.vars=c("pFp","pMp","variable","probability","type"),value.name="prop.")

                                                                       

rep(c("a","b","c"),each=67)


# Model 5 is the best fit group, proportion of primary females and proportion of breeders





predict(house.plr, housing, type = "p")
addterm(house.plr, ~.^2, test = "Chisq")
house.plr2 <- stepAIC(house.plr, ~.^2)
house.plr2$anova
anova(house.plr, house.plr2)

exp(coef(mn5))

head(pp <- fitted(test))


par(mfrow = c(2,2), mar = c(5,5,1.5,0.5))
plot(test, which = c(1:4), add.smooth = F, cex.id = 1)




############## Ordered logistic regression #################
# attempted to do an ordinal logistic regression because the number of juveniles is ordinal
# but our data seems to violate the assumption of proportional odds in a big way
#http://www.ats.ucla.edu/stat/r/dae/ologit.htm

require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
require(matrixStats)

comp$Juv<-as.factor(comp$Juvenile)

m<-polr(Juv~grp.size+pFp+pMp+pFs+pMp,data=na.omit(comp),Hess=T)
m
summary(m, digits=3)



# to estimate p-values
(ctable <- coef(summary(m)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))

# confidence intervals
(ci <- confint(m)) # default method gives profiled CIs

# CIs assuming normality
confint.default(m) # CIs assuming normality

# Odds ratio
exp(coef(m))

# OR and CI
exp(cbind(OR = coef(m), ci))


sf <- function(y) {
  c('Y>=1' = qlogis(mean(y >= 1)),
    'Y>=2' = qlogis(mean(y >= 2)),
    'Y>=3' = qlogis(mean(y >= 3)),
    'Y>=0' = qlogis(mean(y >= 0)))
  
}

(s <- with(na.omit(comp),
           summary(as.numeric(Juv) ~ grp.size+pFp+pMp+pFs+pMp, fun=sf)))

glm(I(as.numeric(Juv) >= 2) ~ grp.size, family="binomial", data = na.omit(comp))

glm(I(as.numeric(Juv) >= 3) ~ grp.size, family="binomial", data = na.omit(comp))

s[, 4] <- s[, 4] - s[, 3]
s[, 3] <- s[, 3] - s[, 3]
s # print

plot(s, which=1:4,pch=1:3, xlab='logit', main=' ', xlim=c(-3.0,0.0))



plot(fitted[1],comp$Juvenile)

predict()

m$fitted.values
m$df.residual


options(contrasts = c("contr.treatment", "contr.poly"))
house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
house.plr
summary(house.plr, digits = 3)
## slightly worse fit from
summary(update(house.plr, method = "probit", Hess = TRUE), digits = 3)
## although it is not really appropriate, can fit
summary(update(house.plr, method = "loglog", Hess = TRUE), digits = 3)
summary(update(house.plr, method = "cloglog", Hess = TRUE), digits = 3)

predict(house.plr, housing, type = "p")
addterm(house.plr, ~.^2, test = "Chisq")
house.plr2 <- stepAIC(house.plr, ~.^2)
house.plr2$anova
anova(house.plr, house.plr2)

house.plr <- update(house.plr, Hess=TRUE)
pr <- profile(house.plr)
confint(pr)
plot(pr)
pairs(pr)


cor(comp[c(7,5,25,26,28)])


comp3<-subset(comp2, comp2$grp.size>4);comp3<-droplevels(comp3)
cor(comp3[c(7,5,13,15,16)])


################## Making model trees

library(rpart)

test <- rpart(as.factor(Juvenile)~grp.size+pFp, data = comp, method = "class",
              minsplit=5, maxdepth=4, cp = 0.001)
plot(test,uniform=T,margin=0.1)
text(test, use.n = T, cex = 1.0)

test1 <- rpart(Juv~grp.size+pFp+pMp+pFs+pMp, data = comp, method = "class",
               minsplit=5, cp = 0.001)

plot(test1,uniform=T,margin=0.1)
text(test1, use.n = T, cex = 1.0)


###############
############# Generalized Additive Modeling ##############################

library(mgcv)
gitHi<-gam(Juvenile~grp.size+pFp+s(pMp)+pFs+pMp,
           data=comp,family=poisson(link=log))
summary(gitHi)

gitHi2<-gam(Juvenile~grp.size+pFp+pFs+pMp,
            data=comp,family=poisson(link=log))
summary(gitHi2)


##################################
#repeate with pMp
pMp <- data.frame(pMp=seq(0,.857,length.out=67), grp.size = mean(comp$grp.size),pFp=mean(comp$pFp))
predict(mn3, newdata = pMp, "probs")
pp.pMp<-cbind(pMp,predict(mn3, newdata = pMp, "probs",se=T));head(pp.pMp)
pp.pMp<-subset(pp.pMp,select=-7);head(pp.pMp)
PpMp <- melt(pp.pMp, id.vars = c("pMp","grp.size","pFp"),
             value.name = "probability");names(PpMp)[4]<-"Offspring"
head(PpMp)
g3<-ggplot(PpMp, aes(x = pMp, y = probability, colour = variable))+
  geom_line()+ggtitle("Proportion of Primary Male Breeders by Number of Infants")+
  labs(color='Infants')+xlab("Proportion of Male Primary Breeders")


# getting coefficients and p-values from models

library(AER)
coeftest(mn5)
coeftest(mn4)

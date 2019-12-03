# Primate Blood Parasite Dynamics
PCR and blood smear infection data

|Ind. Changes Infection Status|Infections by Age Class|
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
### Making the Sample Figures
```rscript

parasite<-read.xlsx("parasite.xlsx",1)

library(grid)

parasite$M.mariae<-ifelse(parasite$M.mariae==1,"Pos","Neg")
parasite$D.gra<-ifelse(parasite$D.gra==1,"Pos","Neg")
parasite$T.min<-ifelse(parasite$T.min==1,"Pos","Neg")

names(parasite)
  pp<-subset(parasite,select=c(2:7,24:26,27,29:30,32))
head(pp)

ppp<-melt(pp, id=c("Species","Sex","Group","AnID","Name","Age_class","year","nWeight","br.dfa"))
head(ppp)
str(ppp)

# fix orders of juveniles
ppp$Age_class = factor(ppp$Age_class,levels(ppp$Age_class)[c(1,3,2)])

#re-name species  and parasites
levels(ppp$variable)<-c("M. mariae", "Dipetalonema spp.","T. minasense","PSR")
levels(ppp$Species)<-c("S. weddelli", "S. imperator")

####################################################################
# run this code before plotting change in infection status over time
library(proto)
# Detect and prevent collisions.
# Powers dodging, stacking and filling.
collidev <- function(data, height = NULL, name, strategy, check.height = TRUE) {
  # Determine height
  if (!is.null(height)) {
    # height set manually
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y - height / 2
      data$ymax <- data$y + height / 2
    }
  } else {
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y
      data$ymax <- data$y
    }
    
    # height determined from data, must be floating point constant
    heights <- unique(data$ymax - data$ymin)
    heights <- heights[!is.na(heights)]
    
    #   # Suppress warning message since it's not reliable
    #     if (!zero_range(range(heights))) {
    #       warning(name, " requires constant height: output may be incorrect",
    #         call. = FALSE)
    #     }
    height <- heights[1]
  }
  
  # Reorder by x position, relying on stable sort to preserve existing
  # ordering, which may be by group or order.
  data <- data[order(data$ymin), ]
  
  # Check for overlap
  intervals <- as.numeric(t(unique(data[c("ymin", "ymax")])))
  intervals <- intervals[!is.na(intervals)]
  
  if (length(unique(intervals)) > 1 & any(diff(scale(intervals)) < -1e-6)) {
    warning(name, " requires non-overlapping y intervals", call. = FALSE)
    # This is where the algorithm from [L. Wilkinson. Dot plots.
    # The American Statistician, 1999.] should be used
  }
  
  if (!is.null(data$xmax)) {
    plyr::ddply(data, "ymin", strategy, height = height)
  } else if (!is.null(data$x)) {
    data$xmax <- data$x
    data <- plyr::ddply(data, "ymin", strategy, height = height)
    data$x <- data$xmax
    data
  } else {
    stop("Neither x nor xmax defined")
  }
}

# Stack overlapping intervals.
# Assumes that each set has the same horizontal position
pos_stackv <- function(df, height) {
  if (nrow(df) == 1) return(df)
  
  n <- nrow(df) + 1
  x <- ifelse(is.na(df$x), 0, df$x)
  if (all(is.na(df$y))) {
    heights <- rep(NA, n)
  } else {
    heights <- c(0, cumsum(x))
  }
  
  df$xmin <- heights[-n]
  df$xmax <- heights[-1]
  df$x <- df$xmax
  df
}

# Stack overlapping intervals and set height to 1.
# Assumes that each set has the same horizontal position.
pos_fillv <- function(df, height) {
  stacked <- pos_stackv(df, height)
  stacked$xmin <- stacked$xmin / max(stacked$xmax)
  stacked$xmax <- stacked$xmax / max(stacked$xmax)
  stacked$x <- stacked$xmax
  stacked
}

# Dodge overlapping interval.
# Assumes that each set has the same horizontal position.
pos_dodgev <- function(df, height) {
  n <- length(unique(df$group))
  if (n == 1) return(df)
  
  if (!all(c("ymin", "ymax") %in% names(df))) {
    df$ymin <- df$y
    df$ymax <- df$y
  }
  
  d_height <- max(df$ymax - df$ymin)
  
  # df <- data.frame(n = c(2:5, 10, 26), div = c(4, 3, 2.666666,  2.5, 2.2, 2.1))
  # ggplot(df, aes(n, div)) + geom_point()
  
  # Have a new group index from 1 to number of groups.
  # This might be needed if the group numbers in this set don't include all of 1:n
  groupidy <- match(df$group, sort(unique(df$group)))
  
  # Find the center for each group, then use that to calculate xmin and xmax
  df$y <- df$y + height * ((groupidy - 0.5) / n - .5)
  df$ymin <- df$y - d_height / n / 2
  df$ymax <- df$y + d_height / n / 2
  
  df
}


#' Adjust position by dodging overlaps to the side.
#'
#' @inheritParams ggplot2::position_identity
#' @param height Dodging height, when different to the height of the individual
#'   elements. This is useful when you want to align narrow geoms with wider
#'   geoms. See the examples for a use case.
#' @family position adjustments
#' @export
#' @examples
#' ggplot(mtcars, aes(factor(cyl), fill = factor(vs))) +
#'   geom_bar(position = "dodge")
#' \donttest{
#' ggplot(diamonds, aes(price, fill = cut)) +
#'   geom_histogram(position="dodge")
#' # see ?geom_boxplot and ?geom_bar for more examples
#'
#' # To dodge items with different heights, you need to be explicit
#' df <- data.frame(x=c("a","a","b","b"), y=2:5, g = rep(1:2, 2))
#' p <- ggplot(df, aes(x, y, group = g)) +
#'   geom_bar(
#'     stat = "identity", position = "dodge",
#'     fill = "grey50", colour = "black"
#'   )
#' p
#'
#' # A line range has no height:
#' p + geom_linerange(aes(ymin = y-1, ymax = y+1), position = "dodge")
#' # You need to explicitly specify the height for dodging
#' p + geom_linerange(aes(ymin = y-1, ymax = y+1),
#'   position = position_dodge(width = 0.9))
#'
#' # Similarly with error bars:
#' p + geom_errorbar(aes(ymin = y-1, ymax = y+1), width = 0.2,
#'   position = "dodge")
#' p + geom_errorbar(aes(ymin = y-1, ymax = y+1, height = 0.2),
#'   position = position_dodge(width = 0.90))
#' }
position_dodgev <- function(height = NULL) {
  ggproto(NULL, PositionDodgeV, height = height)
}



PositionDodgeV <- ggproto("PositionDodgeV", Position,
                          required_aes = "y",
                          height = NULL,
                          setup_params = function(self, data) {
                            if (is.null(data$ymin) && is.null(data$ymax) && is.null(self$height)) {
                              warning("height not defined. Set with `position_dodgev(height = ?)`",
                                      call. = FALSE)
                            }
                            list(height = self$height)
                          },
                          
                          compute_panel = function(data, params, scales) {
                            collidev(data, params$height, "position_dodgev", pos_dodgev, check.height = FALSE)
                          }
)
##################################################################

# now plot the graph of individuals changes in infection status
ggplot(ppp[ppp$variable!="PSR",],aes(x=year,y=value,group=AnID, color=variable))+
  geom_point(position=position_dodgev(height=.2),shape=3)+
  scale_x_discrete(expand=c(0,.2))+scale_y_discrete(expand=c(.04,0))+
  geom_line(position=position_dodge(width=.1),alpha=.2,size=2)+
  geom_line(position=position_dodgev(height=.1),alpha=.2,size=2)+
  facet_grid(~variable, scales="free")+
  theme_set(theme_bw(base_size = 16))+
  labs(x="Year",y="Infection status",
       title="Individual change in infection status over time")

#plot parasite species richness
ppp$Age_class[ppp$Age_class=="Sub-adult"]<-"Juvenile"
ppp<-droplevels(ppp)
levels(ppp$Age_class)<-c("Adult", "Juveniles & Sub-adults")

cbPalette <- c("#999999","#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(ppp[ppp$variable=="PSR",],aes(x=value,fill=Sex))+
  geom_bar()+facet_grid(Species~Age_class)+
  scale_fill_manual(values=c("black","grey"))+
  theme_set(theme_bw(base_size = 14))+
  labs(x="PSR",y="n",
       title="Blood parasite species richness by age class")

```
### Modeling Blood Parasite Infection Across Demographic Variables
```rscript

# Modeling Dipetalonema gracile infection presence/absence
setwd("~/Documents/Work/Research/Data/Analysis/Haemoparasites")

# load packages we need
library(lme4) #one method for doing a GLMM use the glmer function
library(MASS) # another GLMM method useing the glmmPQL function
library(glmmML) # another GLMM method using the glmmML function
library(lmtest) # liklihood ratio test can be used to compare models or not?
library(car)
library(xlsx)

# load the parasite data and create a working version
parasite<-read.xlsx("parasite.xlsx",1); parasite4<-parasite

# load prior infection dataset
para.priors<-read.xlsx("para.priors.xlsx",1)

# load group composition information
grp.comp<-read.xlsx("grp.comp.xlsx",1)

# remove 2009,2010,2011, and 2015 from group comp data
comp1<-subset(grp.comp,grp.comp$year!=2009 & grp.comp$year!=2010 & grp.comp$year!=2011 &
                grp.comp$year!=2015); comp1<-droplevels(comp1)

# sebset the group comp data set to just have group names, years, and group size
comp2<-subset(comp1,select=c("Group","year","grp.size"))

# create a new interaction column for each data set based on year and group name
parasite4$inter<-interaction(parasite4$Group,parasite4$year)
para.priors$inter<-interaction(para.priors$Group,para.priors$year)
comp2$inter<-interaction(comp2$Group,comp2$year)

# using the match function, insert group size into the parasite and para.priors
# datasets; and we will now have group size associated with each parasite infection
# record
parasite4$grp.size<-comp2[match(parasite4$inter,comp2$inter),3]
para.priors$grp.size<-comp2[match(para.priors$inter,comp2$inter),3]

# scale the weight variable in case there are problems modeling
para.priors$Weight<-scale(para.priors$nWeight,scale=F)

# Checking for optimal random structure

# first model has no random structure and reveals that year
# should be a random effect because it significantly effects model
DR1<-glm(D.gra~Species+Sex+Age_class+br.dfa+year+scale(nWeight,scale=F)+
         Group,family=binomial,data=parasite4);DR1
summary(DR1)

# year as random effect, fit is improved
DR2<-glmer(D.gra~Species+Sex+Age_class+br.dfa+scale(nWeight,scale=F)+
                (1|year),
                family=binomial,data=parasite4);DR2
summary(DR2) 

# group as random effect
DR3<-glmer(D.gra~Species+Sex+Age_class+br.dfa+scale(nWeight,scale=F)+
                (1|Group),
                family=binomial,data=parasite4);DR3
summary(DR3) 

# group and year as randome effects
DR4<-glmer(D.gra~Species+Sex+Age_class+br.dfa+scale(nWeight,scale=F)+
                (1|year)+(1|Group),
                family=binomial,data=parasite4);DR4
summary(DR4) 

AIC(DR1,DR2,DR3,DR4)

lrtest(DR3,DR4)

# Including random structure preserved important degrees of freedoom, also it was 
# necessary because year has a siginificant effect when used as fixed factor
# Group did not have a significant effect as a fixed factor, nor did it improve
# model performance when incorporated into the random structure. lrtest confirms
# including group, in addition to year, as a random intercept is not significantly different.


# Dipetalonema using glmer ###############################

#this first model needs to have the weight variable re-scaled
D1.glmer<-glmer(D.gra~Species+Age_class+Sex+br.dfa+scale(nWeight)+(1|year)+(1|AnID),
                family=binomial,data=parasite4);D1.glmer
summary(D1.glmer) #remove Sex

D2.glmer<-glmer(D.gra~Species+Age_class+br.dfa+scale(nWeight)+(1|year)+(1|AnID),
                family=binomial,data=parasite4);D2.glmer
summary(D2.glmer) #remove Breeding status

D3.glmer<-glmer(D.gra~Species+Age_class+scale(nWeight)+(1|year)+(1|AnID),
                family=binomial,data=parasite4);
summary(D3.glmer) #remove Age_class

D4.glmer<-glmer(D.gra~Species+scale(nWeight)+(1|year)+(1|AnID),
                family=binomial,data=parasite4);D4.glmer
summary(D4.glmer) 

D5.glmer<-glmer(D.gra~scale(nWeight)+(1|year)+(1|AnID),
                family=binomial,data=parasite4);D5.glmer
summary(D5.glmer) 

lrtest(D5.glmer,D4.glmer)

Anova(D5.glmer)

cbind(BIC(D1.glmer,D2.glmer,D3.glmer,D4.glmer),
      AIC(D1.glmer,D2.glmer,D3.glmer,D4.glmer))


################Dipetalonema using glmmPQL reveals same result

D1.pql<-glmmPQL(D.gra~Species+Sex+Age_class+br.dfa+nWeight,random=~1|year,
                family=binomial, data=parasite4);D1.pql
summary(D1.pql)#remove Sex

D2.pql<-glmmPQL(D.gra~Species+Age_class+br.dfa+nWeight,random=~1|year,
                family=binomial, data=parasite4);D2.pql
summary(D2.pql)#remove br.dfa

D3.pql<-glmmPQL(D.gra~Species+Age_class+nWeight,random=~1|year,
                family=binomial, data=parasite4);D3.pql
summary(D3.pql)# remove Age_class

D4.pql<-glmmPQL(D.gra~Species+nWeight,random=~1|year,
                family=binomial, data=parasite4);D4.pql
summary(D4.pql)# best model


################################### Adding co-infection
# Begin by add co-infection variabels to the best fit model
summary(D5.glmer) # this was the best fit

D6.glmer<-glmer(D.gra~scale(nWeight)+M.mariae+T.min+(1|year)+(1|AnID),
                family=binomial,data=parasite4);
summary(D6.glmer) #remove T.min 

D7.glmer<-glmer(D.gra~scale(nWeight)+M.mariae+(1|year)+(1|AnID),
                family=binomial,data=parasite4);
summary(D7.glmer) # M.mariae

D8.glmer<-glmer(D.gra~scale(nWeight)+(1|year)+(1|AnID),
                family=binomial,data=parasite4);
summary(D8.glmer)

lrtest(D7.glmer,D5.glmer)

Anova(D7.glmer)



cbind(AIC(D4.glmer,D5.glmer,D6.glmer,D7.glmer),
      BIC(D4.glmer,D5.glmer,D6.glmer,D7.glmer))


############################## Adding infection history

# like before we start with the current best fit model, including co-infection
summary(D7.glmer) # notice n=186, we need to make sure that our reduced n
# for testing influences of prior infections is still picking up the same
# significant factors
D8.glmer<-glmer(D.gra~scale(nWeight,scale=F)+M.mariae+(1|year),
                family=binomial,data=para.priors);
summary(D8.glmer) # n = 74 but same results

# add all prior infection variables
D9.glmer<-glmer(D.gra~scale(nWeight,scale=F)+M.mariae+
                  D.gra_prior+M.mar_prior+T.min_prior+(1|year),
                family=binomial,data=para.priors);
summary(D9.glmer) # remove M.mariae_prior infection

D10.glmer<-glmer(D.gra~Weight+M.mariae+
                  D.gra_prior+T.min_prior+(1|year),
                family=binomial,data=para.priors);
summary(D10.glmer) # Appears to be best model, removal of more terms
# does not significantly imrpove model fit according to AIC, though it
# does according to BIC.  BIC penalizes more for including marginal terms.
# Likelihood ratio tests suggest the inclusion of these terms is appraoching
# significance (p-value = .07)

D11.glmer<-glmer(D.gra~scale(nWeight,scale=F)+M.mariae+
                   D.gra_prior+(1|year), family=binomial, data=para.priors)
summary(D11.glmer) 

D12.glmer<-glmer(D.gra~M.mariae+
                   D.gra_prior+(1|year),
                 family=binomial,data=para.priors);
summary(D12.glmer) 

D13.glmer<-glmer(D.gra~D.gra_prior+(1|year),
                 family=binomial,data=para.priors);
summary(D13.glmer) 

cbind(AIC(D8.glmer,D9.glmer,D10.glmer,D11.glmer,D12.glmer,D13.glmer),
      BIC(D8.glmer,D9.glmer,D10.glmer,D11.glmer,D12.glmer,D13.glmer))

lrtest(D10.glmer,D13.glmer) 

```
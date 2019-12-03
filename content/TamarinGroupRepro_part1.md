# set workind dir
setwd("~/Documents/Work/Research/Data/Analysis/BrStatus&GroupComp")

#load excel packages
require(xlsx)
require(plyr)
require(reshape2)
require(ggplot2)

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

status<-read.xlsx("BreedingStatus3.9.16.xlsx",1) #get the most current br.status file

#subset status to only the columns of interest
data<-subset(status, select=c("AnID","Date","Species","Sex","Group","Name","Age_class",
                              "year","br.dfa"))

# change all disperser categories to secondary
data$F.young<-ifelse(data$Sex=="F" & data$Age_class!="Juvenile" &
                       data$br.dfa=="None","Fy","NA")
data$M.young<-ifelse(data$Sex=="M" & data$Age_class!="Juvenile" &
                       data$br.dfa=="None","My","NA")

data1<-dcast(data, Species+Group+year~Age_class)
data2<-dcast(data, Species+Group+year~Sex+br.dfa)
data3<-dcast(data, Species+Group+year~br.dfa)
data4<-ddply(data,.(Species,Group,year), nrow)
data5<-dcast(data, Species+Group+year~Sex)
data6<-dcast(data, Species+Group+year~M.young)
data7<-dcast(data, Species+Group+year~F.young)
data8<-dcast(data, Species+Group+year~Age_class+Sex)


rm(comp2)
comp2<-merge(data4,data1, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data5, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data2, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data3, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data6, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data7, by=c("Species","Group","year"), all=T)
comp2<-merge(comp2,data8, by=c("Species","Group","year"), all=T)

colnames(comp2)[4]<-"grp.size"
names(comp2)
comp2a<-subset(comp2,select=-c(19,21))



comp2a<-ddply(comp2a,.(Species),mutate,rSex=round(M/F,2),grp.size.new=grp.size-Juvenile,
              Breeders=Secondary+Primary,pBreeder=round(Breeders/grp.size.new,2),
              F.new=F_Primary+F_Secondary+Fy, M.new=M_Primary+M_Secondary)

#remove the callicebus
comp2a<-subset(comp2a, comp2a$Species!="CBRU")


## need to write the data frame to a file and make corrections, then can reimport and model
write.xlsx(comp2a,"grp.comp.xlsx")


# before importing grp.comp.xlsx do the following
# confirm that all 2009 records are removed

#correct FC row 2010 to have 3 (GRC,GPG,GBR) adults, 2 P_females,
# 1 P_Male, 2 Female juvs (GPO and GPBL) the numbers are off because
# we retrapped the juvs many times in 2010

# also, add one primary breeding female to AR6 2013


comp<-read.xlsx("grp.comp.xlsx",1)

comp$grp.size<-comp$grp.size.new
comp$F<-comp$F.new
comp$M<-comp$M.new
names(comp)
comp<-subset(comp, select=-c(28,31,32))

# we need to get proportions because group size is correlated with counts of everything
comp<-ddply(comp,.(Species), mutate,pF=F/grp.size,pM=M/grp.size,
            pFp=F_Primary/grp.size,pMp=M_Primary/grp.size,
            pFs=F_Secondary/grp.size,pMs=M_Secondary/grp.size,rSex=M/F,
            helpers=M_Primary+M_Secondary+F_Secondary-1,pHelp=helpers/grp.size,
            pAM=Adult_M/grp.size,pAF=Adult_F/grp.size,pSM=Sub.adult_M/grp.size,
            pSF=Sub.adult_F/grp.size)

# remove loners that don't have a group
comp<-subset(comp,comp$Group!="Loners")

# remove group MI4 2012 because we didn't actually trap this group, we caught 2
#     individuals that were attempting to disperse
# Remove angels, not really a group
# Remove Orange Emps 2012, not really a group
# remove bees from 2013 becuase it was a newly formed group
comp<-subset(comp, comp$Group!="Angels" & comp$Group!="");comp<-droplevels(comp)
comp<-subset(comp,comp$Group!="OrangeEmps" | comp$year!="2012");comp<-droplevels(comp)
comp<-subset(comp,comp$Group!="MI4" | comp$year!="2012");comp<-droplevels(comp)
comp<-comp[comp$Group!="Bees"|comp$year!=2013,]

#data frame is ready for modeling
write.xlsx(comp, "comp.xlsx")

# then decide which groups are establiedh (old) and unestablished (new)

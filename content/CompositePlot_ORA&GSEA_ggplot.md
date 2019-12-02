|Example|
|---|
|![Example](https://gideonerkenswick.files.wordpress.com/2019/12/ora-gseaplot_example-1.jpg)|
```rscript
# Script Name: CompositePlot_ORA&GSEA_ggplot

## Overview ##

#1) Customize bar graphs
#2) Combine multiple bar graphs into a composite plot

## Require packages ##

library(ggplot2) # data plotting package
library(stringr) # character string manipulation
library(jcolors) # alternative color scheme for ggplot
library(ggpubr) # combine ggplot objects

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# prepare GSEA plot1

#import data
rvBP<-read.table("path/to/enrichment/results.txt",header = T,sep="\t")
rvMF<-read.table("path/to/enrichment/results.txt",header = T,sep="\t")
rvCC<-read.table("path/to/enrichment/results.txt",header = T,sep="\t")

# create GO database variable
rvBP$db<-"BP";rvMF$db<-"MF";rvCC$db<-"CC"

rv<-rbind(rvBP,rvMF,rvCC) # combine results

# shorten GO term description names (add more lines as needed)
rv$description<-gsub('response','resp.',rv$description)
rv$description<-gsub('positive','pos.',rv$description)
rv$description<-gsub('oxygen','oxyg.',rv$description)
rv$description<-gsub('production','produc.',rv$description)
rv$description<-gsub('process','proce.',rv$description)
rv$description<-gsub('reactive','react.',rv$description)
rv$description<-gsub('regulation','reg.',rv$description)
rv$description<-gsub('external','ext.',rv$description)
rv$description<-gsub('metabolic','metabo.',rv$description)

rv<-arrange(rv,normalizedEnrichmentScore) #arrange data be enrichment score
newlevel<-rv$description #extract ordered descriptions # extract ordered descriptions
rv$description<-factor(rv$description,levels=unique(newlevel)) # replace description factor levels with ordered descriptions
rv$reg<-ifelse(rv$normalizedEnrichmentScore>0,"pos","neg") # create directional variable for term enrichment

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot GSEA (plot1)

# provide manual fill colors, so that db color (BP,CC,MF) will be standardized across all subsequent plots
# remove x-axis labels and instead use geom_text to place description label opposite of each bar (on the other side of the x-axis)
# color geom_text to be distinctive for each GO DB (BP,CC,MF)

p1<-ggplot(rv,aes(description,normalizedEnrichmentScore,fill=reg))+
  geom_bar(stat="identity")+coord_flip()+
  scale_fill_manual(values = c("orange3","lightblue3"))+
  geom_text(data=subset(rv,rv$normalizedEnrichmentScore>0 & rv$db=="BP"),
            aes(label=description),y=-.03,hjust=1,size=3,family="serif")+
  geom_text(data=subset(rv,rv$normalizedEnrichmentScore>0 & rv$db=="MF"),
            aes(label=description),y=-.03,hjust=1,size=3,family="serif",color="darkgreen")+
  geom_text(data=subset(rv,rv$normalizedEnrichmentScore>0 & rv$db=="CC"),
            aes(label=description),y=-.03,hjust=1,size=3,family="serif",color="blue4")+
  geom_text(data=subset(rv,rv$normalizedEnrichmentScore<0 & rv$db=="BP"),
            aes(label=description),y=.03,hjust=0,size=3,family="serif")+
  geom_text(data=subset(rv,rv$normalizedEnrichmentScore<0 & rv$db=="MF"),
            aes(label=description),y=.03,hjust=0,size=3,family="serif",color="green4")+
  geom_text(data=subset(rv,rv$normalizedEnrichmentScore<0 & rv$db=="CC"),
            aes(label=description),y=.03,hjust=0,size=3,family="serif",color="blue4")+
  theme(legend.position = "none",axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "grey"),
        text = element_text(size=12))+
  labs(y="Normalized Enrichment",x="",title="");p1


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# prepare ORA (plot2)

#import enrichment results
rvBPoraup<-read.table("path/to/enrichment/results.txt",header = T,sep="\t")
rvMForaup<-read.table("path/to/enrichment/results.txt",header = T,sep="\t")
rvCCoraup<-read.table("path/to/enrichment/results.txt",header = T,sep="\t")
rvCCoradwn<-read.table("path/to/enrichment/results.txt",header = T,sep="\t")

# add GP 'db' variable (BP,CC,MF)
# add negative and positive gene regulation variable 'dir'
rvBPoraup$db<-"BP";rvMForaup$db<-"MF";rvCCoraup$db<-"CC";rvBPoraup$dir<-"up";rvMForaup$dir<-"up";rvCCoraup$dir<-"up"
rvCCoradwn$dir<-"down";rvCCoradwn$db<-"CC";rvCCoradwn$enrichmentRatio<-rvCCoradwn$enrichmentRatio*-1

rvora<-rbind(head(rvBPoraup,20),rvMForaup,rvCCoraup,rvCCoradwn) # combine results

# shorten GO term description names (add more lines as needed)
rvora$description<-gsub('response','resp.',rvora$description)
rvora$description<-gsub('stimulus','stim.',rvora$description)
rvora$description<-gsub('microtubule','microtub.',rvora$description)
rvora$description<-gsub('part','',rvora$description)
rvora$description<-gsub('signaling','signl.',rvora$description)
rvora$description<-gsub('of ','',rvora$description)
rvora$description<-gsub('to ','',rvora$description)
rvora$description<-gsub('molecule','mol.',rvora$description)
rvora$description<-gsub('pathway','path',rvora$description)
rvora$description<-gsub('bacterial','bact.',rvora$description)
rvora$description<-gsub('nitrogen','N',rvora$description)
rvora$description<-gsub('cytokine','cyto.',rvora$description)
rvora$description<-gsub('cellular','cella.',rvora$description)
rvora$description<-gsub('positive','pos.',rvora$description)
rvora$description<-gsub('oxygen','O',rvora$description)


#arrange data by enrichment ratio
rvora<-arrange(rvora,enrichmentRatio)
newlevel<-rvora$description # extract ordered GP descriptions
rvora$description<-factor(rvora$description,levels=unique(newlevel)) # replace description factor levels with ordered descriptions

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot ORA (plot2)

# provide manual fill colors, so that db color (BP,CC,MF) will be standardized across all subsequent plots
# remove x-axis labels and instead use geom_text to place description label opposite of each bar (on the other side of the x-axis)
# color geom_text to be distinctive for each GO DB (BP,CC,MF)

p2<-ggplot(rvora,aes(description,enrichmentRatio,fill=dir))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("orange3","lightblue3"))+
  geom_text(data=subset(rvora,rvora$enrichmentRatio>0 & rvora$db=="BP"),
            aes(label=description),y=-.1,angle=90,hjust=1,size=3,family="serif")+
  geom_text(data=subset(rvora,rvora$enrichmentRatio>0 & rvora$db=="MF"),
            aes(label=description),y=-.1,angle=90,hjust=1,size=3,family="serif",color="darkgreen")+
  geom_text(data=subset(rvora,rvora$enrichmentRatio>0 & rvora$db=="CC"),
            aes(label=description),y=-.1,angle=90,hjust=1,size=3,family="serif",color="blue4")+
  geom_text(data=subset(rvora,rvora$enrichmentRatio<0 & rvora$db=="BP"),
            aes(label=description),y=.1,angle=90,hjust=0,size=3,family="serif")+
  geom_text(data=subset(rvora,rvora$enrichmentRatio<0 & rvora$db=="MF"),
            aes(label=description),y=.1,angle=90,hjust=0,size=3,family="serif",color="green4")+
  geom_text(data=subset(rvora,rvora$enrichmentRatio<0 & rvora$db=="CC"),
            aes(label=description),y=.1,angle=90,hjust=0,size=3,family="serif",color="blue4")+
  theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.y = element_line(colour = "grey"),
        text = element_text(size=12))+
  labs(y="Enrichment Ratio",x="",title="");p2


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# dH  ORA (plot3): REPEAT PREVIOUS STEPS FOR ENRICHMENT RESULTS

#import enrichment results
dhBP<-read.table("path/to/enrichment/results.txt",header = T,sep="\t")
dhMF<-read.table("path/to/enrichment/results.txt",header = T,sep="\t")

#add db variable; pos/neg gene expression variable not necessary (all are upregulated genes)
dhBP$db<-"BP";dhMF$db<-"MF"
dh<-rbind(dhBP,dhMF) # combind the enrichment results

# SAME AS ABOVE
dh$description<-gsub('response','resp.',dh$description)
dh$description<-gsub('positive','pos.',dh$description)
dh$description<-gsub('oxygen','oxyg.',dh$description)
dh$description<-gsub('production','produc.',dh$description)
dh$description<-gsub('process','proce.',dh$description)
dh$description<-gsub('reactive','react.',dh$description)
dh$description<-gsub('regulation','reg.',dh$description)
dh$description<-gsub('external','ext.',dh$description)
dh$description<-gsub('metabolic','metabo.',dh$description)

# SAME AS ABOVE
dh<-arrange(dh,enrichmentRatio)
newlevel<-dh$description
dh$description<-factor(dh$description,levels=unique(newlevel))
dh$reg<-ifelse(dh$enrichmentRatio>0,"pos","neg")

# SAME AS ABOVE
p3<-ggplot(dh,aes(description,enrichmentRatio,fill=reg))+
  geom_bar(stat="identity")+coord_flip()+
  scale_fill_manual(values = c("orange3","lightblue3"))+
  geom_text(data=subset(dh,dh$enrichmentRatio>0 & dh$db=="BP"),
            aes(label=description),y=-1,hjust=1,size=3,family="serif")+
  geom_text(data=subset(dh,dh$enrichmentRatio>0 & dh$db=="MF"),
            aes(label=description),y=-1,hjust=1,size=3,family="serif",color="darkgreen")+
  geom_text(data=subset(dh,dh$enrichmentRatio<0),
            aes(label=description),y=.1,hjust=0,size=3,family="serif")+
  theme(legend.position = "none",axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "grey"),
        text = element_text(size=12))+ylim(-28.5,30)+
  labs(y="Enrichment Ratio",x="",title="");p3
  
  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# write composite plot

g1<-ggarrange(p1,p3,nrow=1,ncol=2,heights = c(3,1.15),labels = list("B","C"))

ggarrange(g1,p2,nrow=2,ncol=1,heights = c(2,1.5),labels = list(" ","A"))

# save composite plot to .svg (fine tune graphic in illustrator)
ggsave("figure2.svg",width=8,height = 10)

```



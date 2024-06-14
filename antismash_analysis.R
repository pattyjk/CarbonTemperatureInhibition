library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)

#read in antismash data
dat<-read.delim('CarbonTemperatureInhibition/isolate_antismash.txt', header=T)

#summarize
anti_sum<-ddply(dat, c("Isolate", "Genus", 'Type'), summarize, no_genes=length(Type))

#plot it
ggplot(anti_sum, aes(Isolate, no_genes, fill=Type))+
  geom_bar(stat='identity')+
  facet_wrap(~Genus, scales='free_x')+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  ylab("Number of BSG")+
  xlab("")+
  theme_bw()


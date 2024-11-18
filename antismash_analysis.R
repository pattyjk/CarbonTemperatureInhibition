library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)

#read in antismash data
dat<-read.delim('Documents/GitHub/CarbonTemperatureInhibition/isolate_antismash.txt', header=T)

#summarize
anti_sum<-ddply(dat, c("Isolate", "Genus", 'Type'), summarize, no_genes=length(Type))

#make best pal
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

#plot it
ggplot(anti_sum, aes(Isolate, no_genes, fill=Type))+
  geom_bar(stat='identity')+
  facet_wrap(~Genus, scales='free_x')+
  scale_fill_manual(values = pal)+
  ylab("Number of BSG")+
  xlab("")+
  theme_bw()+
  scale_y_continuous(expand = c(0,0))


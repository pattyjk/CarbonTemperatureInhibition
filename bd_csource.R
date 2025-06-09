library(ggplot2)

#load pval data
pvals<-read.delim('~/Documents/GitHub/CarbonTemperatureInhibition/bd_csource_pvals.txt', header=T)

#correct pvals
pvals$corrP<-p.adjust(pvals$T.Test, method='BH')

#anything singificant
length(pvals$corrP>0.05)
pvals_not_sig<-pvals[which(pvals$corrP>0.05),]
View(pvals_not_sig)

pvals_sig<-pvals[-which(pvals$corrP>0.05),]
View(pvals_sig)

#plot Bd growth
bd_growth<-read.delim('~/Documents/GitHub/CarbonTemperatureInhibition/bd__csource_growth.txt', header=T)
ggplot(bd_growth, aes(SampleID, Per_growth, fill=SampleID))+
  geom_boxplot()+
  scale_fill_manual(values = c('white','white','white','white','white','white','white', 'blue', 'white', 'white','white','white','white','white','blue','blue','white','red','white','white','white','white','white','white','white','white','white','white','white','white','white', 'blue', 'white', 'white','white','white','white','white','blue','blue','white','red','white','white','white','white','white','white'))+
  facet_wrap(~BD_strain)+
  xlab("")+
  theme_bw()+
  coord_flip()+
  ylab("Percent Growth Relative to Positive Control")

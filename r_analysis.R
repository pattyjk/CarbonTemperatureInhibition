library(ggplot2)
library(ggpubr)

#read in Bd data
bd_data<-read.delim('Documents/GitHub/CarbonTemperatureInhibition/bd_data', header = T)

#split data by Bd strain
bd_split<-split(bd_data, f = bd_data$BdStrain)

#plot OD600 vs inhibition faceted by carbon source,for Jel 423 
ggplot(bd_split$Jel423, aes(OD600, Per_inhib, color=Strain))+
  geom_point()+
  stat_cor(method='pearson')+
  ggtitle('Jel 423')+
  facet_wrap(~Carbon_source)

#plot OD600 vs inhibition faceted by carbon source,for Jel 197
ggplot(bd_split$Jel197, aes(OD600, Per_inhib, color=Strain))+
  geom_point()+
  ggtitle('Jel197')+
  stat_cor(method='pearson')+
  facet_wrap(~Carbon_source)

#plot all data together
ggplot(bd_data, aes(OD600, Per_inhib, color=Strain))+
  geom_point()+
  stat_cor(method='pearson')+
  facet_wrap(~Carbon_source)

#plot OD600 vs inhibition faceted by carbon source,for Jel 423 
ggplot(bd_split$Jel423, aes(OD550, Per_inhib, color=Strain))+
  geom_point()+
  stat_cor(method='pearson')+
  ggtitle('Jel 423')+
  facet_wrap(~Carbon_source)

#plot OD600 vs inhibition faceted by carbon source,for Jel 197
ggplot(bd_split$Jel197, aes(OD550, Per_inhib, color=Strain))+
  geom_point()+
  ggtitle('Jel197')+
  stat_cor(method='pearson')+
  facet_wrap(~Carbon_source)

#plot all data together
ggplot(bd_data, aes(OD550, Per_inhib, color=Strain))+
  geom_point()+
  stat_cor(method='pearson')+
  facet_wrap(~Carbon_source)

#plot histograms to look at distribution of data
ggplot(bd_data, aes(Per_inhib, fill=Strain))+
  geom_histogram()+
  facet_wrap(~Carbon_source)

#plot OD600 vs OD550 to see how biofilms compare to growth
ggplot(bd_data, aes(OD600, OD550, color=Strain))+
  geom_point()+
  theme_bw()+
  stat_cor(method='spearman')+
  facet_wrap(~Carbon_source)

#biofilms over temperature for each isolate
ggplot(bd_data, aes(Temp, OD550, color=Strain))+
  geom_point()+
  theme_bw()+
  stat_cor(method='pearson')+
  facet_wrap(~Carbon_source)

#OD600 over temperature for each isolate
ggplot(bd_data, aes(Temp, OD600, color=Strain))+
  geom_point()+
  theme_bw()+
  geom_line()+
  stat_cor(method='pearson')+
  facet_wrap(~Carbon_source)

bd_split<-split(bd_data, f= bd_data$Strain)


glm_temp_HP<-glm(Temp + OD550 + OD600 ~ Per_inhib, data=bd_split$HP1C)
summary(glm_temp_HP)
anova(glm_temp_HP)

#load libraries 
library(ggplot2)
library(ggpubr)
library(dplyr)

#read in Bd data
bd_data<-read.delim('~/Documents/GitHub/CarbonTemperatureInhibition/bd_data', header = T)

#plot temperature range for each isolate
ggplot(bd_data, aes(Temp, OD600, color=Genus))+
  geom_point()+
  facet_wrap(~Strain)+
  scale_color_manual(values=c('grey', 'black', 'orange'))+
  theme_bw()+
  geom_smooth()+
  ylab("OD600")+
  xlab('Temperature')

#calculate p-values between Per_inhib and temperature
temp_inhib  <- bd_data %>%
  group_by(Strain, Carbon_source, BdStrain) %>%
  summarise(
    N = n(),
    cor_test = list(cor.test(Per_inhib, Temp, method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(
    Spearman_rho = sapply(cor_test, function(x) x$estimate),
    p_value = sapply(cor_test, function(x) x$p.value),
    Significant = p_value < 0.05
  ) %>%
  select(Strain, Carbon_source, BdStrain, Spearman_rho, p_value, Significant, N)
temp_inhib$comparison<-'inhib-temp'


#plot the Rho's
ggplot(temp_inhib, aes(Carbon_source, Spearman_rho, color=Significant, shape=BdStrain))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values=c('black', 'orange'))+
  facet_wrap(~Strain)+
  coord_flip()+
  ylab('Spearman Rho')+
  xlab('')

#calculate correlation between Inhibition and OD600
od600_inhib <- bd_data %>%
  group_by(Strain, Carbon_source, BdStrain) %>%
  summarise(
    N = n(),
    cor_test = list(cor.test(Per_inhib, OD600, method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(
    Spearman_rho = sapply(cor_test, function(x) x$estimate),
    p_value = sapply(cor_test, function(x) x$p.value),
    Significant = p_value < 0.05
  ) %>%
  select(Strain, Carbon_source, BdStrain, Spearman_rho, p_value, Significant, N)

od600_inhib$comparison<-'OD600-Inhib'

#calculate spearman correlation between inhibition and OD550
od550_inhib <- bd_data %>%
  group_by(Strain, Carbon_source, BdStrain) %>%
  summarise(
    N = n(),
    cor_test = list(cor.test(Per_inhib, OD550, method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(
    Spearman_rho = sapply(cor_test, function(x) x$estimate),
    p_value = sapply(cor_test, function(x) x$p.value),
    Significant = p_value < 0.05
  ) %>%
  select(Strain, Carbon_source, BdStrain, Spearman_rho, p_value, Significant, N)
od550_inhib$comparison<-'Inhib550'

#calculate spearman correlation between temp and OD550
od550_temp <- bd_data %>%
  group_by(Strain, Carbon_source, BdStrain) %>%
  summarise(
    N = n(),
    cor_test = list(cor.test(Temp, OD550, method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(
    Spearman_rho = sapply(cor_test, function(x) x$estimate),
    p_value = sapply(cor_test, function(x) x$p.value),
    Significant = p_value < 0.05
  ) %>%
  select(Strain, Carbon_source, BdStrain, Spearman_rho, p_value,Significant, N)
od550_temp$comparison<-'Temp550'

#calculate spearman correlation between temp and OD600
od600_temp <- bd_data %>%
  group_by(Strain, Carbon_source, BdStrain) %>%
  summarise(
    N = n(),
    cor_test = list(cor.test(Temp, OD600, method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(
    Spearman_rho = sapply(cor_test, function(x) x$estimate),
    p_value = sapply(cor_test, function(x) x$p.value),
    Significant = p_value < 0.05
  ) %>%
  select(Strain, Carbon_source, BdStrain, Spearman_rho, p_value,Significant, N)
od600_temp$comparison<-'Temp600'

#calculate spearman correlation between OD600 and OD550
od550_temp <- bd_data %>%
  group_by(Strain, Carbon_source, BdStrain) %>%
  summarise(
    N = n(),
    cor_test = list(cor.test(OD600, OD550, method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(
    Spearman_rho = sapply(cor_test, function(x) x$estimate),
    p_value = sapply(cor_test, function(x) x$p.value),
    Significant = p_value < 0.05
  ) %>%
  select(Strain, Carbon_source, BdStrain, Spearman_rho, p_value,Significant, N)

od550_temp$comparison<-'Temp550'

#bind all together & write to file
combined<-rbind(od600_temp, od550_temp, od600_inhib, od550_inhib, temp_inhib)
write.table(combined, '~/Documents/GitHub/CarbonTemperatureInhibition/correlations.txt', sep='\t', quote=F)




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



#calculate spearman corrs. for temperature and Bd inhibition
library(dplyr)
df<-bd_data
results_list <- list()

# Get all unique combinations
combinations <- df %>%
  distinct(BdStrain, Strain, Carbon_source)

# Loop through each combination and compute Spearman correlation and p-value
for (i in seq_len(nrow(combinations))) {
  sub_data <- df %>%
    filter(
      BdStrain == combinations$BdStrain[i],
      Strain == combinations$Strain[i],
      Carbon_source == combinations$Carbon_source[i]
    )
  
  # Ensure at least 3 non-NA observations for correlation
  if (nrow(sub_data) >= 3 && !all(is.na(sub_data$Temp)) && !all(is.na(sub_data$Per_inhib))) {
    test_result <- cor.test(sub_data$Temp, sub_data$Per_inhib, method = "spearman")
    results_list[[i]] <- data.frame(
      BdStrain = combinations$BdStrain[i],
      Strain = combinations$Strain[i],
      Carbon_source = combinations$Carbon_source[i],
      spearman_corr = test_result$estimate,
      p_value = test_result$p.value,
      n = nrow(sub_data)
    )
  }
}

# Combine results into a single data frame
results_df <- do.call(rbind, results_list)

# Adjust p-values for multiple testing (FDR correction)
results_df$fdr_p_value <- p.adjust(results_df$p_value, method = "fdr")

ggplot(results_df, aes(Carbon_source, spearman_corr, color=BdStrain))+
  geom_point()+
  coord_flip()+
  xlab("")+
  theme_bw()+
  scale_color_manual(values=c('#000', "grey"))+
  ylab("Spearman Correlation (Per. Inhib. - Temp")+
  facet_wrap(~Strain)

########Getting maximum growth with math!
#for OD600, fit loess curve and get max vlaue on curve (aka temperature optimum)
library(dplyr)
library(ggplot2)

df<-bd_data

# Optional: store predictions for each strain
all_loess_predictions <- data.frame()

for (strain in unique(df$Strain)) {
  sub_data <- df %>% filter(Strain == strain)
  
  if (nrow(sub_data) >= 4) {
    loess_fit <- loess(OD600 ~ Temp, data = sub_data, span = 0.75)
    
    temp_seq <- seq(min(sub_data$Temp), max(sub_data$Temp), length.out = 100)
    pred_od <- predict(loess_fit, newdata = data.frame(Temp = temp_seq))
    
    if (!all(is.na(pred_od))) {
      max_index <- which.max(pred_od)
      max_temp <- temp_seq[max_index]
      max_od600 <- pred_od[max_index]
      
      # Save maximum point
      loess_max_df <- rbind(loess_max_df, data.frame(
        Strain = strain,
        Max_Temp = max_temp,
        Max_OD600 = max_od600
      ))
      
      # Save full curve for visualization or output
      strain_preds <- data.frame(
        Strain = strain,
        Temp = temp_seq,
        Predicted_OD600 = pred_od
      )
      all_loess_predictions <- rbind(all_loess_predictions, strain_preds)
    }
  }
}

# View results
print(loess_max_df)

#####Loess probably not best, fit actual model
# Function to fit Brière model and extract parameters and curve max
library(dplyr)
library(nls.multstart)
library(readr)

# Brière model function
fit_briere <- function(sub_data) {
  tryCatch({
    fit <- nls_multstart(
      OD600 ~ a * Temp * (Temp - Tmin) * sqrt(Tmax - Temp),
      data = sub_data,
      iter = 500,
      start_lower = c(a = 0, Tmin = 0, Tmax = 20),
      start_upper = c(a = 1, Tmin = 5, Tmax = 45),
      supp_errors = "Y"
    )
    
    coefs <- coef(fit)
    
    # Predict over temp range
    temp_seq <- seq(min(sub_data$Temp), max(sub_data$Temp), by = 0.1)
    pred_od <- predict(fit, newdata = data.frame(Temp = temp_seq))
    
    max_idx <- which.max(pred_od)
    
    data.frame(
      Strain = unique(sub_data$Strain),
      Carbon_source = unique(sub_data$Carbon_source),
      Max_Temp = temp_seq[max_idx],
      Max_OD600 = pred_od[max_idx],
      a = coefs["a"],
      Tmin = coefs["Tmin"],
      Tmax = coefs["Tmax"]
    )
  }, error = function(e) {
    data.frame(
      Strain = unique(sub_data$Strain),
      Carbon_source = unique(sub_data$Carbon_source),
      Max_Temp = NA,
      Max_OD600 = NA,
      a = NA,
      Tmin = NA,
      Tmax = NA
    )
  })
}

# Run Brière model on each group
results <- df %>%
  group_by(Strain, Carbon_source) %>%
  group_split() %>%
  lapply(fit_briere) %>%
  bind_rows()

# Preview
print(results)

#summarize optimal temperature
library(plyr)

ggplot(results, aes(Strain, Max_Temp))+
  geom_boxplot()

temp_op_sum<-ddply(results, c('Strain'), summarize, mean=mean(Max_Temp), n=length(Max_Temp), sd=sd(Max_Temp), se=sd/n)

#plot with error
ggplot(temp_op_sum, aes(Strain, mean))+
  geom_point()+
  theme_bw()+
  coord_flip()+
  ylab('Average Temperature Optima')+
  xlab("")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)

##Calculate Briere model but ignore carbon source
fit_briere_strain <- function(sub_data) {
  tryCatch({
    fit <- nls_multstart(
      OD600 ~ a * Temp * (Temp - Tmin) * sqrt(Tmax - Temp),
      data = sub_data,
      iter = 500,
      start_lower = c(a = 0, Tmin = 0, Tmax = 20),
      start_upper = c(a = 1, Tmin = 5, Tmax = 45),
      supp_errors = "Y"
    )
    
    coefs <- coef(fit)
    
    # Predict across temperature range for plotting
    temp_seq <- seq(min(sub_data$Temp), max(sub_data$Temp), by = 0.1)
    pred_od <- predict(fit, newdata = data.frame(Temp = temp_seq))
    
    max_idx <- which.max(pred_od)
    max_temp <- temp_seq[max_idx]
    max_od <- pred_od[max_idx]
    
    # Return both prediction curve and summary
    curve_df <- data.frame(
      Strain = unique(sub_data$Strain),
      Temp = temp_seq,
      Predicted_OD600 = pred_od
    )
    
    summary_df <- data.frame(
      Strain = unique(sub_data$Strain),
      Max_Temp = max_temp,
      Max_OD600 = max_od,
      a = coefs["a"],
      Tmin = coefs["Tmin"],
      Tmax = coefs["Tmax"]
    )
    
    list(curve = curve_df, summary = summary_df)
  }, error = function(e) {
    list(
      curve = data.frame(
        Strain = unique(sub_data$Strain),
        Temp = NA,
        Predicted_OD600 = NA
      ),
      summary = data.frame(
        Strain = unique(sub_data$Strain),
        Max_Temp = NA,
        Max_OD600 = NA,
        a = NA,
        Tmin = NA,
        Tmax = NA
      )
    )
  })
}

# Apply to each unique Strain
results_list <- df %>%
  group_by(Strain) %>%
  group_split() %>%
  lapply(fit_briere_strain)

# Combine outputs
curves_df <- bind_rows(lapply(results_list, `[[`, "curve"))
summaries_df <- bind_rows(lapply(results_list, `[[`, "summary"))

# ✅ Plot Brière curves using ggplot2
ggplot(curves_df, aes(x = Temp, y = Predicted_OD600, color = Strain)) +
  geom_line(size = 1) +
  labs(
    title = "Brière Model Fit per Strain",
    x = "Temperature (°C)",
    y = "Predicted OD600"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

#bind back to original data frame for more math fun
bd_with_optima <- df %>%
  left_join(results, by = c("Strain", "Carbon_source"))

#calculate a temperature differential from incubation temperature and optima
bd_with_optima$delta_t<-bd_with_optima$Temp - bd_with_optima$Max_Temp

#histogram for funzies, looks nice
hist(bd_with_optima$delta_t)

#write to file
write.table(bd_with_optima, '~/Documents/GitHub/CarbonTemperatureInhibition/bd_data_optima.txt', row.names=F, quote=F, sep='\t')

#calculate correlations between deviation form optima and inhibition
temp_inhib_optima  <- bd_with_optima %>%
  group_by(Strain, Carbon_source, BdStrain) %>%
  summarise(
    N = n(),
    cor_test = list(cor.test(Per_inhib, delta_t, method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(
    Spearman_rho = sapply(cor_test, function(x) x$estimate),
    p_value = sapply(cor_test, function(x) x$p.value),
    Significant = p_value < 0.05
  ) %>%
  select(Strain, Carbon_source, BdStrain, Spearman_rho, p_value, Significant, N)

#view and plot data
View(temp_inhib_optima)
library(viridis)
ggplot(temp_inhib_optima, aes(Strain, Carbon_source, fill=Spearman_rho))+
  geom_tile()+
  xlab('')+
  scale_fill_viridis(option = "D", name = "Spearman Rho")+
  ylab('')+
  facet_wrap(~BdStrain)


ggplot(temp_inhib_optima, aes(Carbon_source, Spearman_rho, color=Significant, shape=BdStrain))+
  facet_wrap(~Strain)+
  theme_bw()+
  geom_point()+
  scale_color_manual(values = c('black', 'orange'))+
  ylab('Spearman Rho')+
  xlab("")+
  coord_flip()
  

#split data by strain
bd_with_optima_split<-split(bd_with_optima, bd_with_optima$Strain)

library(ggpubr)
ggplot(bd_with_optima_split$`P060114-29-3-1`, aes(delta_t, Per_inhib, color=BdStrain))+
  geom_point()+
  scale_color_manual(values=c('orange', 'black'))+
  theme_bw()+
  ggtitle('P060114-29-3-1 Pseudomonas')+
 stat_cor(method = 'pearson')+
  facet_wrap(~Carbon_source)+
  geom_smooth(method='lm')+
  ylab("Percent Inhibition of Bd")+
  xlab("Change in temperature from calculated optima")

ggplot(bd_with_optima_split$HP1C, aes(delta_t, Per_inhib, color=BdStrain))+
  geom_point()+
  scale_color_manual(values=c('orange', 'black'))+
  theme_bw()+
  ggtitle('HP1C Pseudomonas')+
  stat_cor(method = 'pearson')+
  facet_wrap(~Carbon_source)+
  geom_smooth(method='lm')+
  ylab("Percent Inhibition of Bd")+
  xlab("Change in temperature from calculated optima")

ggplot(bd_with_optima_split$Panama13O, aes(delta_t, Per_inhib, color=BdStrain))+
  geom_point()+
  scale_color_manual(values=c('orange', 'black'))+
  theme_bw()+
  ggtitle('Panama130 Chryseobacterium')+
  stat_cor(method = 'pearson')+
  facet_wrap(~Carbon_source)+
  geom_smooth(method='lm')+
  ylab("Percent Inhibition of Bd")+
  xlab("Change in temperature from calculated optima")

ggplot(bd_with_optima_split$THA3.2, aes(delta_t, Per_inhib, color=BdStrain))+
  geom_point()+
  scale_color_manual(values=c('orange', 'black'))+
  theme_bw()+
  ggtitle('THA3.2 Chryseobacterium')+
  stat_cor(method = 'pearson')+
  facet_wrap(~Carbon_source)+
  geom_smooth(method='lm')+
  ylab("Percent Inhibition of Bd")+
  xlab("Change in temperature from calculated optima")

ggplot(bd_with_optima_split$Panama68B, aes(delta_t, Per_inhib, color=BdStrain))+
  geom_point()+
  scale_color_manual(values=c('orange', 'black'))+
  theme_bw()+
  ggtitle('Panama68B Stenotrophomonas')+
  stat_cor(method = 'pearson')+
  facet_wrap(~Carbon_source)+
  geom_smooth(method='lm')+
  ylab("Percent Inhibition of Bd")+
  xlab("Change in temperature from calculated optima")

ggplot(bd_with_optima_split$RSM3.2, aes(delta_t, Per_inhib, color=BdStrain))+
  geom_point()+
  scale_color_manual(values=c('orange', 'black'))+
  theme_bw()+
  ggtitle('RSM3.2 Stenotrophomonas')+
  stat_cor(method = 'pearson')+
  facet_wrap(~Carbon_source)+
  geom_smooth(method='lm')+
  ylab("Percent Inhibition of Bd")+
  xlab("Change in temperature from calculated optima")



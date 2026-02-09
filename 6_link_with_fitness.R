
# Load packages
library(ggplot2)
library(ggdist)
library(robustHD)
library(brms)
library(rstan)
library(StanHeaders)
library(bayesplot)
library(robustHD)
library(cmdstanr)
library(RColorBrewer)
library(ggplot2)

# set the location of the data
dir_data_BBAL = "C:/Users/bmohring/Documents/GitHub/BBAL-Senescence/Data/"
dir_data = "C:/Users/bmohring/Documents/GitHub/BBAL-Senescence/"


# Open dataset
df_AgeStd_randomSlopesAndInterceptsIDlevel = read.csv( paste0(dir_data, "data_BBAL_senescence_estimates_priors.csv"))
 
df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort = substr(df_AgeStd_randomSlopesAndInterceptsIDlevel$cohortPop,1,4)
df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort =as.numeric(df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort )

# subset of dead individuals
SubsetDead_cohortBefore2000 = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$Year_last_seen <= 2020  &
                                                                                                           df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort <2000),]

table(SubsetDead_cohortBefore2000$population)


SubsetDead_cohortBefore2000_noOutlier  =SubsetDead_cohortBefore2000[which( SubsetDead_cohortBefore2000$OnsetSenesc_from_backTransformedCoefs<40),]
SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs = -SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs

# Correlation between onset and rate of senescence
cor(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen")],
    SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen")])
cor(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island")],
    SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island")])


#### Figure 4 #### 

library(scales)

Fig_4_A = ggplot(data = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island"),], 
                 aes(x =SenesceRate_from_backTransformedCoefs , y=OnsetSenesc_from_backTransformedCoefs , col = LRS))+
  geom_point(size = 1.5, alpha = 0.6)+
  scale_colour_viridis_c(begin = 0.85, end = 0.15, option = "inferno", trans = pseudo_log_trans(base = 10), limits = c(0, 16))+
  theme_classic()+ xlim(c(0.01,0.06))+ ylim(c(17,32))+ 
  xlab("Senescence rate") + ylab("Age at onset of senescence")+
  theme( legend.text = element_text(size=8) , legend.key.size = unit(0.55, "cm"), legend.title = element_text(size=9)  )
  # theme( legend.text = element_text(size=9) , legend.title = element_text(size=11),legend.position = c(0.95, 0.5) )
Fig_4_A

 
Fig_4_B = ggplot(data = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen"  ),], 
                 aes(x =SenesceRate_from_backTransformedCoefs , y=OnsetSenesc_from_backTransformedCoefs , col = LRS))+
  geom_point(size = 1.5, alpha = 0.6)+
  scale_colour_viridis_c(begin = 0.90, end = 0.05, option = "mako", trans = pseudo_log_trans(base = 10))+
  theme_classic()+ xlim(c(0.01,0.06))+ ylim(c(17,32))+
  xlab("Senescence rate") + ylab("Age at onset of senescence")+
  theme( legend.text = element_text(size=8) , legend.key.size = unit(0.55, "cm"), legend.title = element_text(size=9)  )
# theme( legend.text = element_text(size=9), legend.title = element_text(size=11),legend.position = c(0.99, 0.5) )
  # theme( legend.text = element_text(size=9), legend.title = element_text(size=11) )
Fig_4_B

 
cowplot::plot_grid(Fig_4_A, Fig_4_B)


# std variables
SubsetDead_cohortBefore2000_noOutlier$OnsetSenescence_std = robustHD::standardize(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs)
SubsetDead_cohortBefore2000_noOutlier$SenesceRate_std = robustHD::standardize(SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs)

prior_fitness <- c(
  prior(normal(0, 5), class = "Intercept"),
  prior(normal(0, 2), class = "b")
  
)


GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence  = brm(LRS   ~   OnsetSenescence_std * population  ,
                                                 data=SubsetDead_cohortBefore2000_noOutlier ,
                                                 family = poisson   , 
                                                 warmup = 2000,
                                                 iter = 4000,
                                                 chains =8,
                                                 cores = 8,
                                                 prior = prior_fitness,
                                                 backend = "cmdstanr" ,
                                                 control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence) 

bayes_R2(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence)
# pp_check(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence)
# plot(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence)


SubsetDead_cohortBefore2000_noOutlier$population=as.factor(SubsetDead_cohortBefore2000_noOutlier$population)
SubsetDead_cohortBefore2000_noOutlier$population_refKer = relevel(SubsetDead_cohortBefore2000_noOutlier$population , ref=2)

GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer  = brm(LRS   ~   OnsetSenescence_std * population_refKer  ,
                                                 data=SubsetDead_cohortBefore2000_noOutlier ,
                                                 family = poisson   ,
                                                 warmup = 2000,
                                                 iter = 4000,
                                                 chains =8,
                                                 cores = 8,
                                                 prior =prior_fitness,
                                                 backend = "cmdstanr" ,
                                                 control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer) 

# bayes_R2(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer) 
# pp_check(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer)
# plot(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer)

library(dplyr)

ranges_onset <- SubsetDead_cohortBefore2000_noOutlier %>%
  group_by(population) %>%  summarise(    xmin = min(OnsetSenescence_std),    xmax = max(OnsetSenescence_std)  )

newdata_onset <- ranges_onset %>%
  rowwise() %>%
  do(    data.frame(      population = .$population,    OnsetSenescence_std = seq(.$xmin, .$xmax, length.out = 100)    ))
  
preds_onset <- fitted(  GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence,
  newdata = newdata_onset,  summary = TRUE)

newdata_onset$Estimate <- preds_onset[, "Estimate"]
newdata_onset$Q2.5     <- preds_onset[, "Q2.5"]
newdata_onset$Q97.5    <- preds_onset[, "Q97.5"]


newdata_onset$OnsetSenesc_plot <-
  newdata_onset$OnsetSenescence_std *
  sd(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs) +
  mean(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs)


Fig_4_C <- ggplot() +
  geom_ribbon(    data = newdata_onset, aes(x = OnsetSenesc_plot, ymin = Q2.5, ymax = Q97.5, fill = population    ),
    alpha = 0.3  ) +
  geom_line(    data = newdata_onset, aes(x = OnsetSenesc_plot, y = Estimate, color = population    ),
    linewidth = 1  ) +
  geom_jitter(    data = SubsetDead_cohortBefore2000_noOutlier,    aes(      x = OnsetSenesc_from_backTransformedCoefs,      y = LRS,      color = population    ),
    width = 0,    height = 0.2,    alpha = 0.5,    size = 1.5  ) +
  scale_fill_manual(values = c("darkorange", "turquoise4")) +
  scale_color_manual(values = c("darkorange", "turquoise4")) +
  theme_classic() +
  xlab("Age at onset of senescence") +
  ylab("LRS") +
  theme(legend.position = "none")

Fig_4_C

GLM_KerBI_Poisson_BRM_LRS_SenesceRate  = brm(LRS   ~   SenesceRate_std * population  ,
                                             data=SubsetDead_cohortBefore2000_noOutlier ,
                                             family = poisson   ,
                                             warmup = 2000,
                                             iter = 4000,
                                             chains =8,
                                             cores = 8,
                                             prior =prior_fitness,
                                             backend = "cmdstanr" ,
                                             control = list(adapt_delta = 0.85)
)

summary(GLM_KerBI_Poisson_BRM_LRS_SenesceRate) 
bayes_R2(GLM_KerBI_Poisson_BRM_LRS_SenesceRate)


GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer  = brm(LRS   ~   SenesceRate_std * population_refKer  ,
                                             data=SubsetDead_cohortBefore2000_noOutlier ,
                                             family = poisson   ,
                                             warmup = 2000,
                                             iter = 4000,
                                             chains =8,
                                             cores = 8,
                                             prior =prior_fitness,
                                             backend = "cmdstanr" ,
                                             control = list(adapt_delta = 0.85)
)

summary(GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer) 

ranges <- SubsetDead_cohortBefore2000_noOutlier %>%
  group_by(population) %>%
  summarise(    xmin = min(SenesceRate_std),    xmax = max(SenesceRate_std)
  )
newdata <- ranges %>%
  rowwise() %>%
  do(    data.frame(      population = .$population, SenesceRate_std = seq(.$xmin, .$xmax, length.out = 100)
    )
  )
preds <- fitted(  GLM_KerBI_Poisson_BRM_LRS_SenesceRate,  newdata = newdata,  summary = TRUE
)

newdata$Estimate <- preds[, "Estimate"]
newdata$Q2.5     <- preds[, "Q2.5"]
newdata$Q97.5    <- preds[, "Q97.5"]

newdata$SenesceRate_plot <-
  newdata$SenesceRate_std *
  sd(SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs) +
  mean(SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs)


Fig_4_D = ggplot() +
  geom_ribbon(    data = newdata,    aes( x = SenesceRate_plot, ymin = Q2.5, ymax = Q97.5, fill = population    ),
    alpha = 0.3  ) +
  geom_line(    data = newdata,    aes( x = SenesceRate_plot, y = Estimate, color = population    ),
    linewidth = 1  ) +
  geom_jitter(    data = SubsetDead_cohortBefore2000_noOutlier,   
                  aes(  x = SenesceRate_from_backTransformedCoefs, y = LRS, color = population ),
    width = 0,    height = 0.2,    alpha = 0.5,    size = 1.5) +
  scale_fill_manual(values = c("darkorange", "turquoise4")) +
  scale_color_manual(values = c("darkorange", "turquoise4")) +
  theme_classic() +
  xlab("Senescence rate") +
  ylab("LRS") +
  theme(legend.position = "none")
Fig_4_D

cowplot::plot_grid(Fig_4_A,
                   Fig_4_B, 
                   Fig_4_C,
                   Fig_4_D, 
                   align="hv",
                   ncol = 2, labels = c("A", "B", "C", "D"))


ggsave(paste0(dir_data, "Figure4.png" ))




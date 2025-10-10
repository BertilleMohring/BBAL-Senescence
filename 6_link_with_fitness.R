
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
dir_data = "C:/Users/bmohring//"

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



Fig_4_C  = plot(conditional_effects(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence), ask = FALSE)
Fig_4_C=Fig_4_C$`OnsetSenescence_std:population`+  
  theme_classic() + 
  scale_fill_manual(values = c("darkorange", "turquoise4")) + 
  scale_color_manual(values = c("darkorange", "turquoise4")) + 
  geom_jitter(inherit.aes = FALSE, size =1.5,
              data = SubsetDead_cohortBefore2000_noOutlier, 
              width = 0, height = 0.2, alpha = 0.5,
              aes( x= OnsetSenescence_std,  y  = LRS, group  = population , col = population ))+
  # geom_pointrange(inherit.aes = FALSE,
  #                 data = data_summary, 
  #                 aes( x= Age,  y  = mean_ReproS, 
  #                      ymin = mean_ReproS -se_ReproS  , ymax  = mean_ReproS+se_ReproS ,
  #                      group  = population , col = population ), fatten = 2)
  xlab("Age at onset of senescence") + ylab("LRS") +
  # geom_vline(xintercept = min(SubsetDead_cohortBefore2000_noOutlier_Ker$OnsetSenescence), linetype = "dashed", color ="turquoise4" )+
  # geom_vline(xintercept = max(SubsetDead_cohortBefore2000_noOutlier_Ker$OnsetSenescence), linetype = "dashed", color ="turquoise4" )+
  # geom_vline(xintercept = min(SubsetDead_cohortBefore2000_noOutlier_BI$OnsetSenescence), linetype = "dashed", color ="darkorange" )+
  # geom_vline(xintercept = max(SubsetDead_cohortBefore2000_noOutlier_BI$OnsetSenescence), linetype = "dashed", color ="darkorange" )+
  theme(legend.position = "none"  )

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

# pp_check(GLM_KerBI_Poisson_BRM_LRS_SenesceRate)

Fig_4_D  = plot(conditional_effects(GLM_KerBI_Poisson_BRM_LRS_SenesceRate), ask = FALSE)
Fig_4_D = Fig_4_D$`SenesceRate_std:population`+  
  theme_classic() + 
  scale_fill_manual(values = c("darkorange", "turquoise4")) + 
  scale_color_manual(values = c("darkorange", "turquoise4")) + 
  geom_jitter(inherit.aes = FALSE, size =1.5,
              data = SubsetDead_cohortBefore2000_noOutlier, 
              width = 0, height = 0.2, alpha = 0.5,
              aes( x= SenesceRate_std,  y  = LRS, group  = population , col = population ))+
  # geom_pointrange(inherit.aes = FALSE,
  #                 data = data_summary, 
  #                 aes( x= Age,  y  = mean_ReproS, 
  #                      ymin = mean_ReproS -se_ReproS  , ymax  = mean_ReproS+se_ReproS ,
  #                      group  = population , col = population ), fatten = 2)
  xlab("Senescence rate") + ylab("LRS")+
  # geom_vline(xintercept = min(SubsetDead_cohortBefore2000_noOutlier_Ker$SenesceRate), linetype = "dashed", color ="turquoise4" )+
  # geom_vline(xintercept = max(SubsetDead_cohortBefore2000_noOutlier_Ker$SenesceRate), linetype = "dashed", color ="turquoise4" )+
  # geom_vline(xintercept = min(SubsetDead_cohortBefore2000_noOutlier_BI$SenesceRate), linetype = "dashed", color ="darkorange" )+
  # geom_vline(xintercept = max(SubsetDead_cohortBefore2000_noOutlier_BI$SenesceRate), linetype = "dashed", color ="darkorange" )+
  theme(legend.position = "none"  )
Fig_4_D

cowplot::plot_grid(Fig_4_A,
                   Fig_4_B, 
                   Fig_4_C,
                   Fig_4_D, 
                   ncol = 2, labels = c("A", "B", "C", "D"))










Fig_4_C_v2 =  ggplot(data = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen"),], 
                     aes(x= OnsetSenesc_from_backTransformedCoefs,  y  = LRS ))+
  geom_jitter(inherit.aes = FALSE, size =1.5,
              data = SubsetDead_cohortBefore2000_noOutlier, 
              width = 0, height = 0.2, alpha = 0.4,
              aes( x= OnsetSenesc_from_backTransformedCoefs,  y  = LRS, group  = population , col = population ))+
  scale_fill_manual(values = c("darkorange", "turquoise4")) + 
  scale_color_manual(values = c("darkorange", "turquoise4")) + 
  geom_smooth(method = "glm", method.args = list(family = "poisson"), col ="turquoise4", fill ="turquoise4" )+ 
  geom_smooth(inherit.aes = FALSE,  
              data = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island"),],
              aes(x= OnsetSenesc_from_backTransformedCoefs,  y  = LRS ),
              method = "glm", method.args = list(family = "poisson"), col ="darkorange", fill ="darkorange" )+
  xlab("Age at onset of senescence") + ylab("LRS")+
  theme_classic()+
  theme(legend.position = "none"  )

Fig_4_C_v2
 
Fig_4_D_v2 =  ggplot(data = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen"),], 
                     aes(x= SenesceRate_from_backTransformedCoefs,  y  = LRS ))+
  geom_jitter(inherit.aes = FALSE, size =1.5,
              data = SubsetDead_cohortBefore2000_noOutlier, 
              width = 0, height = 0.2, alpha = 0.4,
              aes( x= SenesceRate_from_backTransformedCoefs,  y  = LRS, group  = population , col = population ))+
  scale_fill_manual(values = c("darkorange", "turquoise4")) + 
  scale_color_manual(values = c("darkorange", "turquoise4")) + 
  geom_smooth(method = "glm", method.args = list(family = "poisson"), col ="turquoise4", fill ="turquoise4" )+ 
  geom_smooth(inherit.aes = FALSE,  
              data = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island"),],
              aes(x= SenesceRate_from_backTransformedCoefs,  y  = LRS ),
              method = "glm", method.args = list(family = "poisson"), col ="darkorange", fill ="darkorange" )+
  xlab("Senescence rate") + ylab("LRS")+
  theme_classic()+
  theme(legend.position = "none"  )

Fig_4_D_v2

cowplot::plot_grid(Fig_4_A,
                   Fig_4_B, 
                   Fig_4_C_v2,
                   Fig_4_D_v2, 
                   align="hv",
                   ncol = 2, labels = c("A", "B", "C", "D"),  label_size = 12)



ggsave(paste0(dir_data, "Appendix4_Figure1.png" ))

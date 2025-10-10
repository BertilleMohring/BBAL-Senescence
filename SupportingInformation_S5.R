

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


df_AgeStd_randomSlopesAndInterceptsIDlevel = read.csv(paste0(dir_data, "data_BBAL_senescence_estimates_priors.csv"))

df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort = substr(df_AgeStd_randomSlopesAndInterceptsIDlevel$cohortPop,1,4)
df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort =as.numeric(df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort )

# subset of dead individuals
SubsetDead_cohortBefore2000 = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$Year_last_seen <= 2020  &
                                                                                 df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort <2000),]

table(SubsetDead_cohortBefore2000$population)


SubsetDead_cohortBefore2000_noOutlier  =SubsetDead_cohortBefore2000[which( SubsetDead_cohortBefore2000$OnsetSenesc_from_backTransformedCoefs<40),]

cor(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen")],
    SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen")])
cor(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island")],
    SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island")])


#### Dead and alive ####
df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot = "Alive"
df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$Year_last_seen <= 2020  )]= "Dead"



table(df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot , df_AgeStd_randomSlopesAndInterceptsIDlevel$population)

df_AgeStd_randomSlopesAndInterceptsIDlevel$pop_DeadOrNot = paste0(df_AgeStd_randomSlopesAndInterceptsIDlevel$population,
                                                                  "_",
                                                                  df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot)




Fig_4_C_v2 =  ggplot(data = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs > 0 &
                                                                               df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),] , 
                     aes(x= OnsetSenesc_from_backTransformedCoefs,  y  = LRS, group  = pop_DeadOrNot , col = pop_DeadOrNot ))+
  geom_jitter(inherit.aes = FALSE, size =1.3,
              data = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs > 0 &
                                                                        df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),] , 
              width = 0, height = 0.2, alpha = 0.4,
              aes( x= OnsetSenesc_from_backTransformedCoefs,  y  = LRS, group  = pop_DeadOrNot , col = pop_DeadOrNot ))+
  scale_fill_manual(values = c("darkorange","darkorange4", "turquoise2", "turquoise4")) + 
  scale_color_manual(values = c("darkorange","darkorange4", "turquoise2", "turquoise4")) + 
  geom_smooth(method = "glm", method.args = list(family = "poisson"))+ 
  # geom_smooth(inherit.aes = FALSE,  
  #             data = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population=="Bird island"),],
  #             aes(x= OnsetSenesc_from_backTransformedCoefs,  y  = LRS , group  = pop_DeadOrNot , col = pop_DeadOrNot),
  #             method = "glm", method.args = list(family = "poisson"), col ="darkorange", fill ="darkorange" )+
  xlab("Age at onset of senescence") + ylab("LRS")+
  theme_classic()+
  theme(legend.position = "none"  )

Fig_4_C_v2

Fig_4_D_v2 =  ggplot(data = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs > 0 &
                                                                               df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),] , 
                     aes(x= -SenesceRate_from_backTransformedCoefs,  y  = LRS , group  = pop_DeadOrNot , col = pop_DeadOrNot ))+
  geom_jitter(inherit.aes = FALSE, size =1.3,
              data = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs > 0 &
                                                                        df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),], 
              width = 0, height = 0.2, alpha = 0.4,
              aes( x= -SenesceRate_from_backTransformedCoefs,  y  = LRS, group  = pop_DeadOrNot , col = pop_DeadOrNot ))+
  scale_fill_manual(values = c("darkorange","darkorange4", "turquoise2", "turquoise4")) + 
  scale_color_manual(values = c("darkorange","darkorange4", "turquoise2", "turquoise4")) + 
  geom_smooth(method = "glm", method.args = list(family = "poisson"))+ 
  # geom_smooth(inherit.aes = FALSE,  
  #             data = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population=="Bird island"),],
  #             aes(x= SenesceRate_from_backTransformedCoefs,  y  = LRS , group  = pop_DeadOrNot , col = pop_DeadOrNot),
  #             method = "glm", method.args = list(family = "poisson"), col ="darkorange", fill ="darkorange" )+
  xlab("Senescence rate") + ylab("LRS")+
  theme_classic()+
  theme(legend.position = "none"  )

Fig_4_D_v2


cowplot::plot_grid( Fig_4_C_v2,
                    Fig_4_D_v2, 
                    align="hv",
                    ncol = 2, labels = c("A", "B"),  label_size = 12)

ggsave(paste0(dir_data, "Appendix5_Figure1.png"))





df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs > 0 &
                                                                                                          df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),] 

df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenescence_std = robustHD::standardize(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenesc_from_backTransformedCoefs)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_std = robustHD::standardize(-df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_from_backTransformedCoefs)

table(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot, df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population)

# set prior
prior_fitness <- c(
  prior(normal(0, 5), class = "Intercept"),
  prior(normal(0, 2), class = "b")
  
)


GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence  = brm(LRS   ~   OnsetSenescence_std * population * DeadOrNot ,
                                                        data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
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

df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population=as.factor(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population_refKer = relevel(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population , ref=2)

GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer  = brm(LRS   ~    OnsetSenescence_std * population_refKer * DeadOrNot ,
                                                               data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                               family = poisson   ,
                                                        warmup = 2000,
                                                        iter = 4000,
                                                        chains =8,
                                                        cores = 8,
                                                        prior = prior_fitness,
                                                               backend = "cmdstanr" ,
                                                               control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer) 

df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot=as.factor(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot_refDead = relevel(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot , ref=2)


GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence  = brm(LRS   ~   OnsetSenescence_std * population * DeadOrNot_refDead ,
                                                        data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
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



GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer  = brm(LRS   ~    OnsetSenescence_std * population_refKer * DeadOrNot_refDead ,
                                                               data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                               family = poisson   ,
                                                        warmup = 2000,
                                                        iter = 4000,
                                                        chains =8,
                                                        cores = 8,
                                                        prior = prior_fitness,
                                                               backend = "cmdstanr" ,
                                                               control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer) 


# Senescence rate

GLM_KerBI_Poisson_BRM_LRS_SenesceRate  = brm(LRS   ~   SenesceRate_std * population * DeadOrNot ,
                                                    data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                    family = poisson   ,
                                             warmup = 2000,
                                             iter = 4000,
                                             chains =8,
                                             cores = 8,
                                             prior = prior_fitness,
                                                    backend = "cmdstanr" ,
                                                    control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_LRS_SenesceRate) 
bayes_R2(GLM_KerBI_Poisson_BRM_LRS_SenesceRate) 


df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population=as.factor(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population_refKer = relevel(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population , ref=2)

GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer  = brm(LRS   ~    SenesceRate_std * population_refKer * DeadOrNot ,
                                                           data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                           family = poisson   ,
                                                    warmup = 2000,
                                                    iter = 4000,
                                                    chains =8,
                                                    cores = 8,
                                                    prior = prior_fitness,
                                                           backend = "cmdstanr" ,
                                                           control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer) 

df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot=as.factor(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot_refDead = relevel(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot , ref=2)


GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refDead  = brm(LRS   ~   SenesceRate_std * population * DeadOrNot_refDead ,
                                                            data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                            family = poisson   ,
                                                     warmup = 2000,
                                                     iter = 4000,
                                                     chains =8,
                                                     cores = 8,
                                                     prior = prior_fitness,
                                                            backend = "cmdstanr" ,
                                                            control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refDead) 



GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer_refDead  = brm(LRS   ~    SenesceRate_std * population_refKer * DeadOrNot_refDead ,
                                                                   data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                                   family = poisson   ,
                                                            warmup = 2000,
                                                            iter = 4000,
                                                            chains =8,
                                                            cores = 8,
                                                            prior = prior_fitness,
                                                                   backend = "cmdstanr" ,
                                                                   control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer_refDead) 

pp_check(GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer_refDead)


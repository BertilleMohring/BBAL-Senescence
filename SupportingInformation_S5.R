#### CODE USED TO PRODUCE SUPPORTING INFORMATION S5 ####
#### Supporting Information S5: Robustness test for the fitness outcomes of the life-history strategies exhibited by black-browed albatrosses ####  
#### using all available data #### 

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



df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs > 0 &
                                                                                                          df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),] 

df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenescence_std = robustHD::standardize(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenesc_from_backTransformedCoefs)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_from_backTransformedCoefs= -df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_from_backTransformedCoefs
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_std = robustHD::standardize(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_from_backTransformedCoefs)

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


GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refDead  = brm(LRS   ~   OnsetSenescence_std * population * DeadOrNot_refDead ,
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
summary(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refDead) 



GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer_refDead  = brm(LRS   ~    OnsetSenescence_std * population_refKer * DeadOrNot_refDead ,
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
summary(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer_refDead) 


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

# pp_check(GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer_refDead)



# Make plots

ranges_onset_v2 <- df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier %>%
  group_by(population, DeadOrNot, pop_DeadOrNot) %>%
  summarise(    xmin = min(OnsetSenescence_std),    xmax = max(OnsetSenescence_std),
                .groups = "drop"
  )

newdata_onset_v2 <- ranges_onset_v2 %>%
  rowwise() %>%
  do(data.frame(    population = .$population,    DeadOrNot  = .$DeadOrNot,
                    pop_DeadOrNot = .$pop_DeadOrNot,
                    OnsetSenescence_std = seq(.$xmin, .$xmax, length.out = 100)
  )) 


preds_onset_v2 <- fitted(  GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence,
                           newdata = newdata_onset_v2,  summary = TRUE
)

newdata_onset_v2$Estimate <- preds_onset_v2[, "Estimate"]
newdata_onset_v2$Q2.5     <- preds_onset_v2[, "Q2.5"]
newdata_onset_v2$Q97.5    <- preds_onset_v2[, "Q97.5"]

newdata_onset_v2$OnsetSenesc_plot <-
  newdata_onset_v2$OnsetSenescence_std *
  sd(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenesc_from_backTransformedCoefs) +
  mean(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenesc_from_backTransformedCoefs)

Fig_S5_1_C <- ggplot() +
  geom_jitter(
    data = df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier,
    aes(      x = OnsetSenesc_from_backTransformedCoefs,      y = nb_BA_life,
              color = pop_DeadOrNot    ),    width = 0,    height = 0.2,    alpha = 0.4,    size = 1.5  ) +
  geom_ribbon(    data = newdata_onset_v2,    aes(      x = OnsetSenesc_plot,      ymin = Q2.5,      ymax = Q97.5,      fill = pop_DeadOrNot
  ),    alpha = 0.25  ) +
  geom_line(    data = newdata_onset_v2,    aes(      x = OnsetSenesc_plot,      y = Estimate,      color = pop_DeadOrNot
  ),    linewidth = 1  ) +
  scale_color_manual(values = c(    "darkorange", "darkorange4",    "turquoise2", "turquoise4"  )) +
  scale_fill_manual(values = c(    "darkorange", "darkorange4",    "turquoise2", "turquoise4"  )) +
  theme_classic() +
  xlab("Age at onset of senescence") +  ylab("LRS") +
  theme(legend.position = "none")

Fig_S5_1_C


ranges_rate_v2 <- df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier %>%
  group_by(population, DeadOrNot, pop_DeadOrNot) %>%
  summarise(    xmin = min(SenesceRate_std),    xmax = max(SenesceRate_std),    .groups = "drop"  )
newdata_rate_v2 <- ranges_rate_v2 %>%
  rowwise() %>%
  do(data.frame(    population = .$population,    DeadOrNot  = .$DeadOrNot,    pop_DeadOrNot = .$pop_DeadOrNot,
                    SenesceRate_std = seq(.$xmin, .$xmax, length.out = 100)
  ))

preds_rate_v2 <- fitted(  GLM_KerBI_Poisson_BRM_LRS_SenesceRate,
                          newdata = newdata_rate_v2,
                          summary = TRUE
)

newdata_rate_v2$Estimate <- preds_rate_v2[, "Estimate"]
newdata_rate_v2$Q2.5     <- preds_rate_v2[, "Q2.5"]
newdata_rate_v2$Q97.5    <- preds_rate_v2[, "Q97.5"]

newdata_rate_v2$SenesceRate_plot <-   newdata_rate_v2$SenesceRate_std *
  sd(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_from_backTransformedCoefs) +
  mean(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_from_backTransformedCoefs)


Fig_S5_1_D <- ggplot() +
  geom_jitter(    data = df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier,
                  aes(      x = SenesceRate_from_backTransformedCoefs,      y = nb_BA_life,      color = pop_DeadOrNot    ),
                  width = 0,    height = 0.2,    alpha = 0.4,    size = 1.5
  ) +
  geom_ribbon(    data = newdata_rate_v2,
                  aes(      x = SenesceRate_plot,      ymin = Q2.5,      ymax = Q97.5,      fill = pop_DeadOrNot
                  ),    alpha = 0.25  ) +
  geom_line(    data = newdata_rate_v2,
                aes(      x = SenesceRate_plot,      y = Estimate,      color = pop_DeadOrNot    ),    linewidth = 1  ) +
  scale_color_manual(values = c(    "darkorange", "darkorange4",    "turquoise2", "turquoise4"  )) +
  scale_fill_manual(values = c(    "darkorange", "darkorange4",    "turquoise2", "turquoise4"  )) +
  theme_classic() +
  xlab("Senescence rate") +
  ylab("LRS") +
  theme(legend.position = "none")

Fig_S5_1_D


cowplot::plot_grid( Fig_S5_1_C,
                    Fig_S5_1_D, 
                    align="hv",
                    ncol = 2, labels = c("A", "B"),  label_size = 12)


# ggsave(paste0(dir_data, "SupportingInformation_S4_Figure_S4_2.png"))

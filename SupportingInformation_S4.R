#### CODE USED TO PRODUCE SUPPORTING INFORMATION S4 ####
#### Supporting Information S4: Relationship between number of breeding attempts in life, predicted age at onset of senescence and senescence rates. ####

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

#### Link between the number of breeding attempts in life, predicted age at onset of senescence and senescence rate ####


##### subset of dead individuals ##### 


SubsetDead_cohortBefore2000 = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$Year_last_seen <= 2020  &
                                                                                 df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort <2000),]

table(SubsetDead_cohortBefore2000$population)


SubsetDead_cohortBefore2000_noOutlier  =SubsetDead_cohortBefore2000[which( SubsetDead_cohortBefore2000$OnsetSenesc_from_backTransformedCoefs<40),]

cor(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen")],
    SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen")])
cor(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island")],
    SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island")])


SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs = -SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs

#### Figure S4.1 #### 
library(scales)

Fig_S4_1_A = ggplot(data = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island"),], 
                 aes(x =SenesceRate_from_backTransformedCoefs , y=OnsetSenesc_from_backTransformedCoefs , col = nb_BA_life))+
  geom_point(size = 1.5, alpha = 0.6)+
  scale_colour_viridis_c(begin = 0.92, end = 0.05, option = "inferno", trans = pseudo_log_trans(base = 10), limits = c(0, 30))+
  theme_classic()+ xlim(c(0.01,0.06))+ ylim(c(17,32))+ 
  xlab("Senescence rate") + ylab("Age at onset of senescence")+
  labs(col = "n breeding attempts")+
  theme( legend.text = element_text(size=8) , legend.key.size = unit(0.55, "cm"), legend.title = element_text(size=9))
         
# theme( legend.text = element_text(size=9) , legend.title = element_text(size=11),legend.position = c(0.95, 0.5) )
Fig_S4_1_A
 

Fig_S4_1_B = ggplot(data = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen"  ),], 
                 aes(x =SenesceRate_from_backTransformedCoefs , y=OnsetSenesc_from_backTransformedCoefs , col = nb_BA_life))+
  geom_point(size = 1.5, alpha = 0.6)+
  scale_colour_viridis_c(begin = 0.90, end = 0.02, option = "mako", trans = pseudo_log_trans(base = 10))+
  theme_classic()+ xlim(c(0.01,0.06))+ ylim(c(17,32))+
  xlab("Senescence rate") + ylab("Age at onset of senescence")+
  labs(col = "n breeding attempts")+
  theme( legend.text = element_text(size=8) , legend.key.size = unit(0.55, "cm"), legend.title = element_text(size=9)  )
# theme( legend.text = element_text(size=9), legend.title = element_text(size=11),legend.position = c(0.99, 0.5) )
# theme( legend.text = element_text(size=9), legend.title = element_text(size=11) )
Fig_S4_1_B


cowplot::plot_grid(Fig_S4_1_A, Fig_S4_1_B)


# std variables
SubsetDead_cohortBefore2000_noOutlier$OnsetSenescence_std = robustHD::standardize(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs)
SubsetDead_cohortBefore2000_noOutlier$SenesceRate_std = robustHD::standardize(SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs)


prior_nBA <- c(
  prior(normal(0, 5), class = "Intercept"),
  prior(normal(0, 2), class = "b")
  
)

GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence  = brm(nb_BA_life   ~   OnsetSenescence_std * population  ,
                                                 data=SubsetDead_cohortBefore2000_noOutlier ,
                                                 family = poisson   ,
                                                 warmup = 2000,
                                                 iter = 4000,
                                                 chains =8,
                                                 cores = 8,
                                                 prior = prior_nBA,
                                                 backend = "cmdstanr" ,
                                                 control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence) 
bayes_R2(GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence) 

# pp_check(GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence)

SubsetDead_cohortBefore2000_noOutlier$population=as.factor(SubsetDead_cohortBefore2000_noOutlier$population)
SubsetDead_cohortBefore2000_noOutlier$population_refKer = relevel(SubsetDead_cohortBefore2000_noOutlier$population , ref=2)

GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence_refKer  = brm(nb_BA_life   ~   OnsetSenescence_std * population_refKer  ,
                                                        data=SubsetDead_cohortBefore2000_noOutlier ,
                                                        family = poisson   ,
                                                        warmup = 2000,
                                                        iter = 4000,
                                                        chains =8,
                                                        cores = 8,
                                                        prior = prior_nBA,
                                                        backend = "cmdstanr" ,
                                                        control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence_refKer) 

 

GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate  = brm(nb_BA_life   ~   SenesceRate_std * population  ,
                                             data=SubsetDead_cohortBefore2000_noOutlier ,
                                             family = poisson   ,
                                             warmup = 2000,
                                             iter = 4000,
                                             chains =8,
                                             cores = 8,
                                             prior = prior_nBA,
                                             backend = "cmdstanr" ,
                                             control = list(adapt_delta = 0.85)
)

summary(GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate)  
# pp_check(GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate)
bayes_R2(GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate)
 

GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate_refKer  = brm(nb_BA_life   ~   SenesceRate_std * population_refKer  ,
                                                    data=SubsetDead_cohortBefore2000_noOutlier ,
                                                    family = poisson   ,
                                                    warmup = 2000,
                                                    iter = 4000,
                                                    chains =8,
                                                    cores = 8,
                                                    prior = prior_nBA,
                                                    backend = "cmdstanr" ,
                                                    control = list(adapt_delta = 0.85)
)

summary(GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate_refKer) 



onset_mean <- mean(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs)
onset_sd   <- sd(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs)

rate_mean  <- mean(SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs)
rate_sd    <- sd(SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs)

ranges_onset <- SubsetDead_cohortBefore2000_noOutlier %>%
  group_by(population) %>%
  summarise(
    xmin = min(OnsetSenescence_std),
    xmax = max(OnsetSenescence_std)
  )

newdata_onset <- ranges_onset %>%
  rowwise() %>%
  do(data.frame(
    population = .$population,
    OnsetSenescence_std = seq(.$xmin, .$xmax, length.out = 100)
  ))


newdata_onset$OnsetSenesc_plot <-
  newdata_onset$OnsetSenescence_std * onset_sd + onset_mean

preds_onset <- fitted(
  GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence,
  newdata = newdata_onset,
  summary = TRUE
)

newdata_onset$Estimate <- preds_onset[, "Estimate"]
newdata_onset$Q2.5     <- preds_onset[, "Q2.5"]
newdata_onset$Q97.5    <- preds_onset[, "Q97.5"]

Fig_S4_1_C <- ggplot() +
  geom_ribbon(
    data = newdata_onset,
    aes(x = OnsetSenesc_plot, ymin = Q2.5, ymax = Q97.5, fill = population),
    alpha = 0.3
  ) +
  geom_line(
    data = newdata_onset,
    aes(x = OnsetSenesc_plot, y = Estimate, color = population),
    linewidth = 1
  ) +
  geom_jitter(
    data = SubsetDead_cohortBefore2000_noOutlier,
    aes(
      x = OnsetSenesc_from_backTransformedCoefs,
      y = nb_BA_life,
      color = population
    ),
    width = 0,
    height = 0.2,
    alpha = 0.5,
    size = 1.5
  ) +
  scale_fill_manual(values = c("darkorange", "turquoise4")) +
  scale_color_manual(values = c("darkorange", "turquoise4")) +
  theme_classic() +
  xlab("Age at onset of senescence") +
  ylab("Number of breeding attempts\nin life") +
  theme(legend.position = "none")

Fig_S4_1_C

 

ranges_rate <- SubsetDead_cohortBefore2000_noOutlier %>%
  group_by(population) %>%
  summarise(
    xmin = min(SenesceRate_std),
    xmax = max(SenesceRate_std)
  )

newdata_rate <- ranges_rate %>%
  rowwise() %>%
  do(data.frame(
    population = .$population,
    SenesceRate_std = seq(.$xmin, .$xmax, length.out = 100)
  ))
newdata_rate$SenesceRate_plot <-
  newdata_rate$SenesceRate_std * rate_sd + rate_mean

preds_rate <- fitted(
  GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate,
  newdata = newdata_rate,
  summary = TRUE
)

newdata_rate$Estimate <- preds_rate[, "Estimate"]
newdata_rate$Q2.5     <- preds_rate[, "Q2.5"]
newdata_rate$Q97.5    <- preds_rate[, "Q97.5"]

Fig_S4_1_D <- ggplot() +
  geom_ribbon(
    data = newdata_rate,
    aes(x = SenesceRate_plot, ymin = Q2.5, ymax = Q97.5, fill = population),
    alpha = 0.3
  ) +
  geom_line(
    data = newdata_rate,
    aes(x = SenesceRate_plot, y = Estimate, color = population),
    linewidth = 1
  ) +
  geom_jitter(
    data = SubsetDead_cohortBefore2000_noOutlier,
    aes(
      x = SenesceRate_from_backTransformedCoefs,
      y = nb_BA_life,
      color = population
    ),
    width = 0,
    height = 0.2,
    alpha = 0.5,
    size = 1.5
  ) +
  scale_fill_manual(values = c("darkorange", "turquoise4")) +
  scale_color_manual(values = c("darkorange", "turquoise4")) +
  theme_classic() +
  xlab("Senescence rate") +
  ylab("Number of breeding attempts\nin life") +
  theme(legend.position = "none")

Fig_S4_1_D
cowplot::plot_grid(Fig_S4_1_A,
                   Fig_S4_1_B, 
                   Fig_S4_1_C,
                   Fig_S4_1_D, 
                   align="hv",
                   ncol = 2, labels = c("A", "B", "C", "D"))


# ggsave(paste0(dir_data, "SupportingInformation_S4_Figure_S4_1.png"))



#####  Dead and alive ##### 
#### Figure S4.2 ####  

df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot = "Alive"
df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$Year_last_seen <= 2020  )]= "Dead"

df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs = -df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs

table(df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot , df_AgeStd_randomSlopesAndInterceptsIDlevel$population)

df_AgeStd_randomSlopesAndInterceptsIDlevel$pop_DeadOrNot = paste0(df_AgeStd_randomSlopesAndInterceptsIDlevel$population,
                                                                  "_",
                                                                  df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot)

 


df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs > 0 &
                                                                                                          df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),] 

df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenescence_std = robustHD::standardize(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenesc_from_backTransformedCoefs)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_std = robustHD::standardize(-df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_from_backTransformedCoefs)

table(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot, df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population)
 


prior_nBA <- c(
  prior(normal(0, 5), class = "Intercept"),
  prior(normal(0, 2), class = "b")
  
) 

GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence  = brm(nb_BA_life   ~   OnsetSenescence_std * population * DeadOrNot ,
                                                        data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                        family = poisson   ,
                                                        warmup = 2000,
                                                        iter = 4000,
                                                        chains =8,
                                                        cores = 8,
                                                        prior = prior_nBA,
                                                        backend = "cmdstanr" ,
                                                        control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence) 


df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population=as.factor(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population_refKer = relevel(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population , ref=2)

GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence_refKer  = brm(nb_BA_life   ~    OnsetSenescence_std * population_refKer * DeadOrNot ,
                                                               data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                               family = poisson   ,
                                                               warmup = 2000,
                                                               iter = 4000,
                                                               chains =8,
                                                               cores = 8,
                                                               prior = prior_nBA,
                                                               backend = "cmdstanr" ,
                                                               control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence_refKer) 

bayes_R2(GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence_refKer)


df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot=as.factor(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot_refDead = relevel(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot , ref=2)


GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence_refDead  = brm(nb_BA_life   ~   OnsetSenescence_std * population * DeadOrNot_refDead ,
                                                        data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                        family = poisson   ,
                                                        warmup = 2000,
                                                        iter = 4000,
                                                        chains =8,
                                                        cores = 8,
                                                        prior = prior_nBA,
                                                        backend = "cmdstanr" ,
                                                        control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence_refDead) 



GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence_refDead_refKer  = brm(nb_BA_life   ~    OnsetSenescence_std * population_refKer * DeadOrNot_refDead ,
                                                               data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                               family = poisson   ,
                                                               warmup = 2000,
                                                               iter = 4000,
                                                               chains =8,
                                                               cores = 8,
                                                               prior = prior_nBA,
                                                               backend = "cmdstanr" ,
                                                               control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence_refDead_refKer) 


# Senescence rate

GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate  = brm(nb_BA_life   ~   SenesceRate_std * population * DeadOrNot ,
                                                        data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                        family = poisson   ,
                                                    warmup = 2000,
                                                    iter = 4000,
                                                    chains =8,
                                                    cores = 8,
                                                    prior = prior_nBA,
                                                        backend = "cmdstanr" ,
                                                        control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate) 
bayes_R2(GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate) 

df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population=as.factor(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population_refKer = relevel(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population , ref=2)

GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate_refKer  = brm(nb_BA_life   ~    SenesceRate_std * population_refKer * DeadOrNot ,
                                                               data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                               family = poisson   ,
                                                           warmup = 2000,
                                                           iter = 4000,
                                                           chains =8,
                                                           cores = 8,
                                                           prior = prior_nBA,
                                                               backend = "cmdstanr" ,
                                                               control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate_refKer) 

df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot=as.factor(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot_refDead = relevel(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot , ref=2)


GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate_refDead  = brm(nb_BA_life   ~   SenesceRate_std * population * DeadOrNot_refDead ,
                                                        data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                        family = poisson   ,
                                                        warmup = 2000,
                                                        iter = 4000,
                                                        chains =8,
                                                        cores = 8,
                                                        prior = prior_nBA,
                                                        backend = "cmdstanr" ,
                                                        control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate_refDead) 



GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate_refKer_refDead  = brm(nb_BA_life   ~    SenesceRate_std * population_refKer * DeadOrNot_refDead ,
                                                               data=df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier ,
                                                               family = poisson   ,
                                                               warmup = 2000,
                                                               iter = 4000,
                                                               chains =8,
                                                               cores = 8,
                                                               prior = prior_nBA,
                                                               backend = "cmdstanr" ,
                                                               control = list(adapt_delta = 0.85)
)
summary(GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate_refKer_refDead) 

# pp_check(GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate_refKer_refDead)
 

ranges_onset_v2 <- df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier %>%
  group_by(population, DeadOrNot, pop_DeadOrNot) %>%
  summarise(    xmin = min(OnsetSenescence_std),    xmax = max(OnsetSenescence_std),
                .groups = "drop"
  )

newdata_onset_v2 <- ranges_onset_v2 %>%
  rowwise() %>%
  do(data.frame(    population = .$population,    DeadOrNot  = .$DeadOrNot,    pop_DeadOrNot = .$pop_DeadOrNot,
                    OnsetSenescence_std = seq(.$xmin, .$xmax, length.out = 100)
  )) 
GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence
preds_onset_v2 <- fitted(  GLM_KerBI_Poisson_BRM_nb_BA_life_OnsetSenescence,
                           newdata = newdata_onset_v2,  summary = TRUE
)

newdata_onset_v2$Estimate <- preds_onset_v2[, "Estimate"]
newdata_onset_v2$Q2.5     <- preds_onset_v2[, "Q2.5"]
newdata_onset_v2$Q97.5    <- preds_onset_v2[, "Q97.5"]

newdata_onset_v2$OnsetSenesc_plot <-
  newdata_onset_v2$OnsetSenescence_std *
  sd(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenesc_from_backTransformedCoefs) +
  mean(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenesc_from_backTransformedCoefs)

Fig_S4_2_C <- ggplot() +
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
  xlab("Age at onset of senescence") +  ylab("Number of breeding attempts\nin life") +
  theme(legend.position = "none")

Fig_S4_2_C


ranges_rate_v2 <- df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier %>%
  group_by(population, DeadOrNot, pop_DeadOrNot) %>%
  summarise(    xmin = min(SenesceRate_std),    xmax = max(SenesceRate_std),    .groups = "drop"  )
newdata_rate_v2 <- ranges_rate_v2 %>%
  rowwise() %>%
  do(data.frame(    population = .$population,    DeadOrNot  = .$DeadOrNot,    pop_DeadOrNot = .$pop_DeadOrNot,
                    SenesceRate_std = seq(.$xmin, .$xmax, length.out = 100)
  ))

preds_rate_v2 <- fitted(  GLM_KerBI_Poisson_BRM_nb_BA_life_SenesceRate,
                          newdata = newdata_rate_v2,
                          summary = TRUE
)

newdata_rate_v2$Estimate <- preds_rate_v2[, "Estimate"]
newdata_rate_v2$Q2.5     <- preds_rate_v2[, "Q2.5"]
newdata_rate_v2$Q97.5    <- preds_rate_v2[, "Q97.5"]

newdata_rate_v2$SenesceRate_plot <-   - newdata_rate_v2$SenesceRate_std *
  sd(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_from_backTransformedCoefs) +
  mean(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_from_backTransformedCoefs)


Fig_S4_2_D <- ggplot() +
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
  ylab("Number of breeding attempts\nin life") +
  theme(legend.position = "none")

Fig_S4_2_D


cowplot::plot_grid( Fig_S4_2_C,
                    Fig_S4_2_D, 
                    align="hv",
                    ncol = 2, labels = c("A", "B"),  label_size = 12)


# ggsave(paste0(dir_data, "SupportingInformation_S4_Figure_S4_2.png"))
 
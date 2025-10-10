
# Load packages
library(ggplot2)
library(ggdist)


# set the location of the data
dir_data = "C:/Users/bmohring//"


# Open dataset
 

df_AgeStd_randomSlopesAndInterceptsIDlevel = read.csv( paste0(dir_data, "data_BBAL_senescence_estimates_priors.csv"))
df_plotRes_Senescence = read.csv( paste0(dir_data, "data_BBAL_senescence_posteriorEstimates_priors.csv"))

data_BBAL = read.csv(  paste0(dir_data, "data_BBAL.csv"))
data_BBAL=data_BBAL[,-1]

# Load model output
load(paste0(dir_data, "model_output_GLMM_senescence_priors.RData"))


# calculate age mean and sd
mean_age_forStd =  mean(data_BBAL$Age) 
sd_age_forStd =  sd(data_BBAL$Age) 


#### Figure 2 #### 


# Desired raw age breaks
raw_breaks <- c(5, 10, 15, 20, 25, 30, 35, 40, 45,50)
# Convert raw breaks to standardized scale
std_breaks <- (raw_breaks - mean_age_forStd) / sd_age_forStd

pA  = plot(conditional_effects(GLMM_brms_senescence), ask = FALSE)
pA=pA$`Age_std:population`+  
  theme_classic() + 
  scale_fill_manual(values = c("darkorange", "turquoise4")) + 
  scale_color_manual(values = c("darkorange", "turquoise4")) +
  xlab("Age") + ylab("Reproductive success")+ 
  theme(legend.position = "none"  ,
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 10))+
  ylim(c(0,1))+
  scale_x_continuous(
    breaks =std_breaks,
    labels =raw_breaks
  )



p_onsetSenescence = ggplot(data = df_plotRes_Senescence , 
                           aes(x = OnsetSenesc_from_backTransformedCoefs ,
                               # y =population,
                               group  = population , fill  = population ))+
  stat_halfeye( side = "top",alpha = 0.6 ,
                # size = 3,
                # .width = c(0.66, 0.95),
                # position = "dodge",
                # position = position_dodge(width = 0.1),
                # position = "dodge", 
                height = 0.9)+
  xlim(16, 30)+
  scale_fill_manual(values = c( "darkorange","turquoise4" ))+
  # geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, col = "grey40") +
  ylab("Frequency")+
  xlab("Age at onset of senescence")+
  theme_classic()+ theme(legend.position = "none",
                         axis.title.x = element_text(size = 10),  
                         axis.title.y = element_text(size = 10))
# +
# ylab("Predicted age at onset of senescence")
p_onsetSenescence



p_SenesceRate = ggplot(data = df_plotRes_Senescence , 
                       aes(x = SenesceRate_from_backTransformedCoefs ,
                           # y =Variable, 
                           group  = population , fill  = population ))+
  stat_halfeye( side = "top",alpha = 0.6 ,
                # size = 3,
                # .width = c(0.66, 0.95),
                # position = "dodge",
                # position = position_dodge(width = 0.4),
                # position = "dodge", 
                height = 0.9)+
  # xlim(-1.5, 1.5)+
  scale_fill_manual(values = c( "darkorange","turquoise4" ))+
  # geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, col = "grey40") +
  ylab("Frequency")+
  xlab("Senescence rate")+
  theme_classic()+ theme(legend.position = "none",
                         axis.title.x = element_text(size = 10),  
                         axis.title.y = element_text(size = 10))
# +
# ylab("Predicted age at onset of senescence")
p_SenesceRate



p_PredictedRS_atpeak = ggplot(data = df_plotRes_Senescence , 
                              aes(x = PredictedRS_atpeak ,
                                  # y =Variable, 
                                  group  = population , fill  = population ))+
  stat_halfeye( side = "top",alpha = 0.6 , 
                height = 0.9)+ 
  scale_fill_manual(values = c( "darkorange","turquoise4" ))+ 
  ylab("Frequency")+
  xlab("Probability of successful reproduction\nat the onset of senescence")+
  theme_classic()+ theme(legend.position = "none",
                         axis.title.x = element_text(size = 10),  
                         axis.title.y = element_text(size = 10))
# +
# ylab("Predicted age at onset of senescence")
p_PredictedRS_atpeak


p_PredictedRS_at10yo = ggplot(data = df_plotRes_Senescence , 
                              aes(x = PredictedRS_at10yo ,
                                  # y =Variable, 
                                  group  = population , fill  = population ))+
  stat_halfeye( side = "top",alpha = 0.6 , 
                height = 0.9)+ 
  scale_fill_manual(values = c( "darkorange","turquoise4" ))+ 
  ylab("Frequency")+
  xlab("Early-life probability\nof successful reproduction")+
  theme_classic()+ theme(legend.position = "none",
                         axis.title.x = element_text(size = 10),  
                         axis.title.y = element_text(size = 10) )
p_PredictedRS_at10yo



pB  = cowplot::plot_grid(p_onsetSenescence,
                         p_SenesceRate,
                         p_PredictedRS_at10yo,
                         p_PredictedRS_atpeak,
                         ncol=2,
                         labels = c("B", "C", "D", "E")
)


cowplot::plot_grid(pA,
                   pB, 
                   ncol = 1,
                   labels = c("A"))


ggsave(paste0(dir_data, "Figure2_prior.png"))






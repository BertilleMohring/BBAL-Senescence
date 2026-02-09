#### Code used to produce Figure 3 and calculate the correlation between onset and rate of senescence ####

# Load packages
library(ggplot2)
library(ggdist)

# set the location of the data
dir_data_BBAL = "C:/Users/bmohring/Documents/GitHub/BBAL-Senescence/Data/"
dir_data = "C:/Users/bmohring/Documents/GitHub/BBAL-Senescence/"

# Open datasets

df_AgeStd_randomSlopesAndInterceptsIDlevel = read.csv( paste0(dir_data, "data_BBAL_senescence_estimates_priors.csv"))
# df_plotRes_Senescence = read.csv(paste0(dir_data, "data_BBAL_senescence_posteriorEstimates_priors.csv"))

data_BBAL = read.csv(paste0(dir_data_BBAL, "data_BBAL.csv"))
data_BBAL=data_BBAL[,-1]

# Load model output
# load(paste0(dir_data, "model_output_GLMM_senescence_priors.RData"))


# calculate age mean and sd
mean_age_forStd =  mean(data_BBAL$Age) 
sd_age_forStd =  sd(data_BBAL$Age) 

df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs = -df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs



#### 3. Figure 3: Relationship between individual predicted age at the onset of senescence and senescence rate  #### 

#### __ Outlier removed  #### 

p = ggplot(df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40  ),],
           aes(x = SenesceRate_from_backTransformedCoefs, y =OnsetSenesc_from_backTransformedCoefs , 
               col = population, fill = population, group = population ))+
  geom_point(size = 1.4, alpha = 0.3)+
  scale_fill_manual(values = c("darkorange", "turquoise4")) + 
  scale_color_manual(values = c("darkorange", "turquoise4")) + 
  # geom_smooth(method= lm)+
  # geom_smooth(inherit.aes = FALSE,
  #             data  =df_rawAge_randomSlopesAndInterceptsIDlevel,
  #             aes(x = SenesceRate, y =OnsetSenesc_backTransformed),
  #             method= glm, method.args = list(family = "Gamma"))+
  theme_classic()+ theme(legend.position="none")+
  xlab("Senescence rate") + ylab("Age at onset of senescence")

pmarginal = ggExtra::ggMarginal(p, type = "density",  groupFill  = TRUE, groupColour = TRUE)

pmarginal

ggsave(paste0(dir_data, "Figure3.png"),
       plot = pmarginal)



####  correlation between onset of senescence and senescence rate #### 

#### __ Outlier removed  #### 

# Kerguelen
cor(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Kerguelen" &
                                                                                             df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40 )],
    df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Kerguelen" &
                                                                                             df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40 )])
# Bird Island
cor(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Bird island")],
    df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Bird island")])

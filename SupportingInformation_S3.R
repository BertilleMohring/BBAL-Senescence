#### CODE USED TO PRODUCE SUPPORTING INFORMATION S3 ####
#### Supporting Information S3: Relationship between predicted age at onset of senescence and predicted senescence rate  ####
#### when including potential outliers. ####

# Load packages
library(ggplot2)
library(ggdist)

# set the location of the data
dir_data_BBAL = "C:/Users/bmohring/Documents/GitHub/BBAL-Senescence/Data/"
dir_data = "C:/Users/bmohring/Documents/GitHub/BBAL-Senescence/"


# Open dataset
df_AgeStd_randomSlopesAndInterceptsIDlevel = read.csv( paste0(dir_data, "data_BBAL_senescence_estimates_priors.csv"))
df_plotRes_Senescence = read.csv(paste0(dir_data, "data_BBAL_senescence_posteriorEstimates_priors.csv"))

data_BBAL = read.csv( paste0(dir_data_BBAL, "data_BBAL.csv"))
data_BBAL=data_BBAL[,-1]

# Load model output
# load(paste0(dir_data, "model_output_GLMM_senescence_priors.RData"))


# calculate age mean and sd
mean_age_forStd =  mean(data_BBAL$Age) 
sd_age_forStd =  sd(data_BBAL$Age) 

df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs = - df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs


#### Correlation between onset and rate of senescence #### 

# Kerguelen
cor(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Kerguelen")],
    df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Kerguelen")])

# Bird Island
cor(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Bird island")],
    df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Bird island")])



#### Figure S3.1: Relationship between individual predicted age at the onset of senescence and senescence rate  #### 


p = ggplot(df_AgeStd_randomSlopesAndInterceptsIDlevel,
           aes(x = SenesceRate_from_backTransformedCoefs, y =OnsetSenesc_from_backTransformedCoefs , 
               col = population, fill = population, group = population ))+
  geom_point(size = 1.4, alpha = 0.3)+
  scale_fill_manual(values = c("darkorange", "turquoise4")) + 
  scale_color_manual(values = c("darkorange", "turquoise4")) +  
  theme_classic()+ theme(legend.position="none")+
  xlab("Senescence rate") + ylab("Age at onset of senescence")

pmarginal = ggExtra::ggMarginal(p, type = "density",  groupFill  = TRUE, groupColour = TRUE)

pmarginal

# ggsave(paste0(dir_data, "SupportingInformation_S3_Figure_S3_1.png"),
#        plot = pmarginal)

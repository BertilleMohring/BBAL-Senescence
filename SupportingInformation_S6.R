
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
library(corrplot)


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


SubsetDead_cohortBefore2000_BI = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island"),
                                                                       which(colnames(SubsetDead_cohortBefore2000_noOutlier) %in% c("OnsetSenesc_from_backTransformedCoefs" ,
                                                                                                                                    "SenesceRate_from_backTransformedCoefs", 
                                                                                                                                    "PredictedRS_at10yo"   , "PredictedRS_atpeak" ))]

SubsetDead_cohortBefore2000_Ker = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen"),
                                                                       which(colnames(SubsetDead_cohortBefore2000_noOutlier) %in% c("OnsetSenesc_from_backTransformedCoefs" ,
                                                                                                                                    "SenesceRate_from_backTransformedCoefs", 
                                                                                                                                    "PredictedRS_at10yo"   , "PredictedRS_atpeak" ))]



SubsetDead_cohortBefore2000_Ker$SenesceRate_from_backTransformedCoefs = - SubsetDead_cohortBefore2000_Ker$SenesceRate_from_backTransformedCoefs
SubsetDead_cohortBefore2000_BI$SenesceRate_from_backTransformedCoefs = - SubsetDead_cohortBefore2000_BI$SenesceRate_from_backTransformedCoefs

#### Correlation matrix ####


res_BI <- cor(SubsetDead_cohortBefore2000_BI)
 
colnames(res_BI)
colnames(res_BI)= c(  "RS_onset", "RS_10yo", "Onset", "Rate")
rownames(res_BI)
rownames(res_BI)= c(  "RS_onset", "RS_10yo", "Onset", "Rate")
 
res_Ker <- cor(SubsetDead_cohortBefore2000_Ker)

colnames(res_Ker)
colnames(res_Ker)= c(  "RS_onset", "RS_10yo", "Onset", "Rate")
rownames(res_Ker)
rownames(res_Ker)= c(  "RS_onset", "RS_10yo", "Onset", "Rate")

 
 


par(mfrow = c(1, 2))   

#### __ Bird island  #### 
corrplot(res_BI,  type = "upper", 
         # order = "hclust", 
         method="number",
         tl.col = "black", 
         tl.srt = 45, 
         addCoef.col = "white",      
         number.cex = 1)


#### __ Kerguelen  #### 
corrplot(res_Ker,  
         type = "upper", 
         # order = "hclust",
         method="number",
         tl.col = "black", 
         tl.srt = 45,  
         number.cex = 1)


par(mfrow = c(1, 1)) 
 
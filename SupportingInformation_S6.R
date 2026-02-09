#### CODE USED TO PRODUCE SUPPORTING INFORMATION S6 ####
#### Supporting Information S6: Propagation of uncertainty in the estimation of age at the onset of senescence and senescence rate ####


# Load packages
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
 
# Open datasets
data_BBAL = read.csv( paste0(dir_data_BBAL, "data_BBAL.csv"))
data_BBAL=data_BBAL[,-1]

df_AgeStd_randomSlopesAndInterceptsIDlevel = read.csv( paste0(dir_data, "data_BBAL_senescence_estimates_priors.csv"))
df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort = substr(df_AgeStd_randomSlopesAndInterceptsIDlevel$cohortPop,1,4)
df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort =as.numeric(df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort )

df_plotRes_Senescence = read.csv( paste0(dir_data, "data_BBAL_senescence_posteriorEstimates_priors.csv"))

# Load model output
load( paste0(dir_data, "model_output_GLMM_senescence_priors.RData"))

# calculate age mean and sd
mean_age_forStd =  mean(data_BBAL$Age) 
sd_age_forStd =  sd(data_BBAL$Age) 


# load functions
deriv_logit_quad<-function(xval,quad,lin,const)
{
  output<-(2*quad*xval+lin)*exp(quad*xval^2+lin*xval+const)/((1+exp(quad*xval^2+lin*xval+const))^2)
  return(output)
}
back_trans_intercept<-function(b0, b1, b2, mean, sd){
  output<- (b0 - (b1 * mean)/sd + (b2*I(mean)^2)/I(sd)^2 ) 
  return(output)
}
back_trans_slope<-function(b0, b1, b2, mean, sd){
  output<- (b1/sd )  - ((2 * b2 * mean ) /I(sd)^2 )
  return(output)
}
back_trans_slope2<-function(b0, b1, b2, mean, sd){
  output<-(  b2 / I(sd)^2 ) 
  return(output)
}

back_trans_fun<-function(x, b0, b1, b2, mean, sd){
  b0_new =  (b0 - (b1 * mean)/sd + (b2*I(mean)^2)/I(sd)^2 ) 
  b1_new = (b1/sd )  - ((2 * b2 * mean ) /I(sd)^2 )
  b2_new =   (b2 /  I(sd)^2 ) 
  output = invlogit( b0_new + b1_new * x + b2_new * I(x)^2    )
} 


#### 6. Randomly sample from posterior draws ####

# Extract values for  95% credible intervals of age at the onset of senescence and senescence rate for each population 
# Kerguelen
lower_CI_Ker_onset = quantile(df_plotRes_Senescence$OnsetSenesc_from_backTransformedCoefs[which(df_plotRes_Senescence$population=="Kerguelen")], 0.025)
lower_CI_BI_onset  = quantile(df_plotRes_Senescence$OnsetSenesc_from_backTransformedCoefs[which(df_plotRes_Senescence$population=="Bird island")], 0.025)
upper_CI_Ker_onset = quantile(df_plotRes_Senescence$OnsetSenesc_from_backTransformedCoefs[which(df_plotRes_Senescence$population=="Kerguelen")], 0.975)
upper_CI_BI_onset  = quantile(df_plotRes_Senescence$OnsetSenesc_from_backTransformedCoefs[which(df_plotRes_Senescence$population=="Bird island")], 0.975)
# Bird Island
lower_CI_Ker_rate = quantile(df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs[which(df_plotRes_Senescence$population=="Kerguelen")], 0.025)
lower_CI_BI_rate = quantile(df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs[which(df_plotRes_Senescence$population=="Bird island")], 0.025)
upper_CI_Ker_rate = quantile(df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs[which(df_plotRes_Senescence$population=="Kerguelen")], 0.975)
upper_CI_BI_rate = quantile(df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs[which(df_plotRes_Senescence$population=="Bird island")], 0.975)

# add draw number
df_plotRes_Senescence$draw = rep(c(1:16000), each =2)

# subset datasets within 95% credible intervals of age at the onset of senescence and senescence rate
# Kerguelen
list_draws_to_pick_Ker = df_plotRes_Senescence$draw[which(df_plotRes_Senescence$population=="Kerguelen" &
                                                            df_plotRes_Senescence$OnsetSenesc_from_backTransformedCoefs > lower_CI_Ker_onset&
                                                            df_plotRes_Senescence$OnsetSenesc_from_backTransformedCoefs < upper_CI_Ker_onset &
                                                            df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs > lower_CI_Ker_rate &
                                                            df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs < upper_CI_Ker_rate  )]
length(list_draws_to_pick_Ker)

# Bird Island
list_draws_to_pick_BI = df_plotRes_Senescence$draw[which(df_plotRes_Senescence$population=="Bird island" &
                                                           df_plotRes_Senescence$OnsetSenesc_from_backTransformedCoefs > lower_CI_BI_onset&
                                                           df_plotRes_Senescence$OnsetSenesc_from_backTransformedCoefs < upper_CI_BI_onset &
                                                           df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs > lower_CI_BI_rate &
                                                           df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs < upper_CI_BI_rate  )]
length(list_draws_to_pick_BI)

# combine the data
list_draws_to_pick_Ker_BI = list_draws_to_pick_Ker[which(list_draws_to_pick_Ker %in% list_draws_to_pick_BI)]
length(list_draws_to_pick_Ker_BI)

# Bird island

# sample draws
posterior_samples_finalModel <- posterior_samples(GLMM_brms_senescence)

# randomly sample 2000 draws
list_i = sample(list_draws_to_pick_Ker_BI, 2000)
list_i

df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort=substr(df_AgeStd_randomSlopesAndInterceptsIDlevel$cohortPop, 1,4)

# make a dataset using the first draw
df_draws = data.frame(ID =  df_AgeStd_randomSlopesAndInterceptsIDlevel$ID,
                      population = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                      cohortPop = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                      draw = rep(list_i[1], length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                      Intercept = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                      Slope = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                      Slope2 = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                      Intercept_randID = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                      Slope_randID = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                      Slope2_randID = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                      Intercept_cohortPop = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)) )

for (k in 1:length(df_draws$ID)){
  RING  = df_draws$ID[k]
  DRAW = df_draws$draw[k]
  df_draws$population[which(df_draws$ID == RING)] = as.character(df_AgeStd_randomSlopesAndInterceptsIDlevel$population[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID ==RING)])
  df_draws$cohortPop[which(df_draws$ID == RING)] = df_AgeStd_randomSlopesAndInterceptsIDlevel$cohortPop[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID ==RING)]
  
  POPULATION = df_AgeStd_randomSlopesAndInterceptsIDlevel$population[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID ==RING)]
  cohortPop =  df_AgeStd_randomSlopesAndInterceptsIDlevel$cohortPop[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID ==RING)]
  COHORT =  df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID ==RING)]
  
  if(POPULATION == "Bird island"){
    # cohortPop_POP = paste0(cohortPop, "-Bird island")
    
    df_draws$Intercept[which(df_draws$ID == RING)] = posterior_samples_finalModel$b_Intercept[DRAW]
    df_draws$Slope[which(df_draws$ID == RING)] = posterior_samples_finalModel$b_Age_std[DRAW]
    df_draws$Slope2[which(df_draws$ID == RING)] = posterior_samples_finalModel$b_IAge_stdE2[DRAW]
    
    df_draws$Intercept_cohortPop[which(df_draws$ID == RING)] = posterior_samples_finalModel[ DRAW, which(colnames(posterior_samples_finalModel) == paste0("r_cohortPop[", COHORT , "-Bird.island", ",Intercept]" ))]
    
  }
  if(POPULATION == "Kerguelen"){
    # cohortPop_POP = paste0(cohortPop, "-Kerguelen")
    
    df_draws$Intercept[which(df_draws$ID == RING)] =  posterior_samples_finalModel$b_Intercept[DRAW]  + posterior_samples_finalModel$b_populationKerguelen[DRAW]
    df_draws$Slope[which(df_draws$ID == RING)] = posterior_samples_finalModel$b_Age_std[DRAW] + posterior_samples_finalModel$`b_Age_std:populationKerguelen`[DRAW]
    df_draws$Slope2[which(df_draws$ID == RING)] = posterior_samples_finalModel$b_IAge_stdE2[DRAW] + posterior_samples_finalModel$`b_IAge_stdE2:populationKerguelen`[DRAW]
    
    df_draws$Intercept_cohortPop[which(df_draws$ID == RING)] = posterior_samples_finalModel[ DRAW, which(colnames(posterior_samples_finalModel) == paste0("r_cohortPop[", COHORT , "-Kerguelen", ",Intercept]" ))]
    
  }
  
  
  df_draws$Intercept_randID[which(df_draws$ID == RING)] = posterior_samples_finalModel[ DRAW, which(colnames(posterior_samples_finalModel)== paste0("r_ID[", RING, ",Intercept]" ))]
  df_draws$Slope_randID[which(df_draws$ID == RING)] = posterior_samples_finalModel[ DRAW, which(colnames(posterior_samples_finalModel)== paste0("r_ID[", RING, ",Age_std]" ))]
  df_draws$Slope2_randID[which(df_draws$ID == RING)] = posterior_samples_finalModel[ DRAW, which(colnames(posterior_samples_finalModel)== paste0("r_ID[", RING, ",IAge_stdE2]" ))]
  
  
}


# repeat and bind dataset for the other draws

for (i in 2:length(list_i)){
  # dataset to fill in and bind
  df_draws_toBind = data.frame(ID =  df_AgeStd_randomSlopesAndInterceptsIDlevel$ID,
                               population = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                               cohortPop = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                               draw = rep(list_i[i], length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                               Intercept = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                               Slope = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                               Slope2 = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                               Intercept_randID = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                               Slope_randID = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                               Slope2_randID = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)),
                               Intercept_cohortPop = rep(NA, length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)) )
  
  # fill in for each ID
  for (k in 1:length(df_draws_toBind$ID)){
    RING  = df_draws_toBind$ID[k]
    DRAW = df_draws_toBind$draw[k]
    df_draws_toBind$population[which(df_draws_toBind$ID == RING)] = as.character(df_AgeStd_randomSlopesAndInterceptsIDlevel$population[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID ==RING)])
    df_draws_toBind$cohortPop[which(df_draws_toBind$ID == RING)] = df_AgeStd_randomSlopesAndInterceptsIDlevel$cohortPop[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID ==RING)]
    
    POPULATION = df_AgeStd_randomSlopesAndInterceptsIDlevel$population[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID ==RING)]
    cohortPop =  df_AgeStd_randomSlopesAndInterceptsIDlevel$cohortPop[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID ==RING)]
    COHORT =  df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID ==RING)]
    
    if(POPULATION == "Bird island"){
      
      df_draws_toBind$Intercept[which(df_draws_toBind$ID == RING)] = posterior_samples_finalModel$b_Intercept[DRAW]
      df_draws_toBind$Slope[which(df_draws_toBind$ID == RING)] = posterior_samples_finalModel$b_Age_std[DRAW]
      df_draws_toBind$Slope2[which(df_draws_toBind$ID == RING)] = posterior_samples_finalModel$b_IAge_stdE2[DRAW]
      
      df_draws_toBind$Intercept_cohortPop[which(df_draws_toBind$ID == RING)] = posterior_samples_finalModel[ DRAW, which(colnames(posterior_samples_finalModel) == paste0("r_cohortPop[", COHORT , "-Bird.island", ",Intercept]" ))]
      
    }
    if(POPULATION == "Kerguelen"){
      
      
      df_draws_toBind$Intercept[which(df_draws_toBind$ID == RING)] =  posterior_samples_finalModel$b_Intercept[DRAW]  + posterior_samples_finalModel$b_populationKerguelen[DRAW]
      df_draws_toBind$Slope[which(df_draws_toBind$ID == RING)] = posterior_samples_finalModel$b_Age_std[DRAW] + posterior_samples_finalModel$`b_Age_std:populationKerguelen`[DRAW]
      df_draws_toBind$Slope2[which(df_draws_toBind$ID == RING)] = posterior_samples_finalModel$b_IAge_stdE2[DRAW] + posterior_samples_finalModel$`b_IAge_stdE2:populationKerguelen`[DRAW]
      
      df_draws_toBind$Intercept_cohortPop[which(df_draws_toBind$ID == RING)] = posterior_samples_finalModel[ DRAW, which(colnames(posterior_samples_finalModel) == paste0("r_cohortPop[", COHORT , "-Kerguelen", ",Intercept]" ))]
      
    }
    
    
    df_draws_toBind$Intercept_randID[which(df_draws_toBind$ID == RING)] = posterior_samples_finalModel[ DRAW, which(colnames(posterior_samples_finalModel)== paste0("r_ID[", RING, ",Intercept]" ))]
    df_draws_toBind$Slope_randID[which(df_draws_toBind$ID == RING)] = posterior_samples_finalModel[ DRAW, which(colnames(posterior_samples_finalModel)== paste0("r_ID[", RING, ",Age_std]" ))]
    df_draws_toBind$Slope2_randID[which(df_draws_toBind$ID == RING)] = posterior_samples_finalModel[ DRAW, which(colnames(posterior_samples_finalModel)== paste0("r_ID[", RING, ",IAge_stdE2]" ))]
    
  }
  
  # bind datasets
  df_draws  =rbind(df_draws,df_draws_toBind )
}


df_draws$Intercept_corrected = df_draws$Intercept + df_draws$Intercept_cohortPop + df_draws$Intercept_randID
df_draws$Slope_corrected = df_draws$Slope + df_draws$Slope_randID  
df_draws$Slope2_corrected = df_draws$Slope2 + df_draws$Slope2_randID 

# recalculate parameters of interest
# onset of senescence
df_draws$OnsetSenescence = - df_draws$Slope_corrected / (2* df_draws$Slope2_corrected )

# senescence rate
# Find minimum of the derivative 
xinfl<-c()
gradmins<-c()
for(index in 1:length(df_draws$Slope2_corrected)){
  print(index)
  
  minval<-optimize(deriv_logit_quad,c(0,100),df_draws$Slope2_corrected[index],df_draws$Slope_corrected[index],df_draws$Intercept_corrected[index]) 
  xinfl<-c(xinfl,minval$minimum)
  gradmins<-c(gradmins,minval$objective)
}

# Add to df_draws  
df_draws$SenesceRate<-gradmins
df_draws$ageInfl<-xinfl
df_draws$SenesceRate = -df_draws$SenesceRate


# same from back-transformed values
df_draws$Intercept_backtransformed = back_trans_intercept(df_draws$Intercept_corrected, df_draws$Slope_corrected, df_draws$Slope2_corrected, mean_age_forStd, sd_age_forStd)
df_draws$Slope_backtransformed = back_trans_slope(df_draws$Intercept_corrected, df_draws$Slope_corrected, df_draws$Slope2_corrected, mean_age_forStd, sd_age_forStd)
df_draws$Slope2_backtransformed = back_trans_slope2(df_draws$Intercept_corrected, df_draws$Slope_corrected, df_draws$Slope2_corrected, mean_age_forStd, sd_age_forStd)

df_draws$OnsetSenescence_backtransformed = - df_draws$Slope_backtransformed / (2* df_draws$Slope2_backtransformed )

xinfl<-c()
gradmins<-c()
for(index in 1:length(df_draws$Slope2_backtransformed)){
  print(index)
  minval<-optimize(deriv_logit_quad,c(0,100),df_draws$Slope2_backtransformed[index],df_draws$Slope_backtransformed[index],df_draws$Intercept_backtransformed[index]) 
  xinfl<-c(xinfl,minval$minimum)
  gradmins<-c(gradmins,minval$objective)
}

# Add to df_draws  
df_draws$SenesceRate_backtransformed<-gradmins
df_draws$ageInfl_backtransformed<-xinfl
df_draws$SenesceRate_backtransformed = -df_draws$SenesceRate_backtransformed


# Make a summary of the results of the draws
df_draws_summary = data.frame(ID = unique(df_draws$ID),
                              population = rep(NA,length( unique(df_draws$ID))),
                              
                              SenesceRate = rep(NA,length( unique(df_draws$ID))),
                              SenesceRate_backtransformed = rep(NA,length( unique(df_draws$ID))),
                              sd_SenesceRate_backtransformed = rep(NA,length( unique(df_draws$ID))),
                              se_SenesceRate_backtransformed = rep(NA,length( unique(df_draws$ID))),
                              
                              OnsetSenescence = rep(NA,length( unique(df_draws$ID))),
                              OnsetSenescence_backtransformed = rep(NA,length( unique(df_draws$ID))),
                              sd_OnsetSenescence_backtransformed = rep(NA,length( unique(df_draws$ID))),
                              se_OnsetSenescence_backtransformed = rep(NA,length( unique(df_draws$ID))),
                              
                              SenesceRate_backtransformed_fullModel = rep(NA,length( unique(df_draws$ID))),
                              OnsetSenescence_backtransformed_fullModel = rep(NA,length( unique(df_draws$ID))),
                              SenesceRate_fullModel = rep(NA,length( unique(df_draws$ID))),
                              OnsetSenescence_fullModel = rep(NA,length( unique(df_draws$ID))),
                              
                              OnsetSenescence_backtransformed_mean = rep(NA,length( unique(df_draws$ID))),
                              SenesceRate_backtransformed_mean = rep(NA,length( unique(df_draws$ID))),
                              
                              
                              Intercept_fullModel = rep(NA,length( unique(df_draws$ID))),
                              Slope_fullModel = rep(NA,length( unique(df_draws$ID))),
                              Slope2_fullModel = rep(NA,length( unique(df_draws$ID))), 
                              
                              Intercept_median = rep(NA,length( unique(df_draws$ID))),
                              Slope_median = rep(NA,length( unique(df_draws$ID))),
                              Slope2_median = rep(NA,length( unique(df_draws$ID))), 
                              
                              Intercept_mean = rep(NA,length( unique(df_draws$ID))),
                              Slope_mean = rep(NA,length( unique(df_draws$ID))),
                              Slope2_mean = rep(NA,length( unique(df_draws$ID)))  
                              
)

for (k in 1:length(df_draws_summary$ID)){
  RING = df_draws_summary$ID[k]
  
  df_draws_summary$population[k]=unique(df_draws$population[which(df_draws$ID==RING)])
  
  df_draws_summary$SenesceRate[k]=median(df_draws$SenesceRate[which(df_draws$ID==RING)])
  df_draws_summary$OnsetSenescence[k]=median(df_draws$OnsetSenescence[which(df_draws$ID==RING)])
  
  df_draws_summary$SenesceRate_backtransformed[k]=median(df_draws$SenesceRate_backtransformed[which(df_draws$ID==RING)])
  df_draws_summary$OnsetSenescence_backtransformed[k]=median(df_draws$OnsetSenescence_backtransformed[which(df_draws$ID==RING)])
  
  df_draws_summary$SenesceRate_backtransformed_mean[k]=mean(df_draws$SenesceRate_backtransformed[which(df_draws$ID==RING)])
  df_draws_summary$OnsetSenescence_backtransformed_mean[k]=mean(df_draws$OnsetSenescence_backtransformed[which(df_draws$ID==RING)])
  
  
  df_draws_summary$sd_SenesceRate_backtransformed[k]=sd(df_draws$SenesceRate_backtransformed[which(df_draws$ID==RING)])
  df_draws_summary$sd_OnsetSenescence_backtransformed[k]=sd(df_draws$OnsetSenescence_backtransformed[which(df_draws$ID==RING)])
  
  df_draws_summary$se_SenesceRate_backtransformed[k]= df_draws_summary$sd_SenesceRate_backtransformed[k] /sqrt(length(df_draws$SenesceRate_backtransformed[which(df_draws$ID==RING)]))
  df_draws_summary$se_OnsetSenescence_backtransformed[k]=sd(df_draws$OnsetSenescence_backtransformed[which(df_draws$ID==RING)]) /sqrt(length(df_draws$SenesceRate_backtransformed[which(df_draws$ID==RING)]))
  
  
  df_draws_summary$SenesceRate_backtransformed_fullModel[k]=unique(df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID==RING)])
  df_draws_summary$OnsetSenescence_backtransformed_fullModel[k]=unique(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID==RING)])
  
  df_draws_summary$SenesceRate_fullModel[k]=unique(df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID==RING)])
  df_draws_summary$OnsetSenescence_fullModel[k]=unique(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID==RING)])
  
  
  df_draws_summary$Intercept_fullModel[k]=unique(df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID==RING)])
  df_draws_summary$Slope_fullModel[k]=unique(df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID==RING)])
  df_draws_summary$Slope2_fullModel[k]=unique(df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID==RING)])
  
  df_draws_summary$Intercept_median[k]=median(df_draws$Intercept_corrected[which(df_draws$ID==RING)])
  df_draws_summary$Slope_median[k]=median(df_draws$Slope_corrected[which(df_draws$ID==RING)])
  df_draws_summary$Slope2_median[k]=median(df_draws$Slope2_corrected[which(df_draws$ID==RING)])
  
  
  df_draws_summary$Intercept_mean[k]=mean(df_draws$Intercept_corrected[which(df_draws$ID==RING)])
  df_draws_summary$Slope_mean[k]=mean(df_draws$Slope_corrected[which(df_draws$ID==RING)])
  df_draws_summary$Slope2_mean[k]=mean(df_draws$Slope2_corrected[which(df_draws$ID==RING)])
  
}

 

#### _ df summary of iteration values ####

df_iterations_results = data.frame(draw = unique(df_draws$draw),
                                   corr_BI  = rep(NA,length(unique(df_draws$draw))),
                                   corr_Ker  = rep(NA,length(unique(df_draws$draw))),
                                   corr_OnsetInf100_BI  = rep(NA,length(unique(df_draws$draw))),
                                   corr_OnsetInf100_Ker  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_onset_BI  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_onset_Ker  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_onset_Diff  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_rate_BI  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_rate_Ker  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_rate_Diff  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_OnsetInf100_onset_BI  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_OnsetInf100_onset_Ker  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_OnsetInf100_onset_Diff  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_OnsetInf100_rate_BI  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_OnsetInf100_rate_Ker  = rep(NA,length(unique(df_draws$draw))),
                                   E_LRS_OnsetInf100_rate_Diff  = rep(NA,length(unique(df_draws$draw)))
)

#### _ correlations between onset and rate of senescence ####

for ( k in 1:length(df_iterations_results$draw)){
  DRAW = df_iterations_results$draw[k]
  
  df_iterations_results$corr_BI[which(df_iterations_results$draw==DRAW)] = cor(df_draws$OnsetSenescence_backtransformed[which(df_draws$population=="Bird island" &
                                                                                                                                df_draws$draw ==DRAW)],
                                                                               df_draws$SenesceRate_backtransformed[which(df_draws$population=="Bird island" &
                                                                                                                            df_draws$draw ==DRAW)])
  
  df_iterations_results$corr_Ker[which(df_iterations_results$draw==DRAW)] = cor(df_draws$OnsetSenescence_backtransformed[which(df_draws$population=="Kerguelen" &
                                                                                                                                 df_draws$draw ==DRAW)],
                                                                                df_draws$SenesceRate_backtransformed[which(df_draws$population=="Kerguelen" &
                                                                                                                             df_draws$draw ==DRAW)])
  df_iterations_results$corr_OnsetInf100_BI[which(df_iterations_results$draw==DRAW)] = cor(df_draws$OnsetSenescence_backtransformed[which(df_draws$population=="Bird island" &
                                                                                                                                            df_draws$draw ==DRAW &
                                                                                                                                            df_draws$OnsetSenescence_backtransformed < 100)],
                                                                                           df_draws$SenesceRate_backtransformed[which(df_draws$population=="Bird island" &
                                                                                                                                        df_draws$draw ==DRAW &
                                                                                                                                        df_draws$OnsetSenescence_backtransformed < 100)])
  
  df_iterations_results$corr_OnsetInf100_Ker[which(df_iterations_results$draw==DRAW)] = cor(df_draws$OnsetSenescence_backtransformed[which(df_draws$population=="Kerguelen" &
                                                                                                                                             df_draws$draw ==DRAW &
                                                                                                                                             df_draws$OnsetSenescence_backtransformed < 100)],
                                                                                            df_draws$SenesceRate_backtransformed[which(df_draws$population=="Kerguelen" &
                                                                                                                                         df_draws$draw ==DRAW &
                                                                                                                                         df_draws$OnsetSenescence_backtransformed < 100)])
  
}

hist(df_iterations_results$corr_BI, 100)
hist(df_iterations_results$corr_Ker, 100)

median(df_iterations_results$corr_BI)
median(df_iterations_results$corr_Ker)

mean(df_iterations_results$corr_BI)
mean(df_iterations_results$corr_Ker)


hist(df_iterations_results$corr_OnsetInf100_BI, 100)
hist(df_iterations_results$corr_OnsetInf100_Ker, 100)

median(df_iterations_results$corr_OnsetInf100_BI)
median(df_iterations_results$corr_OnsetInf100_Ker)

mean(df_iterations_results$corr_OnsetInf100_BI)
mean(df_iterations_results$corr_OnsetInf100_Ker)



#### _ LRS ~ onset and rate of senescence, potential outliers removed ####



for ( k in 1:length(df_iterations_results$draw)){
  print(k)
  
  DRAW = df_iterations_results$draw[k]
  
  
  df_draws_subsetDead_subsetIteration = df_draws_subsetDead[which(df_draws_subsetDead$draw ==DRAW &
                                                                    df_draws_subsetDead$OnsetSenescence_backtransformed < 40&
                                                                    df_draws_subsetDead$OnsetSenescence_backtransformed >0),]
  
  df_draws_subsetDead_subsetIteration$OnsetSenescence_std = robustHD::standardize(df_draws_subsetDead_subsetIteration$OnsetSenescence_backtransformed)
  df_draws_subsetDead_subsetIteration$SenesceRate_std = robustHD::standardize(df_draws_subsetDead_subsetIteration$SenesceRate_backtransformed)
  
  # ONSET OF SENESCENCE
  GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence  = brm(LRS   ~   OnsetSenescence_std * population  ,
                                                   data=df_draws_subsetDead_subsetIteration ,
                                                   family = poisson   ,
                                                   # warmup = 1500, 
                                                   # iter = 8000, 
                                                   warmup = 2000,
                                                   iter = 4000,
                                                   chains =4,
                                                   refresh = 0,
                                                   cores =4,
                                                   backend = "cmdstanr" ,
                                                   control = list(adapt_delta = 0.85)
  )
  
  df_iterations_results$E_LRS_OnsetInf100_onset_BI[which(df_iterations_results$draw==DRAW )] =  fixef(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence)[2]
  df_iterations_results$E_LRS_OnsetInf100_onset_Diff[which(df_iterations_results$draw==DRAW )] =  fixef(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence)[4]
  
  df_draws_subsetDead_subsetIteration$population=as.factor(df_draws_subsetDead_subsetIteration$population)
  df_draws_subsetDead_subsetIteration$population_refKer = relevel(df_draws_subsetDead_subsetIteration$population , ref=2)
  
  GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer  = brm(LRS   ~   OnsetSenescence_std * population_refKer  ,
                                                          data=df_draws_subsetDead_subsetIteration ,
                                                          family = poisson   ,
                                                          # warmup = 1500, 
                                                          # iter = 8000, 
                                                          warmup = 2000,
                                                          iter = 4000,
                                                          refresh = 0,
                                                          chains =4,
                                                          cores =4,
                                                          backend = "cmdstanr" ,
                                                          control = list(adapt_delta = 0.85)
  )
  
  df_iterations_results$E_LRS_OnsetInf100_onset_Ker[which(df_iterations_results$draw==DRAW )] =  fixef(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer)[2]
  
  # SENESCENCE RATE
  
  GLM_KerBI_Poisson_BRM_LRS_SenesceRate  = brm(LRS   ~   SenesceRate_std * population  ,
                                               data=df_draws_subsetDead_subsetIteration ,
                                               family = poisson   ,
                                               # warmup = 1500, 
                                               # iter = 8000, 
                                               warmup = 2000,
                                               iter = 4000,
                                               chains =4,
                                               refresh = 0,
                                               cores =4,
                                               backend = "cmdstanr" ,
                                               control = list(adapt_delta = 0.85)
  )
  
  df_iterations_results$E_LRS_OnsetInf100_rate_BI[which(df_iterations_results$draw==DRAW )] =  fixef(GLM_KerBI_Poisson_BRM_LRS_SenesceRate)[2]
  df_iterations_results$E_LRS_OnsetInf100_rate_Diff[which(df_iterations_results$draw==DRAW )] =  fixef(GLM_KerBI_Poisson_BRM_LRS_SenesceRate)[4]
  
  df_draws_subsetDead_subsetIteration$population=as.factor(df_draws_subsetDead_subsetIteration$population)
  df_draws_subsetDead_subsetIteration$population_refKer = relevel(df_draws_subsetDead_subsetIteration$population , ref=2)
  
  GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer  = brm(LRS   ~   SenesceRate_std * population_refKer  ,
                                                      data=df_draws_subsetDead_subsetIteration ,
                                                      family = poisson   ,
                                                      # warmup = 1500, 
                                                      # iter = 8000, 
                                                      warmup = 2000,
                                                      iter = 4000,
                                                      chains =4,
                                                      cores =4,
                                                      refresh = 0,
                                                      backend = "cmdstanr" ,
                                                      control = list(adapt_delta = 0.85)
  )
  
  df_iterations_results$E_LRS_OnsetInf100_rate_Ker[which(df_iterations_results$draw==DRAW )] =  fixef(GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer)[2] 
  
}



# EFFECT OF ONSET BI 
hist( df_iterations_results$E_LRS_OnsetInf100_onset_BI,100)

mean(df_iterations_results$E_LRS_OnsetInf100_onset_BI)
median(df_iterations_results$E_LRS_OnsetInf100_onset_BI)
quantile(df_iterations_results$E_LRS_OnsetInf100_onset_BI, 0.025)
quantile(df_iterations_results$E_LRS_OnsetInf100_onset_BI, 0.975)

# EFFECT OF ONSET KER 
hist( df_iterations_results$E_LRS_OnsetInf100_onset_Ker, 100)

mean(df_iterations_results$E_LRS_OnsetInf100_onset_Ker)
median(df_iterations_results$E_LRS_OnsetInf100_onset_Ker)
quantile(df_iterations_results$E_LRS_OnsetInf100_onset_Ker, 0.025)
quantile(df_iterations_results$E_LRS_OnsetInf100_onset_Ker, 0.975)

# ONSET DIFFERENCE BETWEEN SLOPES
hist( df_iterations_results$E_LRS_OnsetInf100_onset_Diff, 100)

mean(df_iterations_results$E_LRS_OnsetInf100_onset_Diff)
median(df_iterations_results$E_LRS_OnsetInf100_onset_Diff)
quantile(df_iterations_results$E_LRS_OnsetInf100_onset_Diff, 0.025)
quantile(df_iterations_results$E_LRS_OnsetInf100_onset_Diff, 0.975)



# EFFECT OF rate BI 
hist( df_iterations_results$E_LRS_OnsetInf100_rate_BI,100)

mean(df_iterations_results$E_LRS_OnsetInf100_rate_BI)
median(df_iterations_results$E_LRS_OnsetInf100_rate_BI)
quantile(df_iterations_results$E_LRS_OnsetInf100_rate_BI, 0.025)
quantile(df_iterations_results$E_LRS_OnsetInf100_rate_BI, 0.975)

# EFFECT OF rate KER 
hist( df_iterations_results$E_LRS_OnsetInf100_rate_Ker, 100)

mean(df_iterations_results$E_LRS_OnsetInf100_rate_Ker, 0.025)
median(df_iterations_results$E_LRS_OnsetInf100_rate_Ker, 0.975)
quantile(df_iterations_results$E_LRS_OnsetInf100_rate_Ker, 0.025)
quantile(df_iterations_results$E_LRS_OnsetInf100_rate_Ker, 0.975)

# rate DIFFERENCE BETWEEN SLOPES
hist( df_iterations_results$E_LRS_OnsetInf100_rate_Diff, 100)

mean(df_iterations_results$E_LRS_OnsetInf100_rate_Diff, 0.025)
median(df_iterations_results$E_LRS_OnsetInf100_rate_Diff, 0.975)
quantile(df_iterations_results$E_LRS_OnsetInf100_rate_Diff, 0.025)
quantile(df_iterations_results$E_LRS_OnsetInf100_rate_Diff, 0.975)


# Plot the results

df_summary_LRS_Onset_Rate_inf60 = data.frame(population = c("Kerguelen" , "Kerguelen", "Bird Island", "Bird Island"),
                                             parameter = c("onset", "rate","onset", "rate"),
                                             mean= rep(NA,4 ),
                                             se = rep(NA,4 ),
                                             sd = rep(NA,4 ),
                                             lower_ci = rep(NA,4 ),
                                             upper_Ci = rep(NA,4 ) )

df_summary_LRS_Onset_Rate_inf60$mean[which(df_summary_LRS_Onset_Rate_inf60$population =="Kerguelen" & df_summary_LRS_Onset_Rate_inf60$parameter =="onset")] = 
  mean(df_iterations_results$E_LRS_OnsetInf100_onset_Ker)
df_summary_LRS_Onset_Rate_inf60$mean[which(df_summary_LRS_Onset_Rate_inf60$population =="Bird Island" & df_summary_LRS_Onset_Rate_inf60$parameter =="onset")] = 
  mean(df_iterations_results$E_LRS_OnsetInf100_onset_BI)
df_summary_LRS_Onset_Rate_inf60$mean[which(df_summary_LRS_Onset_Rate_inf60$population =="Kerguelen" & df_summary_LRS_Onset_Rate_inf60$parameter =="rate")] = 
  mean(df_iterations_results$E_LRS_OnsetInf100_rate_Ker)
df_summary_LRS_Onset_Rate_inf60$mean[which(df_summary_LRS_Onset_Rate_inf60$population =="Bird Island" & df_summary_LRS_Onset_Rate_inf60$parameter =="rate")] = 
  mean(df_iterations_results$E_LRS_OnsetInf100_rate_BI)

df_summary_LRS_Onset_Rate_inf60$sd[which(df_summary_LRS_Onset_Rate_inf60$population =="Kerguelen" & df_summary_LRS_Onset_Rate_inf60$parameter =="onset")] = 
  sd(df_iterations_results$E_LRS_OnsetInf100_onset_Ker)
df_summary_LRS_Onset_Rate_inf60$sd[which(df_summary_LRS_Onset_Rate_inf60$population =="Bird Island" & df_summary_LRS_Onset_Rate_inf60$parameter =="onset")] = 
  sd(df_iterations_results$E_LRS_OnsetInf100_onset_BI)
df_summary_LRS_Onset_Rate_inf60$sd[which(df_summary_LRS_Onset_Rate_inf60$population =="Kerguelen" & df_summary_LRS_Onset_Rate_inf60$parameter =="rate")] = 
  sd(df_iterations_results$E_LRS_OnsetInf100_rate_Ker)
df_summary_LRS_Onset_Rate_inf60$sd[which(df_summary_LRS_Onset_Rate_inf60$population =="Bird Island" & df_summary_LRS_Onset_Rate_inf60$parameter =="rate")] = 
  sd(df_iterations_results$E_LRS_OnsetInf100_rate_BI)


df_summary_LRS_Onset_Rate_inf60$lower_ci[which(df_summary_LRS_Onset_Rate_inf60$population =="Kerguelen" & df_summary_LRS_Onset_Rate_inf60$parameter =="onset")] = 
  quantile(df_iterations_results$E_LRS_OnsetInf100_onset_Ker, 0.025)
df_summary_LRS_Onset_Rate_inf60$lower_ci[which(df_summary_LRS_Onset_Rate_inf60$population =="Bird Island" & df_summary_LRS_Onset_Rate_inf60$parameter =="onset")] = 
  quantile(df_iterations_results$E_LRS_OnsetInf100_onset_BI, 0.025)
df_summary_LRS_Onset_Rate_inf60$lower_ci[which(df_summary_LRS_Onset_Rate_inf60$population =="Kerguelen" & df_summary_LRS_Onset_Rate_inf60$parameter =="rate")] = 
  quantile(df_iterations_results$E_LRS_OnsetInf100_rate_Ker, 0.025)
df_summary_LRS_Onset_Rate_inf60$lower_ci[which(df_summary_LRS_Onset_Rate_inf60$population =="Bird Island" & df_summary_LRS_Onset_Rate_inf60$parameter =="rate")] = 
  quantile(df_iterations_results$E_LRS_OnsetInf100_rate_BI, 0.025)

df_summary_LRS_Onset_Rate_inf60$upper_ci[which(df_summary_LRS_Onset_Rate_inf60$population =="Kerguelen" & df_summary_LRS_Onset_Rate_inf60$parameter =="onset")] = 
  quantile(df_iterations_results$E_LRS_OnsetInf100_onset_Ker, 0.975)
df_summary_LRS_Onset_Rate_inf60$upper_ci[which(df_summary_LRS_Onset_Rate_inf60$population =="Bird Island" & df_summary_LRS_Onset_Rate_inf60$parameter =="onset")] = 
  quantile(df_iterations_results$E_LRS_OnsetInf100_onset_BI, 0.975)
df_summary_LRS_Onset_Rate_inf60$upper_ci[which(df_summary_LRS_Onset_Rate_inf60$population =="Kerguelen" & df_summary_LRS_Onset_Rate_inf60$parameter =="rate")] = 
  quantile(df_iterations_results$E_LRS_OnsetInf100_rate_Ker, 0.975)
df_summary_LRS_Onset_Rate_inf60$upper_ci[which(df_summary_LRS_Onset_Rate_inf60$population =="Bird Island" & df_summary_LRS_Onset_Rate_inf60$parameter =="rate")] = 
  quantile(df_iterations_results$E_LRS_OnsetInf100_rate_BI, 0.975)

ggplot(data = df_summary_LRS_Onset_Rate_inf60,
       aes(x = mean, y = parameter, col = population,
           xmin = mean - sd, xmax = mean + sd))+
  geom_pointrange(position = position_dodge(width = 0.3),
                  fatten = 2, linewidth = 1.2)+
  geom_pointrange(inherit.aes = FALSE,
                  aes(x = mean, y = parameter, col = population,
                      xmin = lower_ci, xmax = upper_ci),
                  position = position_dodge(width = 0.3))+
  geom_vline(xintercept = 0, linetype = "dashed", col ="red")+
  theme_classic()


ggsave(paste0(dir_data, "SupportingInformation_S6_Figure_S6_1.png"),
       plot = pmarginal)



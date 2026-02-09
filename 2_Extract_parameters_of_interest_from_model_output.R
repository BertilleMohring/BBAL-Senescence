#### Extract parameters of interest (senescence rate and onset of senescence) from model output ####

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


# Load functions used later
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



# Open dataset

data_BBAL = read.csv( paste0(dir_data_BBAL, "data_BBAL.csv"))
data_BBAL=data_BBAL[,-1]

# Load model output
load( paste0(dir_data, "model_output_GLMM_senescence_priors.RData"))


# calculate age mean and sd
mean_age_forStd =  mean(data_BBAL$Age) 
sd_age_forStd =  sd(data_BBAL$Age) 

length(data_BBAL$population)


#### 2. Extract senescence patterns from the model ####

#### __ Extract random slopes and intercepts #### 

# Extract random intercepts for cohorts (for each cohort from each population)
Random_AgeStd_Intercept_cohort=  ranef(GLMM_brms_senescence)$cohortPop[, , "Intercept"] 
Random_AgeStd_Intercept_cohort=as.data.frame(Random_AgeStd_Intercept_cohort)
Random_AgeStd_Intercept_cohort$cohortPop= rownames(Random_AgeStd_Intercept_cohort)
Random_AgeStd_Intercept_cohort$population=substr(Random_AgeStd_Intercept_cohort$cohortPop, 6, nchar(Random_AgeStd_Intercept_cohort$cohortPop))
Random_AgeStd_Intercept_cohort$cohort=substr(Random_AgeStd_Intercept_cohort$cohortPop, 1,4)

head(Random_AgeStd_Intercept_cohort)

# Extract random slopes and intercepts on individuals

# random intercepts
Random_AgeStd_Intercept_ID=  ranef(GLMM_brms_senescence)$ID[, , "Intercept"] 

Random_AgeStd_Intercept_ID = as.data.frame(Random_AgeStd_Intercept_ID)
Random_AgeStd_Intercept_ID$ID= rownames(Random_AgeStd_Intercept_ID)
data_BBAL$ID = as.character(data_BBAL$ID)

# add information on population and cohort to the dataset
Random_AgeStd_Intercept_ID=merge(Random_AgeStd_Intercept_ID,
                                   unique(data_BBAL[,which(colnames(data_BBAL) %in% c("ID", "population", "cohortPop"))]),
                                   by = "ID",
                                   all.x = TRUE, all.y = FALSE) 
# extract intercepts
InterceptModel_ref=fixef(GLMM_brms_senescence)[1]
InterceptModel_Ker=fixef(GLMM_brms_senescence)[5]

# Calculate the intercept for each individual (fixed effect + fixed effect of population effect for Kerguelen + random intercept on individual)
Random_AgeStd_Intercept_ID$Estimate_corrected = Random_AgeStd_Intercept_ID$Estimate + fixef(GLMM_brms_senescence)[1]
Random_AgeStd_Intercept_ID$Estimate_corrected[which(Random_AgeStd_Intercept_ID$population =="Kerguelen" )] =   Random_AgeStd_Intercept_ID$Estimate_corrected[which(Random_AgeStd_Intercept_ID$population =="Kerguelen" )] + fixef(GLMM_brms_senescence)[5]

# add cohort effect to the intercept value
for (k in 1:length(unique(Random_AgeStd_Intercept_ID$ID))){
  ID = unique(Random_AgeStd_Intercept_ID$ID)[k]
  COHORT = Random_AgeStd_Intercept_ID$cohortPop[which(Random_AgeStd_Intercept_ID$ID==ID)]
  
  Random_AgeStd_Intercept_ID$Estimate_corrected[which(Random_AgeStd_Intercept_ID$ID==ID)]=  Random_AgeStd_Intercept_ID$Estimate_corrected[which(Random_AgeStd_Intercept_ID$ID==ID)]+ Random_AgeStd_Intercept_cohort$Estimate[which(Random_AgeStd_Intercept_cohort$cohortPop==COHORT)]
  
}
 
# random slopes on linear age effect
# extract the slope for each individual
Random_AgeStd_Slope_ID=  ranef(GLMM_brms_senescence)$ID[, , "Age_std"] 

Random_AgeStd_Slope_ID = as.data.frame(Random_AgeStd_Slope_ID)
Random_AgeStd_Slope_ID$ID= rownames(Random_AgeStd_Slope_ID)

# add information on population
Random_AgeStd_Slope_ID=merge(Random_AgeStd_Slope_ID,
                               unique(data_BBAL[,which(colnames(data_BBAL) %in% c("ID", "population"))]),
                               by = "ID",
                               all.x = TRUE, all.y = FALSE)

# extract the slope values for each population
SlopeModel_ref=fixef(GLMM_brms_senescence)[2]
SlopeModel_Ker=fixef(GLMM_brms_senescence)[6]

# calculate the linear age effect for each individaul (fixed effect + random effect + fixed population effect)
Random_AgeStd_Slope_ID$Estimate_corrected = Random_AgeStd_Slope_ID$Estimate + fixef(GLMM_brms_senescence)[2]
Random_AgeStd_Slope_ID$Estimate_corrected[which(Random_AgeStd_Slope_ID$population =="Kerguelen" )] =   Random_AgeStd_Slope_ID$Estimate_corrected[which(Random_AgeStd_Slope_ID$population =="Kerguelen" )] + fixef(GLMM_brms_senescence)[6]
 


# random slopes on quadratic age effect
Random_AgeStd_Slope_sq_ID=  ranef(GLMM_brms_senescence)$ID[, , "IAge_stdE2"] 
# extract the slope for each individaul
Random_AgeStd_Slope_sq_ID = as.data.frame(Random_AgeStd_Slope_sq_ID)
Random_AgeStd_Slope_sq_ID$ID= rownames(Random_AgeStd_Slope_sq_ID)
# add information on population
Random_AgeStd_Slope_sq_ID=merge(Random_AgeStd_Slope_sq_ID,
                                  unique(data_BBAL[,which(colnames(data_BBAL) %in% c("ID", "population"))]),
                                  by = "ID",
                                  all.x = TRUE, all.y = FALSE)

SlopeSqModel_ref=fixef(GLMM_brms_senescence)[3]
SlopeSqModel_Ker=fixef(GLMM_brms_senescence)[7]
# calculate the quadratic age effect for each individaul (fixed effect + random effect + fixed population effect)
Random_AgeStd_Slope_sq_ID$Estimate_corrected = Random_AgeStd_Slope_sq_ID$Estimate + fixef(GLMM_brms_senescence)[3]
Random_AgeStd_Slope_sq_ID$Estimate_corrected[which(Random_AgeStd_Slope_sq_ID$population =="Kerguelen" )] =   Random_AgeStd_Slope_sq_ID$Estimate_corrected[which(Random_AgeStd_Slope_sq_ID$population =="Kerguelen" )] + fixef(GLMM_brms_senescence)[7]

 
#### __  Make a dataframe with individual-level values for each parameter, predicted by the model #### 

df_AgeStd_randomSlopesAndInterceptsIDlevel =data.frame(Intercept = rep(NA, length( Random_AgeStd_Slope_ID$ID)),
                                                       Slope = Random_AgeStd_Slope_ID$Estimate_corrected ,
                                                       Slope2 =rep(NA, length( Random_AgeStd_Slope_ID$ID)),   
                                                       OnsetSenescence =  rep(NA, length( Random_AgeStd_Slope_ID$ID)),
                                                       ID = Random_AgeStd_Slope_ID$ID,
                                                       population = Random_AgeStd_Slope_ID$population,
                                                       cohort =rep(NA, length( Random_AgeStd_Slope_ID$ID)),
                                                       cohortPop =rep(NA, length( Random_AgeStd_Slope_ID$ID)),
                                                       Year_last_seen =rep(NA, length( Random_AgeStd_Slope_ID$ID)),
                                                       LRS = rep(NA, length( Random_AgeStd_Slope_ID$ID)),
                                                       
                                                       nb_BA_life = rep(NA, length( Random_AgeStd_Slope_ID$ID)) )

# fill in missing data
for (k in 1:length(df_AgeStd_randomSlopesAndInterceptsIDlevel$ID)){ # for each ID
  ID = df_AgeStd_randomSlopesAndInterceptsIDlevel$ID[k]
  df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] = Random_AgeStd_Intercept_ID$Estimate_corrected[which(Random_AgeStd_Intercept_ID$ID == ID)]
  df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] = Random_AgeStd_Slope_sq_ID$Estimate_corrected[which(Random_AgeStd_Slope_sq_ID$ID == ID)]
  
  df_AgeStd_randomSlopesAndInterceptsIDlevel$Year_last_seen[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] = unique(data_BBAL$Year_last_seen[which(data_BBAL$ID == ID)]) 
  df_AgeStd_randomSlopesAndInterceptsIDlevel$LRS[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] = unique(data_BBAL$lifetime_repro_output[which(data_BBAL$ID == ID)]) 
  df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] =  unique(data_BBAL$cohort[which(data_BBAL$ID == ID)]) 
  df_AgeStd_randomSlopesAndInterceptsIDlevel$cohortPop[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] =  unique(data_BBAL$cohortPop[which(data_BBAL$ID == ID)]) 
  df_AgeStd_randomSlopesAndInterceptsIDlevel$nb_BA_life[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] = unique(data_BBAL$nb_breeding_attempts[which(data_BBAL$ID == ID)]) 
  
}


#### __ calculate parameters of interest for each individual #### 

# Onset of senescence
df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenescence = - df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope / (2*df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2 )


# predicted RS values at the onset of senescence and at 10yo
df_AgeStd_randomSlopesAndInterceptsIDlevel$PredictedRS_atpeak = exp(df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2 * df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenescence ^2 + df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope * df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenescence  + df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept) / (1 + exp(df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2 * I(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenescence)^2 + df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope *  df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenescence + df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept))
df_AgeStd_randomSlopesAndInterceptsIDlevel$PredictedRS_at10yo = exp(df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2 * I((10-mean_age_forStd )/sd_age_forStd)^2 + df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope * ((10-mean_age_forStd )/sd_age_forStd)   + df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept) / (1 + exp(df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2 * I( ((10-mean_age_forStd )/sd_age_forStd) )^2 + df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope *   ((10-mean_age_forStd )/sd_age_forStd)  + df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept))


# extract rate of senescence
# Function that calculates the derivative of 
# exp(quad*xval^2+lin*xval+const)/(1+exp(quad*xval^2+lin*xval+const))
deriv_logit_quad<-function(xval,quad,lin,const)
{
  output<-(2*quad*xval+lin)*exp(quad*xval^2+lin*xval+const)/((1+exp(quad*xval^2+lin*xval+const))^2)
  return(output)
}
# Find minimum of the derivative 
xinfl<-c()
gradmins<-c()
for(index in 1:length(df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2))
{
  minval<-optimize(deriv_logit_quad,c(0,100),df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2[index],df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope[index],df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept[index]) 
  xinfl<-c(xinfl,minval$minimum)
  gradmins<-c(gradmins,minval$objective)
}


# Add to dataframe data and put out to CSV
df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate<-gradmins
df_AgeStd_randomSlopesAndInterceptsIDlevel$ageInfl<-xinfl

# multiply senescence rate value by -1 so that higher values correspond to stronger senescence rate
df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate = -df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate



#### ____ Back-transform parameters #### 

# back-transform quadratic slope term
df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2_backtransformed = back_trans_slope2(df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept, 
                                                                                      df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope, 
                                                                                      df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2 , 
                                                                                      mean_age_forStd, sd_age_forStd)


# back-transform linear slope term
df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope_backtransformed = back_trans_slope(df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept, 
                                                                                    df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope, 
                                                                                    df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2 , 
                                                                                    mean_age_forStd, sd_age_forStd)

# back-transform intercept term
df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercet_backtransformed = back_trans_intercept(df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept, 
                                                                                           df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope, 
                                                                                           df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2 , 
                                                                                           mean_age_forStd, sd_age_forStd)

# re-calculate onset of senescence using back-transformed coefficients
df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs = - df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope_backtransformed / (2*df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2_backtransformed)


df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc=as.numeric(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenescence)
df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_backTransformed = (df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc * sd_age_forStd ) + mean_age_forStd


# extract senescence rate from backtransformed coef
# Find minimum of the derivative 
xinfl<-c()
gradmins<-c()
for(index in 1:length(df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2))
{
  minval<-optimize(deriv_logit_quad,c(0,100),
                   df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope2_backtransformed[index],
                   df_AgeStd_randomSlopesAndInterceptsIDlevel$Slope_backtransformed[index],
                   df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercet_backtransformed[index]) 
  xinfl<-c(xinfl,minval$minimum)
  gradmins<-c(gradmins,minval$objective)
}

# Add to df_AgeStd_randomSlopesAndInterceptsIDlevelframe and put out to CSV
df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs<-gradmins

 

df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept = as.numeric(as.character(df_AgeStd_randomSlopesAndInterceptsIDlevel$Intercept))
df_AgeStd_randomSlopesAndInterceptsIDlevel$LRS = as.numeric(as.character(df_AgeStd_randomSlopesAndInterceptsIDlevel$LRS))
df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_backTransformed = as.numeric(as.character(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_backTransformed))
df_AgeStd_randomSlopesAndInterceptsIDlevel$PredictedRS_atpeak = as.numeric(as.character(df_AgeStd_randomSlopesAndInterceptsIDlevel$PredictedRS_atpeak))
df_AgeStd_randomSlopesAndInterceptsIDlevel$PredictedRS_at10yo = as.numeric(as.character(df_AgeStd_randomSlopesAndInterceptsIDlevel$PredictedRS_at10yo))
df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate = as.numeric(as.character(df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate))



#### ____ Save output #### 

write.csv(df_AgeStd_randomSlopesAndInterceptsIDlevel,
          paste0(dir_data, "data_BBAL_senescence_estimates_priors.csv"))





#### __ calculate parameters of interest from posterior distribution #### 

#extract posterior values from the model
posterior_samples_finalModel <- posterior_samples(GLMM_brms_senescence)

# prepare dataframe for intercept
df_plotRes_Intercept = data.frame(Value = rep(NA, length(posterior_samples_finalModel$b_Intercept)*2 ),
                                  Variable = rep("Intercept", length(posterior_samples_finalModel$b_Intercept)*2 ),
                                  population  = rep(c("Kerguelen", "Bird island"), length(posterior_samples_finalModel$b_Intercept)) )
# fill in dataframe
df_plotRes_Intercept$Value[which(df_plotRes_Intercept$population ==  "Bird island" )] = posterior_samples_finalModel$b_Intercept
df_plotRes_Intercept$Value[which(df_plotRes_Intercept$population ==  "Kerguelen"  )] = posterior_samples_finalModel$b_Intercept + posterior_samples_finalModel$b_populationKerguelen

# prepare dataframe for slope on linear age effect
df_plotRes_Age = data.frame(Value = rep(NA, length(posterior_samples_finalModel$b_Age_std)*2 ),
                            Variable = rep("Age", length(posterior_samples_finalModel$b_Intercept)*2 ),
                            population  = rep(c("Kerguelen", "Bird island"), length(posterior_samples_finalModel$b_Age_std)) )
# fill in dataframe
df_plotRes_Age$Value[which(df_plotRes_Age$population ==  "Bird island"  )] = posterior_samples_finalModel$b_Age_std
df_plotRes_Age$Value[which(df_plotRes_Age$population ==  "Kerguelen" )] = posterior_samples_finalModel$b_Age_std + posterior_samples_finalModel$`b_Age_std:populationKerguelen`

# prepare dataframe for slope on quadratic age effect
df_plotRes_IAgeE2 = data.frame(Value = rep(NA, length(posterior_samples_finalModel$b_IAge_stdE2)*2 ),
                               Variable = rep("IAgeE2", length(posterior_samples_finalModel$b_Intercept)*2 ),
                               population  = rep(c("Kerguelen", "Bird island"), length(posterior_samples_finalModel$b_IAge_stdE2)) )
# fill in dataframe
df_plotRes_IAgeE2$Value[which(df_plotRes_IAgeE2$population ==  "Bird island")] = posterior_samples_finalModel$b_IAge_stdE2
df_plotRes_IAgeE2$Value[which(df_plotRes_IAgeE2$population ==  "Kerguelen" )] = posterior_samples_finalModel$b_IAge_stdE2 + posterior_samples_finalModel$`b_IAge_stdE2:populationKerguelen`

# bind these datasets into one
df_plotRes_All = rbind(df_plotRes_Intercept,
                       df_plotRes_Age,
                       df_plotRes_IAgeE2)

# make a new dataframe for senescence rate and onset of senescence
df_plotRes_Senescence= data.frame(population =df_plotRes_Intercept$population,
                                  SenesceRate = rep(NA, length(df_plotRes_Intercept$Value)) ,
                                  OnsetSenesc= rep(NA, length(df_plotRes_Intercept$Value)) )
# calculate onset of senescence
df_plotRes_Senescence$OnsetSenesc = - df_plotRes_Age$Value / (2* df_plotRes_IAgeE2$Value)

# Calculate senescence rate
# Find minimum of the derivative 
xinfl<-c()
gradmins<-c()
for(index in 1:length(df_plotRes_Intercept$population))
{
  minval<-optimize(deriv_logit_quad,c(0,100),
                   df_plotRes_IAgeE2$Value[index],
                   df_plotRes_Age$Value[index],
                   df_plotRes_Intercept$Value[index]) 
  xinfl<-c(xinfl,minval$minimum)
  gradmins<-c(gradmins,minval$objective)
}
df_plotRes_Senescence$SenesceRate<-gradmins

# calculate predicted values of reproductive success at the onset of senescence and at 10 yo
df_plotRes_Senescence$PredictedRS_atpeak = exp(df_plotRes_IAgeE2$Value * df_plotRes_Senescence$OnsetSenesc ^2 + df_plotRes_Age$Value * df_plotRes_Senescence$OnsetSenesc  + df_plotRes_Intercept$Value) / (1 + exp(df_plotRes_IAgeE2$Value * I(df_plotRes_Senescence$OnsetSenesc)^2 + df_plotRes_Age$Value *  df_plotRes_Senescence$OnsetSenesc + df_plotRes_Intercept$Value))
df_plotRes_Senescence$PredictedRS_at10yo = exp(df_plotRes_IAgeE2$Value * ((10-mean_age_forStd )/sd_age_forStd)^2 + df_plotRes_Age$Value * ((10-mean_age_forStd )/sd_age_forStd)  + df_plotRes_Intercept$Value) / (1 + exp(df_plotRes_IAgeE2$Value * I(((10-mean_age_forStd )/sd_age_forStd))^2 + df_plotRes_Age$Value *  ((10-mean_age_forStd )/sd_age_forStd) + df_plotRes_Intercept$Value))

# multiply senescence rate value by -1 so that higher values correspond to stronger senescence rate
df_plotRes_Senescence$SenesceRate = - df_plotRes_Senescence$SenesceRate
 


#### ____ Back-transform parameters #### 
 
df_plotRes_Senescence$OnsetSenesc_backTransformed = (df_plotRes_Senescence$OnsetSenesc * sd_age_forStd) + mean_age_forStd


# back-transform quadratic slope term
df_plotRes_Senescence$Slope2_backtransformed = back_trans_slope2(df_plotRes_Intercept$Value, 
                                                                                      df_plotRes_Age$Value, 
                                                                                      df_plotRes_IAgeE2$Value , 
                                                                                      mean_age_forStd, sd_age_forStd)


# back-transform linear slope term
df_plotRes_Senescence$Slope_backtransformed = back_trans_slope(df_plotRes_Intercept$Value, 
                                                                                    df_plotRes_Age$Value, 
                                                                                    df_plotRes_IAgeE2$Value , 
                                                                                    mean_age_forStd, sd_age_forStd)

# back-transform intercept term
df_plotRes_Senescence$Intercept_backtransformed = back_trans_intercept(df_plotRes_Intercept$Value, 
                                                                                           df_plotRes_Age$Value, 
                                                                                           df_plotRes_IAgeE2$Value , 
                                                                                           mean_age_forStd, sd_age_forStd)

# re-calculate onset of senescence using back-transformed coefficients
df_plotRes_Senescence$OnsetSenesc_from_backTransformedCoefs = - df_plotRes_Senescence$Slope_backtransformed / (2*df_plotRes_Senescence$Slope2_backtransformed)



# extract senescence rate from backtransformed coef
# Find minimum of the derivative 
xinfl<-c()
gradmins<-c()
for(index in 1:length(df_plotRes_Senescence$Slope2))
{
  minval<-optimize(deriv_logit_quad,c(0,100),
                   df_plotRes_Senescence$Slope2_backtransformed[index],
                   df_plotRes_Senescence$Slope_backtransformed[index],
                   df_plotRes_Senescence$Intercept_backtransformed[index]) 
  xinfl<-c(xinfl,minval$minimum)
  gradmins<-c(gradmins,minval$objective)
}

# Add to df_plotRes_Senescenceframe and put out to CSV
df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs<-gradmins

# multiply senescence rate value by -1 so that higher values correspond to stronger senescence rate
df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs = -df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs


#### ____ Save output #### 
write.csv(df_plotRes_Senescence,
          paste0(dir_data, "data_BBAL_senescence_posteriorEstimates_priors.csv"))



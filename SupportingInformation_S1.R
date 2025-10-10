
#### Supporting Information S1  : Robustness test for variation in life-history strategies between populations and individuals  #### 
#### using a subset of individual with at least four breeding attempts ####

# Load packages

library(brms)
library(rstan)
library(StanHeaders)
library(bayesplot)
library(robustHD)
library(cmdstanr)

# Open dataset
data_BBAL = read.csv( paste0(dir_data, "data_BBAL.csv"))
data_BBAL=data_BBAL[,-1]

data_BBAL$nbBA = NA
for (k in 1:length(unique(data_BBAL$ID))){
  INDIVIDUAL = unique(data_BBAL$ID)[k]
  
  data_BBAL$nbBA[which(data_BBAL$ID == INDIVIDUAL)] =length( data_BBAL$ReproS[which(data_BBAL$ID == INDIVIDUAL)])
}
table(data_BBAL$nbBA, data_BBAL$population)


data_BBAL_atLeast4BA = data_BBAL[which(data_BBAL$nbBA >=4),]

table(data_BBAL_atLeast4BA$population)
length(unique(data_BBAL_atLeast4BA$ID[which(data_BBAL_atLeast4BA$population=="Kerguelen")]))
length(unique(data_BBAL_atLeast4BA$ID[which(data_BBAL_atLeast4BA$population=="Bird island")]))

# standardise variables
data_BBAL_atLeast4BA$Age_std = robustHD::standardize(data_BBAL_atLeast4BA$Age)

mean_age_forStd =  mean(data_BBAL_atLeast4BA$Age) 
sd_age_forStd =  sd(data_BBAL_atLeast4BA$Age) 



#### .####

#### PART 1: Senescence GLMM ####


prior_senescence <- c(
  # Fixed effects
  prior(normal(0, 5), class = "Intercept"),
  prior(normal(0, 2), class = "b"),
  
  # Group-level standard deviations
  prior(student_t(3, 0, 1), class = "sd"),
  
  # Correlation matrices
  prior(lkj(2), class = "cor")
)

# run model
GLMM_brms_senescence_atLeast4BA=brm(ReproS   ~   Age_std  + 
                           I(Age_std ^2)  +
                           First_breeding_attempt +  
                           population+
                           population:Age_std + 
                           population:I(Age_std ^2)  +    
                           (1 + Age_std +  I(Age_std ^2)|gr(ID, by = population)) +
                           (1| gr(YearPop, by = population))+
                           (1|gr(cohortPop, by = population)) ,
                         data=data_BBAL_atLeast4BA,
                         family = bernoulli(link = "logit"), 
                         warmup = 2000,
                         iter = 4000,
                         chains =8,
                         cores = 8,
                         prior = prior_senescence,
                         backend = "cmdstanr"  ,
                         control = list(adapt_delta = 0.90)
)

summary(GLMM_brms_senescence_atLeast4BA)


# save output
save(GLMM_brms_senescence_atLeast4BA, 
     file = paste0(dir_data, "model_output_GLMM_senescence_appendix_atLeast4BA_priors.RData"))





#### .####

#### PART 2:Extract parameters of interest (senescence rate and onset of senescence) from model output####


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



#### _2.1. Extract senescence patterns from the model ####

#### __ Extract random slopes and intercepts #### 

# Extract random intercepts for cohorts (for each cohort from each population)
Random_AgeStd_Intercept_cohort=  ranef(GLMM_brms_senescence_atLeast4BA)$cohortPop[, , "Intercept"] 
Random_AgeStd_Intercept_cohort=as.data.frame(Random_AgeStd_Intercept_cohort)
Random_AgeStd_Intercept_cohort$cohortPop= rownames(Random_AgeStd_Intercept_cohort)
Random_AgeStd_Intercept_cohort$population=substr(Random_AgeStd_Intercept_cohort$cohortPop, 6, nchar(Random_AgeStd_Intercept_cohort$cohortPop))
Random_AgeStd_Intercept_cohort$cohort=substr(Random_AgeStd_Intercept_cohort$cohortPop, 1,4)

head(Random_AgeStd_Intercept_cohort)

# Extract random slopes and intercepts on individuals

# random intercepts
Random_AgeStd_Intercept_ID=  ranef(GLMM_brms_senescence_atLeast4BA)$ID[, , "Intercept"] 

Random_AgeStd_Intercept_ID = as.data.frame(Random_AgeStd_Intercept_ID)
Random_AgeStd_Intercept_ID$ID= rownames(Random_AgeStd_Intercept_ID)
data_BBAL_atLeast4BA$ID = as.character(data_BBAL_atLeast4BA$ID)

# add information on population and cohort to the dataset
Random_AgeStd_Intercept_ID=merge(Random_AgeStd_Intercept_ID,
                                 unique(data_BBAL_atLeast4BA[,which(colnames(data_BBAL_atLeast4BA) %in% c("ID", "population", "cohortPop"))]),
                                 by = "ID",
                                 all.x = TRUE, all.y = FALSE) 
# extract intercepts
InterceptModel_ref=fixef(GLMM_brms_senescence_atLeast4BA)[1]
InterceptModel_Ker=fixef(GLMM_brms_senescence_atLeast4BA)[5]

# Calculate the intercept for each individual (fixed effect + fixed effect of population effect for Kerguelen + random intercept on individual)
Random_AgeStd_Intercept_ID$Estimate_corrected = Random_AgeStd_Intercept_ID$Estimate + fixef(GLMM_brms_senescence_atLeast4BA)[1]
Random_AgeStd_Intercept_ID$Estimate_corrected[which(Random_AgeStd_Intercept_ID$population =="Kerguelen" )] =   Random_AgeStd_Intercept_ID$Estimate_corrected[which(Random_AgeStd_Intercept_ID$population =="Kerguelen" )] + fixef(GLMM_brms_senescence_atLeast4BA)[5]

# add cohort effect to the intercept value
for (k in 1:length(unique(Random_AgeStd_Intercept_ID$ID))){
  ID = unique(Random_AgeStd_Intercept_ID$ID)[k]
  COHORT = Random_AgeStd_Intercept_ID$cohortPop[which(Random_AgeStd_Intercept_ID$ID==ID)]
  
  Random_AgeStd_Intercept_ID$Estimate_corrected[which(Random_AgeStd_Intercept_ID$ID==ID)]=  Random_AgeStd_Intercept_ID$Estimate_corrected[which(Random_AgeStd_Intercept_ID$ID==ID)]+ Random_AgeStd_Intercept_cohort$Estimate[which(Random_AgeStd_Intercept_cohort$cohortPop==COHORT)]
  
}

# random slopes on linear age effect
# extract the slope for each individaul
Random_AgeStd_Slope_ID=  ranef(GLMM_brms_senescence_atLeast4BA)$ID[, , "Age_std"] 

Random_AgeStd_Slope_ID = as.data.frame(Random_AgeStd_Slope_ID)
Random_AgeStd_Slope_ID$ID= rownames(Random_AgeStd_Slope_ID)

# add information on population
Random_AgeStd_Slope_ID=merge(Random_AgeStd_Slope_ID,
                             unique(data_BBAL_atLeast4BA[,which(colnames(data_BBAL_atLeast4BA) %in% c("ID", "population"))]),
                             by = "ID",
                             all.x = TRUE, all.y = FALSE)

# extract the slope values for each population
SlopeModel_ref=fixef(GLMM_brms_senescence_atLeast4BA)[2]
SlopeModel_Ker=fixef(GLMM_brms_senescence_atLeast4BA)[6]

# calculate the linear age effect for each individaul (fixed effect + random effect + fixed population effect)
Random_AgeStd_Slope_ID$Estimate_corrected = Random_AgeStd_Slope_ID$Estimate + fixef(GLMM_brms_senescence_atLeast4BA)[2]
Random_AgeStd_Slope_ID$Estimate_corrected[which(Random_AgeStd_Slope_ID$population =="Kerguelen" )] =   Random_AgeStd_Slope_ID$Estimate_corrected[which(Random_AgeStd_Slope_ID$population =="Kerguelen" )] + fixef(GLMM_brms_senescence_atLeast4BA)[6]



# random slopes on quadratic age effect
Random_AgeStd_Slope_sq_ID=  ranef(GLMM_brms_senescence_atLeast4BA)$ID[, , "IAge_stdE2"] 
# extract the slope for each individaul
Random_AgeStd_Slope_sq_ID = as.data.frame(Random_AgeStd_Slope_sq_ID)
Random_AgeStd_Slope_sq_ID$ID= rownames(Random_AgeStd_Slope_sq_ID)
# add information on population
Random_AgeStd_Slope_sq_ID=merge(Random_AgeStd_Slope_sq_ID,
                                unique(data_BBAL_atLeast4BA[,which(colnames(data_BBAL_atLeast4BA) %in% c("ID", "population"))]),
                                by = "ID",
                                all.x = TRUE, all.y = FALSE)

SlopeSqModel_ref=fixef(GLMM_brms_senescence_atLeast4BA)[3]
SlopeSqModel_Ker=fixef(GLMM_brms_senescence_atLeast4BA)[7]
# calculate the quadratic age effect for each individaul (fixed effect + random effect + fixed population effect)
Random_AgeStd_Slope_sq_ID$Estimate_corrected = Random_AgeStd_Slope_sq_ID$Estimate + fixef(GLMM_brms_senescence_atLeast4BA)[3]
Random_AgeStd_Slope_sq_ID$Estimate_corrected[which(Random_AgeStd_Slope_sq_ID$population =="Kerguelen" )] =   Random_AgeStd_Slope_sq_ID$Estimate_corrected[which(Random_AgeStd_Slope_sq_ID$population =="Kerguelen" )] + fixef(GLMM_brms_senescence_atLeast4BA)[7]


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
  
  df_AgeStd_randomSlopesAndInterceptsIDlevel$Year_last_seen[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] = unique(data_BBAL_atLeast4BA$Year_last_seen[which(data_BBAL_atLeast4BA$ID == ID)]) 
  df_AgeStd_randomSlopesAndInterceptsIDlevel$LRS[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] = unique(data_BBAL_atLeast4BA$lifetime_repro_output[which(data_BBAL_atLeast4BA$ID == ID)]) 
  df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] =  unique(data_BBAL_atLeast4BA$cohort[which(data_BBAL_atLeast4BA$ID == ID)]) 
  df_AgeStd_randomSlopesAndInterceptsIDlevel$cohortPop[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] =  unique(data_BBAL_atLeast4BA$cohortPop[which(data_BBAL_atLeast4BA$ID == ID)]) 
  df_AgeStd_randomSlopesAndInterceptsIDlevel$nb_BA_life[which( df_AgeStd_randomSlopesAndInterceptsIDlevel$ID == ID)] = unique(data_BBAL_atLeast4BA$nb_breeding_attempts[which(data_BBAL_atLeast4BA$ID == ID)]) 
  
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

# write.csv(df_AgeStd_randomSlopesAndInterceptsIDlevel,
#           paste0(dir_data, "data_BBAL_atLeast4BA_senescence_estimates.csv"))







#### __ calculate parameters of interest from posterior distribution #### 

#extract posterior values from the model
posterior_samples_finalModel <- posterior_samples(GLMM_brms_senescence_atLeast4BA)

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

# #### ____ Save output #### 
# write.csv(df_plotRes_Senescence,
#           "C:/Users/bmohring/OneDrive - The University of Liverpool/Current work folder/Code/Code_BBAL_senescence_paper/data_BBAL_atLeast4BA_senescence_posteriorEstimates.csv")




#### 3. Calculate variance in all parameters ####
data_BBAL_atLeast4BA$population=as.factor(data_BBAL_atLeast4BA$population)


#### ___ Year effect ####

# Extract the standard deviations for the random effects
sd_samples <- posterior_samples_finalModel[, grepl("sd_YearPop", names(posterior_samples_finalModel))]

# Calculate variances
variance_samples <- sd_samples^2

# Get the mean variances and 95% credible intervals
mean_variances <- colMeans(variance_samples)
ci_variances <- apply(variance_samples, 2, quantile, probs = c(0.025, 0.975))

# Combine results into a data frame
results <- data.frame(
  z_level = c("Bird island", "Kerguelen"),
  Variance = mean_variances,
  CI_lower = ci_variances[1, ],
  CI_upper = ci_variances[2, ]
)

print(results)


results$Variance[which(results$z_level =="Bird island")]/
  results$Variance[which(results$z_level =="Kerguelen")]  




#### ___ Cohort effect ####

# Extract the standard deviations for the random effects
sd_samples <- posterior_samples_finalModel[, grepl("sd_cohortPop", names(posterior_samples_finalModel))]

# Calculate variances
variance_samples <- sd_samples^2

# Get the mean variances and 95% credible intervals
mean_variances <- colMeans(variance_samples)
ci_variances <- apply(variance_samples, 2, quantile, probs = c(0.025, 0.975))

# Combine results into a data frame
results <- data.frame(
  z_level = c("Bird island", "Kerguelen"),
  Variance = mean_variances,
  CI_lower = ci_variances[1, ],
  CI_upper = ci_variances[2, ]
)

print(results)


results$Variance[which(results$z_level =="Bird island")]/
  results$Variance[which(results$z_level =="Kerguelen")]  




#### ___ random intercept ID effect ####

# Extract the standard deviations for the random effects
sd_samples <- posterior_samples_finalModel[, grepl("sd_ID__Intercept", names(posterior_samples_finalModel))]

# Calculate variances
variance_samples <- sd_samples^2

# Get the mean variances and 95% credible intervals
mean_variances <- colMeans(variance_samples)
ci_variances <- apply(variance_samples, 2, quantile, probs = c(0.025, 0.975))

# Combine results into a data frame
results <- data.frame(
  z_level = c("Bird island", "Kerguelen"),
  Variance = mean_variances,
  CI_lower = ci_variances[1, ],
  CI_upper = ci_variances[2, ]
)

print(results)


results$Variance[which(results$z_level =="Bird island")]/
  results$Variance[which(results$z_level =="Kerguelen")]  




#### ___ random age slope (linear) per ID effect ####

# Extract the standard deviations for the random effects
sd_samples <- posterior_samples_finalModel[, grepl("sd_ID__Age_std", names(posterior_samples_finalModel))]

# Calculate variances
variance_samples <- sd_samples^2

# Get the mean variances and 95% credible intervals
mean_variances <- colMeans(variance_samples)
ci_variances <- apply(variance_samples, 2, quantile, probs = c(0.025, 0.975))

# Combine results into a data frame
results <- data.frame(
  z_level = c("Bird island", "Kerguelen"),
  Variance = mean_variances,
  CI_lower = ci_variances[1, ],
  CI_upper = ci_variances[2, ]
)

print(results)


results$Variance[which(results$z_level =="Bird island")]/
  results$Variance[which(results$z_level =="Kerguelen")]  



#### ___ random age slope (quadratic) per ID effect ####

# Extract the standard deviations for the random effects
sd_samples <- posterior_samples_finalModel[, grepl("sd_ID__IAge_stdE2", names(posterior_samples_finalModel))]

# Calculate variances
variance_samples <- sd_samples^2

# Get the mean variances and 95% credible intervals
mean_variances <- colMeans(variance_samples)
ci_variances <- apply(variance_samples, 2, quantile, probs = c(0.025, 0.975))

# Combine results into a data frame
results <- data.frame(
  z_level = c("Bird island", "Kerguelen"),
  Variance = mean_variances,
  CI_lower = ci_variances[1, ],
  CI_upper = ci_variances[2, ]
)

print(results)


results$Variance[which(results$z_level =="Bird island")]/
  results$Variance[which(results$z_level =="Kerguelen")]  


#### ___ onset of senescence ####


var(df_plotRes_Senescence$OnsetSenesc[which(df_plotRes_Senescence$population=="Kerguelen")])
var(df_plotRes_Senescence$OnsetSenesc[which(df_plotRes_Senescence$population=="Bird island")])


var(df_plotRes_Senescence$OnsetSenesc[which(df_plotRes_Senescence$population=="Bird island")])/
  var(df_plotRes_Senescence$OnsetSenesc[which(df_plotRes_Senescence$population=="Kerguelen")]) 
 
#### ____ Var senescence rates #### 

var(df_plotRes_Senescence$SenesceRate[which(df_plotRes_Senescence$population=="Kerguelen")])
var(df_plotRes_Senescence$SenesceRate[which(df_plotRes_Senescence$population=="Bird island")])



var(df_plotRes_Senescence$SenesceRate[which(df_plotRes_Senescence$population=="Bird island")])/
  var(df_plotRes_Senescence$SenesceRate[which(df_plotRes_Senescence$population=="Kerguelen")]) 


 

#### ____ Var RS at 10 yo #### 

var(df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population=="Kerguelen")])
var(df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population=="Bird island")])


var(df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population=="Bird island")])/
  var(df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population=="Kerguelen")]) 

 

#### ____ Var RS at onset of senescence #### 

var(df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population=="Kerguelen")])
var(df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population=="Bird island")])


var(df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population=="Bird island")])/
  var(df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population=="Kerguelen")]) 

 

# Load packages
library(ggplot2)
library(ggdist)

 
#### PART 3: Figure 2 #### 
 

# Desired raw age breaks
raw_breaks <- c(5, 10, 15, 20, 25, 30, 35, 40, 45,50)
# Convert raw breaks to standardized scale
std_breaks <- (raw_breaks - mean_age_forStd) / sd_age_forStd

pA  = plot(conditional_effects(GLMM_brms_senescence_atLeast4BA), ask = FALSE)
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
  xlab("Early-life probability\nof successful reproduction")+
  theme_classic()+ theme(legend.position = "none",
                         axis.title.x = element_text(size = 10),  
                         axis.title.y = element_text(size = 10) )
# +
# ylab("Predicted age at onset of senescence")
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


# ggsave(paste0(dir_data, "Appendix1_Figure1_atleast4BA.png"))






#### ___ Overlap between distributions  #### 

library(overlapping)

# Onset of senescence
x_OnsetSenescence = list(Ker = df_plotRes_Senescence$OnsetSenesc_from_backTransformedCoefs[which(df_plotRes_Senescence$population =="Kerguelen")],
                         BI = df_plotRes_Senescence$OnsetSenesc_from_backTransformedCoefs[which(df_plotRes_Senescence$population =="Bird island")]) 
## bootstrapping 
out_OnsetSenescence  <- boot.overlap( x_OnsetSenescence, B = 100 ) 
out_OnsetSenescence$OVboot_stats 


# Rate of senescence
x_SenesceRate = list(Ker = df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs[which(df_plotRes_Senescence$population =="Kerguelen")],
                     BI=df_plotRes_Senescence$SenesceRate_from_backTransformedCoefs[which(df_plotRes_Senescence$population =="Bird island")]) 

## bootstrapping 
out_SenesceRate  <- boot.overlap( x_SenesceRate, B = 100 ) 
out_SenesceRate$OVboot_stats 

# PredictedRS_atpeak
x_PredictedRS_atpeak = list(Ker = df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population =="Kerguelen")],
                            BI=df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population =="Bird island")]) 

## bootstrapping 
out_PredictedRS_atpeak  <- boot.overlap( x_PredictedRS_atpeak, B = 100 ) 
out_PredictedRS_atpeak$OVboot_stats 



# PredictedRS_at10yo
x_PredictedRS_at10yo = list(Ker = df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population =="Kerguelen")],
                            BI=df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population =="Bird island")]) 

## bootstrapping 
out_PredictedRS_at10yo  <- boot.overlap( x_PredictedRS_at10yo, B = 100 ) 
out_PredictedRS_at10yo$OVboot_stats 





#### Figure 3: Relationship between individual predicted age at the onset of senescence and senescence rate  #### 


#### __ All data  #### 

df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs = -df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs

p = ggplot(df_AgeStd_randomSlopesAndInterceptsIDlevel,
           aes(x = SenesceRate_from_backTransformedCoefs, y =OnsetSenesc_from_backTransformedCoefs , 
               col = population, fill = population, group = population ))+
  geom_point(size = 1, alpha = 0.3)+
  scale_fill_manual(values = c("darkorange", "turquoise4")) + 
  scale_color_manual(values = c("darkorange", "turquoise4")) + 
  # scale_linetype_manual(values = c("darkorange", "turquoise4")) + 
  # geom_smooth(method= lm)+
  # geom_smooth(inherit.aes = FALSE,
  #             data  =df_rawAge_randomSlopesAndInterceptsIDlevel,
  #             aes(x = SenesceRate, y =OnsetSenesc_backTransformed),
  #             method= glm, method.args = list(family = "Gamma"))+
  theme_classic()+ theme(legend.position="none")+
  xlab("Senescence rate") + ylab("Age at onset of senescence")

pmarginal = ggExtra::ggMarginal(p, type = "density",  groupFill  = TRUE, groupColour = TRUE)

pmarginal
 


#### __ Outliers removed  #### 

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

# ggsave("C:/Users/bmohring/OneDrive - The University of Liverpool/Current work folder/Code/Code_BBAL_senescence_paper/Appendix_Figure3_atLeast4BA.png",
#        plot = pmarginal)

ggsave(paste0(dir_data, "Appendix1_Figure2_atleast4BA.png"),
       plot = pmarginal)


####  correlation between onset of senescence and senescence rate #### 

#### __ All data  #### 

# Kerguelen
cor(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Kerguelen")],
    df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Kerguelen")])
# -0.5810209

# Bird Island
cor(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Bird island")],
    df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Bird island")])
# -0.6162027



#### __ Outliers removed  #### 
cor(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Kerguelen" &
                                                                                             df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40 )],
    df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Kerguelen" &
                                                                                             df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40 )])
                                                                                             
cor(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Bird island" &
                                                                                             df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40 )],
    df_AgeStd_randomSlopesAndInterceptsIDlevel$SenesceRate_from_backTransformedCoefs[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$population =="Bird island" &
                                                                                             df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40 )])


#### PART 4: Link with LRS #### 

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


# Open dataset


df_AgeStd_randomSlopesAndInterceptsIDlevel = read.csv(paste0(dir_data, "data_BBAL_atLeast4BA_senescence_estimates.csv"))

df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort = substr(df_AgeStd_randomSlopesAndInterceptsIDlevel$cohortPop,1,4)
df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort =as.numeric(df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort )

# subset of dead individuals
SubsetDead_cohortBefore2000 = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$Year_last_seen <= 2020  &
                                                                                 df_AgeStd_randomSlopesAndInterceptsIDlevel$cohort <2000),]

table(SubsetDead_cohortBefore2000$population)


SubsetDead_cohortBefore2000_noOutlier  =SubsetDead_cohortBefore2000[which( SubsetDead_cohortBefore2000$OnsetSenesc_from_backTransformedCoefs<40 &
                                                                             SubsetDead_cohortBefore2000$OnsetSenesc_from_backTransformedCoefs >10),]

cor(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen")],
    SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen")])
cor(SubsetDead_cohortBefore2000_noOutlier$OnsetSenesc_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island")],
    SubsetDead_cohortBefore2000_noOutlier$SenesceRate_from_backTransformedCoefs[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island")])


#### Figure 4 #### 
library(scales)

Fig_4_A = ggplot(data = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Bird island"),], 
                 aes(x =SenesceRate_from_backTransformedCoefs , y=OnsetSenesc_from_backTransformedCoefs , col = LRS))+
  geom_point(size = 1.5, alpha = 0.6)+
  scale_colour_viridis_c(begin = 0.85, end = 0.15, option = "inferno", trans = pseudo_log_trans(base = 10), limits = c(0, 17))+
  theme_classic()+ xlim(c(0.005,0.06))+ ylim(c(16,32))+ 
  xlab("Senescence rate") + ylab("Age at onset of senescence")+
  theme( legend.text = element_text(size=8) , legend.key.size = unit(0.55, "cm"), legend.title = element_text(size=9)  )
# theme( legend.text = element_text(size=9) , legend.title = element_text(size=11),legend.position = c(0.95, 0.5) )
Fig_4_A

 
Fig_4_B = ggplot(data = SubsetDead_cohortBefore2000_noOutlier[which(SubsetDead_cohortBefore2000_noOutlier$population=="Kerguelen"  ),], 
                 aes(x =SenesceRate_from_backTransformedCoefs , y=OnsetSenesc_from_backTransformedCoefs , col = LRS))+
  geom_point(size = 1.5, alpha = 0.6)+
  scale_colour_viridis_c(begin = 0.90, end = 0.05, option = "mako", trans = pseudo_log_trans(base = 10))+
  theme_classic()+ xlim(c(0.005,0.06))+ ylim(c(15,32))+
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

pp_check(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence)

SubsetDead_cohortBefore2000_noOutlier$population=as.factor(SubsetDead_cohortBefore2000_noOutlier$population)
SubsetDead_cohortBefore2000_noOutlier$population_refKer = relevel(SubsetDead_cohortBefore2000_noOutlier$population , ref=2)

GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer  = brm(LRS   ~   OnsetSenescence_std * population_refKer  ,
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
summary(GLM_KerBI_Poisson_BRM_LRS_OnsetSenescence_refKer) 





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
                                             prior = prior_fitness,
                                             backend = "cmdstanr" ,
                                             control = list(adapt_delta = 0.85)
)

summary(GLM_KerBI_Poisson_BRM_LRS_SenesceRate) 


GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer  = brm(LRS   ~   SenesceRate_std * population_refKer  ,
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

summary(GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer) 

pp_check(GLM_KerBI_Poisson_BRM_LRS_SenesceRate)

Fig_4_D  = plot(conditional_effects(GLM_KerBI_Poisson_BRM_LRS_SenesceRate), ask = FALSE)
Fig_4_D = Fig_4_D$`SenesceRate_std:population`+  
  theme_classic() + 
  scale_fill_manual(values = c("darkorange", "turquoise4")) + 
  scale_color_manual(values = c("darkorange", "turquoise4")) + 
  geom_jitter(inherit.aes = FALSE, size =1.5,
              data = SubsetDead_cohortBefore2000_noOutlier, 
              width = 0, height = 0.2, alpha = 0.,
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

ggsave(paste0(dir_data, "Appendix_FigS2_3_AtLeast4BA.png"))





#### PART 5: Link with LRS dead or alive #### 



df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot = "Alive"
df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$Year_last_seen <= 2020  )]= "Dead"



table(df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot , df_AgeStd_randomSlopesAndInterceptsIDlevel$population)

df_AgeStd_randomSlopesAndInterceptsIDlevel$pop_DeadOrNot = paste0(df_AgeStd_randomSlopesAndInterceptsIDlevel$population,
                                                                  "_",
                                                                  df_AgeStd_randomSlopesAndInterceptsIDlevel$DeadOrNot)




Fig_4_C_v2 =  ggplot(data = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(   df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),] , 
                     aes(x= OnsetSenesc_from_backTransformedCoefs,  y  = LRS, group  = pop_DeadOrNot , col = pop_DeadOrNot ))+
  geom_jitter(inherit.aes = FALSE, size =1.5,
              data = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(   df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),] , 
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

Fig_4_D_v2 =  ggplot(data = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(  df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),] , 
                     aes(x= SenesceRate_from_backTransformedCoefs,  y  = LRS , group  = pop_DeadOrNot , col = pop_DeadOrNot ))+
  geom_jitter(inherit.aes = FALSE, size =1.5,
              data = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(  df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),], 
              width = 0, height = 0.2, alpha = 0.4,
              aes( x=  SenesceRate_from_backTransformedCoefs,  y  = LRS, group  = pop_DeadOrNot , col = pop_DeadOrNot ))+
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
 
 
ggsave(paste0(dir_data, "Appendix1_Figure4_atleast4BA.png"))








df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier = df_AgeStd_randomSlopesAndInterceptsIDlevel[which(df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs > 10 &
                                                                                                          df_AgeStd_randomSlopesAndInterceptsIDlevel$OnsetSenesc_from_backTransformedCoefs < 40   ),] 

df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenescence_std = robustHD::standardize(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$OnsetSenesc_from_backTransformedCoefs)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_std = robustHD::standardize(-df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$SenesceRate_from_backTransformedCoefs)
df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot = as.factor(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot)

table(df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$DeadOrNot, df_AgeStd_randomSlopesAndInterceptsIDlevel_noOutlier$population)

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

pp_check(GLM_KerBI_Poisson_BRM_LRS_SenesceRate_refKer_refDead)





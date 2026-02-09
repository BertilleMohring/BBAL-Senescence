#### Extract the variance in parameters of interest  ####
#### This code is used to produce Table 2  ####

# set the location of the data
dir_data_BBAL = "C:/Users/bmohring/Documents/GitHub/BBAL-Senescence/Data/"
dir_data = "C:/Users/bmohring/Documents/GitHub/BBAL-Senescence/"

# Open datasets
df_AgeStd_randomSlopesAndInterceptsIDlevel = read.csv( paste0(dir_data, "data_BBAL_senescence_estimates_priors.csv"))
df_plotRes_Senescence = read.csv(paste0(dir_data, "data_BBAL_senescence_posteriorEstimates_priors.csv"))

data_BBAL = read.csv(  paste0(dir_data_BBAL, "data_BBAL.csv"))
data_BBAL=data_BBAL[,-1]

# Load model output
load(paste0(dir_data, "model_output_GLMM_senescence_priors.RData"))


# calculate age mean and sd
mean_age_forStd =  mean(data_BBAL$Age) 
sd_age_forStd =  sd(data_BBAL$Age) 



#### 4. Calculate variance in all parameters ####
data_BBAL$population=as.factor(data_BBAL$population)


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




#### ___ random age slope (linear age effect) per ID effect ####

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



#### ___ random age slope (quadratic age effect) per ID effect ####

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


#### ___ Onset of senescence ####


var(df_plotRes_Senescence$OnsetSenesc[which(df_plotRes_Senescence$population=="Kerguelen")])
var(df_plotRes_Senescence$OnsetSenesc[which(df_plotRes_Senescence$population=="Bird island")])


var(df_plotRes_Senescence$OnsetSenesc[which(df_plotRes_Senescence$population=="Bird island")])/
  var(df_plotRes_Senescence$OnsetSenesc[which(df_plotRes_Senescence$population=="Kerguelen")]) 



#### ___ Senescence rate #### 

var(df_plotRes_Senescence$SenesceRate[which(df_plotRes_Senescence$population=="Kerguelen")])
var(df_plotRes_Senescence$SenesceRate[which(df_plotRes_Senescence$population=="Bird island")])

var(df_plotRes_Senescence$SenesceRate[which(df_plotRes_Senescence$population=="Bird island")])/
  var(df_plotRes_Senescence$SenesceRate[which(df_plotRes_Senescence$population=="Kerguelen")]) 




#### ___  Early-life probability of successful reproduction (Probability of successful reproduction at 10 years old) #### 

var(df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population=="Kerguelen")])
var(df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population=="Bird island")])


var(df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population=="Bird island")])/
  var(df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population=="Kerguelen")]) 

 
#### ___ Peak performance (Probability of successful reproduction at the onset of senescence) #### 

var(df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population=="Kerguelen")])
var(df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population=="Bird island")])


var(df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population=="Bird island")])/
  var(df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population=="Kerguelen")]) 

 


#### Calculate the overlap between the posterior distributions of the four life-history parameters for the two populations  #### 

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

# Probability of successful reproduction at the onset of senescence 
x_PredictedRS_atpeak = list(Ker = df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population =="Kerguelen")],
                            BI=df_plotRes_Senescence$PredictedRS_atpeak[which(df_plotRes_Senescence$population =="Bird island")]) 

## bootstrapping 
out_PredictedRS_atpeak  <- boot.overlap( x_PredictedRS_atpeak, B = 100 ) 
out_PredictedRS_atpeak$OVboot_stats 



# Early-life probability of successful reproduction 
x_PredictedRS_at10yo = list(Ker = df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population =="Kerguelen")],
                            BI=df_plotRes_Senescence$PredictedRS_at10yo[which(df_plotRes_Senescence$population =="Bird island")]) 

## bootstrapping 
out_PredictedRS_at10yo  <- boot.overlap( x_PredictedRS_at10yo, B = 100 ) 
out_PredictedRS_at10yo$OVboot_stats 




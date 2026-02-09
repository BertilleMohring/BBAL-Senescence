#### Run the main model for the paper ####

# Load packages

library(brms)
library(rstan)
library(StanHeaders)
library(bayesplot)
library(robustHD)
library(cmdstanr)

# set the location of the data
dir_data_BBAL = "C:/Users/bmohring/Documents/GitHub/BBAL-Senescence/Data/"

# Open dataset
data_BBAL = read.csv( paste0(dir_data_BBAL, "data_BBAL.csv"))
data_BBAL=data_BBAL[,-1]

# scale variables
data_BBAL$Age_std = robustHD::standardize(data_BBAL$Age)

mean_age_forStd =  mean(data_BBAL$Age) 
sd_age_forStd =  sd(data_BBAL$Age) 

# set priors for the model
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
GLMM_brms_senescence=brm(ReproS   ~   Age_std  + 
                           I(Age_std ^2)  +
                           First_breeding_attempt +  
                           population+
                           population:Age_std + 
                           population:I(Age_std ^2)  +    
                           (1 + Age_std +  I(Age_std ^2)|gr(ID, by = population)) +
                           (1| gr(YearPop, by = population))+
                           (1|gr(cohortPop, by = population)) ,
                         data=data_BBAL,
                         family = bernoulli(link = "logit"), 
                         warmup = 2000,
                         iter = 4000,
                         chains =8,
                         cores = 2, 
                         prior = prior_senescence,
                         backend = "cmdstanr"  ,
                         control = list(adapt_delta = 0.90)
)

summary(GLMM_brms_senescence)


# save output
save(GLMM_brms_senescence, 
     file = paste0(dir_data, "model_output_GLMM_senescence_priors.RData"))




# run model with the other population as a reference
# change reference population
data_BBAL$population=as.factor(data_BBAL$population)
data_BBAL$population_refKer = relevel(data_BBAL$population , ref=2)

# run model
GLMM_brms_senescence_refKer=brm(ReproS   ~   Age_std  + 
                           I(Age_std ^2)  +
                           First_breeding_attempt +  
                             population_refKer+
                             population_refKer:Age_std + 
                             population_refKer:I(Age_std ^2)  +    
                           (1 + Age_std +  I(Age_std ^2)|gr(ID, by = population_refKer)) +
                           (1| gr(YearPop, by = population_refKer))+
                           (1|gr(cohortPop, by = population_refKer)) ,
                         data=data_BBAL,
                         family = bernoulli(link = "logit"), 
                         warmup = 2000,
                         iter = 4000,
                         chains =8,
                         cores = 2,
                         prior = prior_senescence,
                         backend = "cmdstanr"  ,
                         control = list(adapt_delta = 0.90)
)

summary(GLMM_brms_senescence_refKer)


# save output
save(GLMM_brms_senescence_refKer, 
     file = paste0(dir_data, "model_output_GLMM_senescence_priors_refKer.RData"))

########################################
#           BEGINNER TUTORIAL          #
########################################


# Load Packages -----------------------------------------------------------

library(lme4)
library(lmerTest)

# Load Data ---------------------------------------------------------------

data <- read.csv('data/Steyversetal2019_subset_trim.csv')
head(data)

data$isswitch <- factor(data$isswitch)
data$compatible <- factor(data$compatible)

# Fit the null model ------------------------------------------------------

null_model <- lmer(response_time ~ 1 + (1|uid), data=data)

summary(null_model)

# Fit a one-predictor model

one_pred_model <- lmer(response_time ~ trial_num + (1|uid), data=data, REML = F)

summary(one_pred_model)


# Fit a one-predictor model with random slopes ----------------------------

one_pred_model_slopes <- lmer(response_time ~ trial_num + (trial_num|uid), data=data, REML = F)

summary(one_pred_model_slopes)

# Compare models with and without random slopes

anova(one_pred_model, one_pred_model_slopes)


# Fit a one-predictor with random slopes but no tau covariance -------------

one_pred_model_slopes_no_cov <- lmer(response_time ~ trial_num + (trial_num||uid), data=data, REML=F)

summary(one_pred_model_slopes_no_cov)

# Compare models with and without tau10

anova(one_pred_model_slopes, one_pred_model_slopes_no_cov)


# Fit a two-predictor model -----------------------------------------------

# center within-clusters

data$trial_num_WCC <- data$trial_num - ave(data$trial_num, data$uid, FUN=mean) 

# center between-clusters

data$trial_num_BCC <- data$trial_num - mean(data$trial_num)

# Model with centered variables 

two_pred_model <- lmer(response_time ~ isswitch + trial_num_WCC + (1|uid), data=data)

summary(two_pred_model)


# Fit a model with an level-1 interaction term ------------------------------------

interaction_model_l1 <- lmer(response_time ~ isswitch*trial_num_WCC + (1|uid), data=data)

summary(interaction_model)


# Fit a model with a cross-level interaction term -------------------------

interaction_model_cross <- lmer(response_time ~ agebin*trial_num_WCC + (1|uid), data=data)

summary(interaction_model_cross)


# Pseudo-R2 ---------------------------------------------------------------
# A reminder:
# null_model <- lmer(response_time ~ 1 + (1|uid), data=data)
# one_pred_model <- lmer(response_time ~ trial_num + (1|uid), data=data, REML = F)
# two_pred_model <- lmer(response_time ~ isswitch + trial_num_WCC + (1|uid), data=data)

null_model_sigma2 <-  as.data.frame(VarCorr(null_model))$vcov[2]
one_pred_model_sigma2 <- as.data.frame(VarCorr(one_pred_model))$vcov[2]
two_pred_model_sigma2 <- as.data.frame(VarCorr(two_pred_model))$vcov[2]

pseudoR2_null_one <- (null_model_sigma2-one_pred_model_sigma2)/null_model_sigma2
pseudoR2_null_two <- (null_model_sigma2-two_pred_model_sigma2)/null_model_sigma2
pseudoR2_one_two <- (one_pred_model_sigma2-two_pred_model_sigma2)/one_pred_model_sigma2

pseudoR2_null_one
pseudoR2_null_two
pseudoR2_one_two


# Useful Functions from Packages ------------------------------------------

library(sjPlot)

full_model <- lmer(response_time ~ agebin*isswitch + trial_num_WCC + trial_type + compatible + 
                     (1|uid), data=data)

tab_model(full_model)

plot_model(full_model, type='eff')

plot_model(full_model, type='int')


########################################
#             ADVANCED TUTORIAL        #
########################################

# Obtain profile confidence intervals -------------------------------------

# First, fit a model with random slopes

random_slopes_model <- lmer(response_time ~ compatible + isswitch + trial_num_WCC + 
                              (trial_num_WCC||uid), data=data)

summary(random_slopes_model)

# Obtain profile confidence intervals

profile_CI_for_random_slopes <- confint(random_slopes_model, oldNames=F)

profile_CI_for_random_slopes

# Empirical Bayes' estimates ----------------------------------------------

eb <- coef(random_slopes_model)$uid

head(eb)

ukj = ranef(random_slopes_model) # just deviation from gamma

head(ukj$uid)

lattice::dotplot(ukj)


# Categorical coding ------------------------------------------------------

# Dummy coding

contrasts(data$isswitch) = contr.treatment(length(unique(data$isswitch)))

contrasts(data$isswitch)

simple_cat_dummy <- lmer(response_time ~ isswitch + (1|uid), data=data)

summary(simple_cat_dummy)

eb <- coef(simple_cat_dummy)$uid

head(eb)

# Effects coding

contrasts(data$isswitch) <- contr.sum(length(unique(data$isswitch)))

simple_cat_effects <- lmer(response_time ~ isswitch + (1|uid), data=data)

summary(simple_cat_effects)

# Multiple categorical variables

contrasts(data$isswitch) <- contr.treatment(length(unique(data$isswitch)))
contrasts(data$trial_type) <- contr.treatment(length(unique(data$isswitch)))

complex_cat_effects <- lmer(response_time ~ isswitch + trial_type + (1|uid), data=data)

summary(complex_cat_effects)


# Fully-nested model ------------------------------------------------------

FN_model_null <- lmer(response_time ~ 1 + (1|uid/game_result_id), data=data)

summary(FN_model_null)

# Compute ICC 

tau2_t00 <- as.data.frame(VarCorr(FN_model_null))[1, 'vcov']
tau2_p00 <- as.data.frame(VarCorr(FN_model_null))[2, 'vcov']
sigma2 <-   as.data.frame(VarCorr(FN_model_null))[3, 'vcov']

ICC_game <- tau2_t00/sum(tau2_t00, tau2_p00, sigma2)
ICC_uid  <- tau2_p00/sum(tau2_t00, tau2_p00, sigma2)
ICC_res  <-sigma2/sum(tau2_t00, tau2_p00, sigma2)

ICC_game
ICC_uid
ICC_res


# Non-continuous outcomes -------------------------------------------------

binary_mlm_null <- glmer(accuracy ~ 1 + (1|uid), data=data, family='binomial')

summary(binary_mlm_null)

# convert G00 to probability
OR <- exp(3.0881)
pY <- OR/(1+OR)

eb <- coef(binary_mlm_null)$uid

# add column for probabilities

eb[,'pY'] = exp(eb$`(Intercept)`)/(1+exp(eb$`(Intercept)`)) 

head(eb)

# add a predictor to the model

contrasts(data$isswitch) <- contr.sum(length(unique(data$isswitch))) # effects coding

binary_mlm_one_pred <- glmer(accuracy ~ isswitch + (1|uid), data=data, family='binomial')

summary(binary_mlm_one_pred)

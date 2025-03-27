# Model code for CSTAT 06/27/2024
# M. Clark

# We are interested in the effect of Fgrm, a measure of genomic inbreeding, on the number of offspring measured from a 
# pedigree of rattlesnakes. We also think that there are likely impacts of sex and site on reproduction, as well as 
# kinship within the site. We are using principal components 2-6 to represent kinship structure. PC1 mostly captures
# the difference between sites and was dropped. Many individuals do not have offspring, so the data are zero-inflated. 
# We expect that some zeros are because the individual did not reproduce, and some zeros are because an individual's 
# offspring were not captured/part of the pedigree.

# load required libraries
library(ggeffects)
library(glmmTMB)
library(lmtest)
library(GGally)
library(effects) 
library(overdisp)
library(pscl)
library(emmeans)
library(MASS)

# load data object
load("../inbreeding_models/data_for_analyses_07132024.Robj", verbose = T) # contains data frame with relevant data ("data_flt")

# Filter data object
#------------------------------------------------------------------------------------------------------------------------------
dim(data) # 1037
plot(offspring~Fgrm, data = data, pch = 19, col = alpha("black", alpha = 0.25))

# remove individuals without sex information
data_flt <- data[-which(is.na(data$sex)),] # 1032
# convert sex and site to factors
data_flt$sex <- as.factor(data_flt$sex) 
data_flt$site <- as.factor(data_flt$site)
# remove individuals for whom we do not have year estimates
data_flt <- data_flt[-which(is.na(data_flt$yearsContribute2Ped)),] # removes 46 individuals
data_flt <- data_flt[-which(data_flt$yearsContribute2Ped == 0),] # removes 1 individual

# 977 individuals in data_flt, 240 from Barry County, 737 from Cass County 

# subset data
data_flt <- subset(data_flt, select = c(offspring, Fgrm, sex, site, PC1, PC2, PC3, PC4, PC5, PC6, yearsContribute2Ped))

plot(offspring~Fgrm, data = data_flt, pch = 19, col = alpha("black", alpha = 0.25))

# remove outliers in yearsContribute2Ped?
plot(offspring~yearsContribute2Ped, data = data_flt)
# data_flt <- subset(data_flt, subset = yearsContribute2Ped <= 16) # lets keep them in for now, can compare models and hopefully they won't 
# change results

head(data_flt)
#------------------------------------------------------------------------------------------------------------------------------

# visualize predictors
#------------------------------------------------------------------------------------------------------------------------------
ggpairs(data_flt)

# example correlations with Fgrm and PCs
plot(Fgrm~PC3, data = data_flt)
abline(lm(Fgrm~PC3, data = data_flt))

plot(Fgrm~PC5, data = data_flt)
abline(lm(Fgrm~PC5, data = data_flt))

plot(Fgrm~PC6, data = data_flt)
abline(lm(Fgrm~PC6, data = data_flt))

# We are primarily interested in the relationship between number of offspring and Fgrm, a measure of genomic inbreeding
plot(offspring~Fgrm, data = data_flt, pch = 19, col = alpha("black", alpha = 0.25))
plot(offspring~Fgrm, data = data_flt, pch = 19, col = data_flt$sex)
plot(offspring~Fgrm, data = data_flt[which(data_flt$sex == 1),], pch = 19, col = alpha("black", alpha = 0.25))
plot(offspring~Fgrm, data = data_flt[which(data_flt$sex == 2),], pch = 19, col = alpha("black", alpha = 0.25))

# I now also have estimates of the number of years each individual was alive in the population and could contribute 
# offspring to the pedigree. I am calling this "yearsContribute2Ped". I would like to use this as an offset in the models. 

hist(data_flt$yearsContribute2Ped)
plot(offspring~yearsContribute2Ped, data = data_flt, pch = 19, col = alpha("black", alpha = 0.25))

# years correlated with inbreeding?
plot(Fgrm~yearsContribute2Ped, data = data_flt, pch = 19, col = alpha("black", alpha = 0.25))
abline(lm(Fgrm~yearsContribute2Ped, data = data_flt))
summary(lm(Fgrm~yearsContribute2Ped, data = data_flt))

#------------------------------------------------------------------------------------------------------------------------------

# scale all predictors
#------------------------------------------------------------------------------------------------------------------------------
# standardized  covariates by subtracting each value by its mean and dividing it
# by the standard deviation so that continuous and binary regression coefficients were
# approximately on the same scale (Gelman 2008)
# Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped

# I am not centering predictors because Fgrm = 0 is meaningful 
data_scaled <- data.frame(matrix(data = NA, nrow = nrow(data_flt), ncol = 10))
colnames(data_scaled) <- c("offspring", "Fgrm", "sex", "site", "PC2", "PC3", "PC4", "PC5", "PC6", "yearsContribute2Ped")
data_scaled$offspring <- data_flt$offspring
data_scaled$Fgrm <- scale(data_flt$Fgrm, center = FALSE, scale = sd(data_flt$Fgrm))
data_scaled$sex <- as.factor(data_flt$sex)
data_scaled$site <- as.factor(data_flt$site)
data_scaled$PC2 <- scale(data_flt$PC2, center = FALSE, scale = sd(data_flt$PC2))
data_scaled$PC3 <- scale(data_flt$PC3, center = FALSE, scale = sd(data_flt$PC3))
data_scaled$PC4 <- scale(data_flt$PC4, center = FALSE, scale = sd(data_flt$PC4))
data_scaled$PC5 <- scale(data_flt$PC5, center = FALSE, scale = sd(data_flt$PC5))
data_scaled$PC6 <- scale(data_flt$PC6, center = FALSE, scale = sd(data_flt$PC6))
data_scaled$yearsContribute2Ped <- scale(data_flt$yearsContribute2Ped, center = FALSE, scale = sd(data_flt$yearsContribute2Ped))
# scale data the same as Eric to keep methods consistent
# data_scaled <- data.frame(matrix(data = NA, nrow = nrow(data_flt), ncol = 10))
# colnames(data_scaled) <- c("offspring", "Fgrm", "sex", "site", "PC2", "PC3", "PC4", "PC5", "PC6", "yearsContribute2Ped")
# data_scaled$offspring <- data_flt$offspring
# data_scaled$Fgrm <- (data_flt$Fgrm - mean(data_flt$Fgrm)) / (2*sd(data_flt$Fgrm))
# data_scaled$sex <- as.factor(data_flt$sex)
# data_scaled$site <- as.factor(data_flt$site)
# data_scaled$PC2 <- (data_flt$PC2 - mean(data_flt$PC2)) / (2*sd(data_flt$PC2))
# data_scaled$PC3 <- (data_flt$PC3 - mean(data_flt$PC3)) / (2*sd(data_flt$PC3))
# data_scaled$PC4 <- (data_flt$PC4 - mean(data_flt$PC4)) / (2*sd(data_flt$PC4))
# data_scaled$PC5 <- (data_flt$PC5 - mean(data_flt$PC5)) / (2*sd(data_flt$PC5))
# data_scaled$PC6 <- (data_flt$PC6 - mean(data_flt$PC6)) / (2*sd(data_flt$PC6))
# data_scaled$yearsContribute2Ped <- (data_flt$yearsContribute2Ped - mean(data_flt$yearsContribute2Ped)) / (2*sd(data_flt$yearsContribute2Ped))

#------------------------------------------------------------------------------------------------------------------------------

# Use modeling selection approach from Fávero et al 2021
#------------------------------------------------------------------------------------------------------------------------------

# Do we have overdispersion?
#------------------------------------------------------------------------------------------------------------------------------
overdisp(data_scaled, dependent.position = 1, predictor.position = c(2:10)) 

# yes, data is overdispersed 
mod1_pois <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                     data = data_scaled, 
                     family = poisson)
DHARMa::testDispersion(mod1_pois) # DHARMa test agrees 

# conclusion: use NB instead of POIS

#------------------------------------------------------------------------------------------------------------------------------

# is the data zero-inflated? 
#------------------------------------------------------------------------------------------------------------------------------
# how much of dataset is 0s? 
sum(data_scaled$offspring == 0)/length(data_scaled$offspring) # 77% 

# AIC comparison
mod1_NB <- glmmTMB(offspring~Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                   data = data_scaled, 
                   family = nbinom2)
plot(DHARMa::simulateResiduals(mod1_NB)) # look funky 

mod1_ZINB <- glmmTMB(offspring~Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                   ziformula = ~Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                   data = data_scaled, 
                   family = nbinom2)
plot(DHARMa::simulateResiduals(mod1_ZINB)) # look good 

mod1_NB1 <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                    data = data_scaled, 
                    family = nbinom1)
plot(DHARMa::simulateResiduals(mod1_NB1)) # look good 

mod1_ZINB1 <- glmmTMB(offspring ~ Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                    ziformula = ~Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                    data = data_scaled, 
                    family = nbinom1)
plot(DHARMa::simulateResiduals(mod1_NB1)) # look good 

mod_list <- list(mod1_NB, mod1_ZINB, mod1_NB1, mod1_ZINB1)
AICcmodavg::aictab(mod_list) # ZINB model most supported

# ZIP vs Poisson as suggested: 
# https://stats.stackexchange.com/questions/396336/r-glmm-for-unbalanced-zero-inflated-data-glmmtmb

mod1_P <- glmmTMB(offspring~Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                   data = data_scaled, 
                   family = poisson)
plot(DHARMa::simulateResiduals(mod1_P)) # bad 

mod1_ZIP <- glmmTMB(offspring~Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                     ziformula = ~Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                     data = data_scaled, 
                     family = poisson)
plot(DHARMa::simulateResiduals(mod1_ZIP)) # better, still funky 

mod_list <- list(mod1_P, mod1_ZIP)
AICcmodavg::aictab(mod_list) # ZIP model most supported

# Vuong test to compare NB to ZINB
# Keep code below for records, Vuong test is not approriate to test for zero inflation as reported by 
# Wilson 2015
# won't work with glmmTMB output, recreate models using glm.nb() and zeroinfl() 
# nb1 <- MASS::glm.nb(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, data=data_scaled)
# zinb <- zeroinfl(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, data = data_scaled, dist = "negbin", EM = TRUE)
#vuong(nb1, zinb)
# significant p-values (indicate ZINB is most appropriate) for raw and AIC-corrected, but not BIC-corrected

# what if I get rid of the interaction term?
# nb1 <- MASS::glm.nb(offspring ~ Fgrm + sex  + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, data=data_scaled)
# zinb <- zeroinfl(offspring ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, data = data_scaled, dist = "negbin", EM = TRUE)
# vuong(nb1, zinb)
# significant p-values (indicate ZINB is most appropriate) for raw and AIC-corrected, but not BIC-corrected

# Test with DHARMa 
DHARMa::testZeroInflation(mod1_NB) # p = 0.68

# from the DHARMa testZeroInflation doc: 

# Zero-inflation tests after fitting the model are crucial to see if you have zero-inflation. 
# Just because there are a lot of zeros doesn't mean you have zero-inflation, see Warton, D. I. (2005). 
# Many zeros does not mean zero inflation: comparing the goodness-of-fit of parametric models to 
# multivariate abundance data. Environmetrics 16(3), 275-289.

# That being said, zero-inflation tests are often not a reliable guide to decide wheter to add a zi term 
# or not. In general, model structures should be decided on ideally a priori, if that is not possible via 
# model selection techniques (AIC, BIC, WAIC, Bayes Factor). A zero-inflation test should only be run 
# after that decision, and to validate the decision that was taken.

# How does output from ZINB compare to NB? 
drop1(mod1_NB, test = "Chisq")
drop1(mod1_ZINB, test = "Chisq")

plot(allEffects(mod1_NB))
summary(mod1_NB)

plot(allEffects(mod1_ZINB))
summary(mod1_ZINB)

# compare predicted values
plot(x = predict(mod1_ZINB, type = "response", re.form = NA), 
     y = predict(mod1_NB, type = "response", re.form = NA), 
     xlab = "ZINB", 
     ylab = "NB")
abline(a = 0, b = 1)

# Does using ZINB instead of NB change results? 
ZINB_sigs <- drop1(mod1_ZINB, test = "Chisq")
NB_sigs <- drop1(mod1_NB, test = "Chisq")

plot(ZINB_sigs$`Pr(>Chi)`, NB_sigs$`Pr(>Chi)`)
abline(a = 0, b = 1)

comp_p_vals <- cbind.data.frame(rownames(NB_sigs)[-1], 
                 round(as.numeric(format(NB_sigs$`Pr(>Chi)`[-1]), scientific = F), digits = 4), 
                 round(as.numeric(format(ZINB_sigs$`Pr(>Chi)`[-1]), scientific = F), digits = 4))
colnames(comp_p_vals) <- c("effect", "NB", "ZINB")

# what about effect sizes? 
ZINB_beta <- summary(mod1_ZINB)$coefficients$cond
NB_beta <- summary(mod1_NB)$coefficients$cond

cbind.data.frame(ZINB_beta[,"Estimate"], NB_beta[,"Estimate"])
# direction and general magnitude of effect sizes in the conditional model stay the same, except for PC2

#conclusion: most lines of evidence, including our a priori expectations indicate that we should use a ZINB model
#------------------------------------------------------------------------------------------------------------------------------

# Model building and comparison
#------------------------------------------------------------------------------------------------------------------------------

# I have three models built to test specific hypotheses: 
# mod1: offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6
    # is there an effect of Fgrm and Fgrm*sex on number of offspring?
# mod2: offspring ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6
    # is there an effect of Fgrm on number of offspring?
# mod3: offspring ~ sex + site + PC2 + PC3 + PC4 + PC5 + PC6
    # number of offspring is well explained by a combination of sex, site, and kinship (approximated using principal components)
    # I initially wanted to use a genetic relatedness matrix instead of PCs, but I couldn't figure out how to set a custom variance matrix 
    # using glmmTMB()

#------------------------------------------------------------------------------------------------------------------------------

# yearsContribute2Ped as an offset term?
#------------------------------------------------------------------------------------------------------------------------------

# Should I treat yearsContribute2Ped as an offset? 

mod1_compois <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + offset(log(yearsContribute2Ped)), 
                ziformula = ~., 
                data = data_scaled, 
                family = compois)

mod1_compois_nooff <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                        ziformula = ~., 
                        data = data_scaled, 
                        family = compois)

mod1_pois <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + offset(log(yearsContribute2Ped)), 
                ziformula = ~., 
                data = data_scaled, 
                family = poisson)

mod1_pois_noof <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                     ziformula = ~., 
                     data = data_scaled, 
                     family = poisson)

mod1_negbinom <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + offset(log(yearsContribute2Ped)), 
                ziformula = ~., 
                data = data_flt, 
                family = nbinom2)
mod1_negbinom_nooff <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                         ziformula = ~., 
                         data = data_flt, 
                         family = nbinom2)


plot(DHARMa::simulateResiduals(mod1_compois)) # not great
plot(DHARMa::simulateResiduals(mod1_pois)) # not great
plot(DHARMa::simulateResiduals(mod1_negbinom)) # not great

plot(DHARMa::simulateResiduals(mod1_compois_nooff)) # looks good
plot(DHARMa::simulateResiduals(mod1_pois_noof)) # not great
plot(DHARMa::simulateResiduals(mod1_negbinom_nooff)) # looks good

# using AIC to compare the models with different distributions 
mod_list <- list(mod1_compois, mod1_pois, mod1_negbinom, mod1_compois_nooff, mod1_pois_noof, mod1_negbinom_nooff) # , modnames = c("compois", "pois", "negbinom"))
AICcmodavg::aictab(mod_list) 

# The negative binomial model with no offset is the best fit to the data. 

# compare predicted values from different modes: 
plot(x = predict(mod1_negbinom_nooff, type = "response", re.form = NA), 
     y = predict(mod1_negbinom, type = "response", re.form = NA), 
     xlab = "negbinom, no off", 
     ylab = "negbinom, off")

plot(x = predict(mod1_negbinom_nooff, type = "response", re.form = NA), 
     y = predict(mod1_compois_nooff, type = "response", re.form = NA), 
     xlab = "negbinom, no off", 
     ylab = "compois, off")

plot(x = predict(mod1_negbinom_nooff, type = "response", re.form = NA), 
     y = predict(mod1_pois_noof, type = "response", re.form = NA), 
     xlab = "negbinom, no off", 
     ylab = "pois, off")

plot(y = round(predict(mod1_negbinom_nooff, type = "response", re.form = NA)), x = data_scaled$Fgrm, 
     xlab = "Fgrm", ylab = "predicted offspring")
plot(offspring~Fgrm, data = data_scaled, xlab = "Fgrm", ylab = "offspring")

#------------------------------------------------------------------------------------------------------------------------------

# build three models using ZINB and no offset: 
#------------------------------------------------------------------------------------------------------------------------------

mod1 <- glmmTMB(offspring ~ Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                ziformula = ~., 
                data = data_scaled, 
                family = nbinom2)
# check model output
plot(DHARMa::simulateResiduals(mod1))
diagnose(mod1) # Unusually large Z-statistics
DHARMa::testDispersion(mod1)

# look at important effects
drop1(mod1, test = "Chisq") # p-value tells you term is explaining significant variation in the model
summary(mod1) # p-values not super informative, but effect sizes are (but they are scaled)
plot(allEffects(mod1)) # just the conditional model 
emm_mod1_zi <- emmeans(mod1, ~ Fgrm*sex, at = list(Fgrm = round(seq(-2.4, 6, length.out = 10), digits = 2)), component = "zi")
plot(emm_mod1_zi, horizontal = FALSE, 
     xlab = "emmean (predicted offspring on the log scale), conditional model", 
     ylab = "scaled Fgrm", 
     comparisons = FALSE)

mod2 <- glmmTMB(offspring ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                ziformula = ~., 
                data = data_scaled, 
                family = nbinom2)

# check model output
plot(DHARMa::simulateResiduals(mod2))
diagnose(mod2) # Unusually large Z-statistics
DHARMa::testDispersion(mod2)

# look at important effects
drop1(mod2, test = "Chisq") # p-value tells you term is explaining significant variation in the model
summary(mod2) # p-values not super informative, but effect sizes are (but they are scaled)
plot(allEffects(mod2)) 
emm_mod2_zi <- emmeans(mod2, ~ Fgrm, at = list(Fgrm = round(seq(-2.4, 6, length.out = 10), digits = 2)), component = "zi")
plot(emm_mod2_zi, horizontal = FALSE, 
     xlab = "emmean (predicted offspring on the log scale), conditional model", 
     ylab = "scaled Fgrm", 
     comparisons = FALSE)

mod3 <- glmmTMB(offspring ~ sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                ziformula = ~., 
                data = data_scaled, 
                family = nbinom2)

# check model output
plot(DHARMa::simulateResiduals(mod3))
diagnose(mod3) # Unusually large Z-statistics
DHARMa::testDispersion(mod3)

# look at important effects
drop1(mod3, test = "Chisq") # p-value tells you term is explaining significant variation in the model
summary(mod3) # p-values not super informative, but effect sizes are (but they are scaled)
plot(allEffects(mod3)) 

# use AICc to evaluate the best supported model: 
mod_list <- list(mod1, mod2, mod3)
AICcmodavg::aictab(mod_list) # mod2, mod1, mod3

# average models? it doesn't make sense to me to average an interaction term 
# MuMIn::model.avg(mod1, mod2)

#------------------------------------------------------------------------------------------------------------------------------

# Evaluate impact of outliers 
#------------------------------------------------------------------------------------------------------------------------------
# Four individuals with large number of offspring because including them results in 
# extremely high odds ratios in the zero-inflated model

data_noHighOff <- subset(data_scaled, offspring < 15) # dropping 4 individuals 

plot(offspring~Fgrm, data = data_noHighOff, pch = 19, col = alpha("black", alpha = 0.25))

mod1_noHighOff <- glmmTMB(offspring ~ Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, ziformula = ~., 
                          data = data_noHighOff, family = nbinom2)
mod2_noHighOff <- glmmTMB(offspring ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, ziformula = ~., 
                          data = data_noHighOff, family = nbinom2)
mod3_noHighOff <- glmmTMB(offspring ~ sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, ziformula = ~., 
                          data = data_noHighOff, family = nbinom2)

mod_list_noHighOff <- list(mod1_noHighOff, mod2_noHighOff, mod3_noHighOff)
AICcmodavg::aictab(mod_list_noHighOff) # mod2_noHighOff, mod1_noHighOff, mod3_noHighOff, does not change model support 

# Two individuals with very long years contributing

data_noHighCont <- subset(data_scaled, yearsContribute2Ped < 6) # scaled data! 
par(mfrow=c(2,2))
plot(offspring~yearsContribute2Ped, data = data_scaled, pch = 19, col = alpha("black", alpha = 0.25))
plot(offspring~yearsContribute2Ped, data = data_noHighCont, pch = 19, col = alpha("black", alpha = 0.25))

plot(offspring~Fgrm, data = data_scaled, pch = 19, col = alpha("black", alpha = 0.25))
plot(offspring~Fgrm, data = data_noHighCont, pch = 19, col = alpha("black", alpha = 0.25))

mod1_noHighCont <- glmmTMB(offspring ~ Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, ziformula = ~., 
                          data = data_noHighCont, family = nbinom2)
mod2_noHighCont <- glmmTMB(offspring ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, ziformula = ~., 
                          data = data_noHighCont, family = nbinom2)
mod3_noHighCont <- glmmTMB(offspring ~ sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, ziformula = ~., 
                          data = data_noHighCont, family = nbinom2)

mod_list_noHighCont <- list(mod1_noHighCont, mod2_noHighCont, mod3_noHighCont)
AICcmodavg::aictab(mod_list_noHighCont) # mod1_noHighCont, mod2_noHighCont, mod3_noHighCont, does change model support! 

summary(mod1)
summary(mod1_noHighCont)
plot(allEffects(mod1))
plot(allEffects(mod1_noHighCont))

# dropping two individuals with high longevity (we’re talking 19 and 25 YEARS), did make model 1 the 
# top model over model 2 (again, very small delta AICc between the models). These two individuals 
# happen to have low inbreeding, but 0 and 2 offspring respectively. I had initially dropped these 
# individuals in previous runs because both were captured as large adults, which increases error around  
# birth year estimates, and the longevity seemed unrealistic. Eric’s response was “the two… are a bit out 
# there, but PCC does have higher survival rates so it isn’t entirely implausible that a couple 
# individuals could live to 20 and 25". The effect of “yearsContribute2Ped” is very significant with 
# similar effect sizes regardless of whether these two individual are included. I’m in favor of keeping 
# them I think. We have no reason to doubt their measure of Fgrm or offspring, and it’s likely that there 
# is error around many of the estimates of “yearsContribute2Ped” that we can’t detect, plus the fact that 
# one of the emr experts didn’t think it was totally nuts to have individuals that old.

# Explore models 
#------------------------------------------------------------------------------------------------------------------------------

# Do I see similar model AIC comparison with the NB instead of ZINB? 
mod1_NB <- glmmTMB(offspring ~ Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                data = data_scaled, 
                family = nbinom2)
mod2_NB <- glmmTMB(offspring ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                data = data_scaled, 
                family = nbinom2)
mod3_NB <- glmmTMB(offspring ~ sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                   data = data_scaled, 
                   family = nbinom2)
mod_list_NB <- list(mod1_NB, mod2_NB, mod3_NB)
AICcmodavg::aictab(mod_list) # mod2_NB, mod1_NB, mod3_NB

# including a zero inflation term make the interaction between Fgrm and sex more important to the fit of the model, though the 
# order of models is the same 

#------------------------------------------------------------------------------------------------------------------------------
# Visualize 0s and non-0s separately
#------------------------------------------------------------------------------------------------------------------------------
data_noZero <- subset(data_scaled, offspring > 0)
plot(offspring~Fgrm, data_noZero, pch = 19, col = alpha("black", alpha = 0.25))
abline(lm(offspring~Fgrm, data_noZero))

mod2_noZero <- glmmTMB(offspring ~ Fgrm + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                       data = data_noZero, 
                       family = nbinom2)
drop1(mod2_noZero, test = "Chisq")
summary(mod2_noZero)

# probability of having 0 offspring
# number of individuals with Fgrm > 0 with no offspring divided by the number of individuals wiht Fgrm > 0
sum(subset(data_scaled, Fgrm > 0, select = offspring) == 0) / nrow(subset(data_scaled, Fgrm > 0, select = offspring))
# 80% of individuals with Fgrm > 0 have no offspring 

sum(subset(data_scaled, Fgrm <= 0, select = offspring) == 0) / nrow(subset(data_scaled, Fgrm <= 0, select = offspring))
# 73% of individuals with Fgrm <= 0 have no offspring 

# what is the average number of offspring for different catagories of inbred? 
mean(subset(data_scaled, select = offspring)$offspring) # overall: 0.7635619
mean(subset(data_scaled, Fgrm > 0, select = offspring)$offspring) # inbred: 0.5875
mean(subset(data_scaled, Fgrm <= 0, select = offspring)$offspring) # not inbred: 0.9336016

# males
mean(subset(data_scaled, sex == 2, select = offspring)$offspring) # overall: 0.682243
mean(subset(data_scaled, Fgrm > 0 & sex == 2, select = offspring)$offspring) # inbred: 0.6165803
mean(subset(data_scaled, Fgrm <= 0 & sex ==2, select = offspring)$offspring) # not inbred: 0.7361702

# females
mean(subset(data_flt, sex == 1, select = offspring)$offspring) # overall: 0.8269581
mean(subset(data_flt, Fgrm > 0 & sex == 1, select = offspring)$offspring) # inbred: 0.5679443
mean(subset(data_flt, Fgrm <= 0 & sex == 1, select = offspring)$offspring) # not inbred: 1.110687
# this looks like there is a higher impact of inbreeding on number of offspring in females, matches cond model with interaction term (mod1)

high_Fgrm_off <- subset(data_scaled, Fgrm > 0, select = offspring)$offspring
low_Fgrm_off <- subset(data_scaled, Fgrm <= 0, select = offspring)$offspring

dat <- cbind.data.frame(c(rep("high_Fgrm", times = length(high_Fgrm_off)), rep("low_Fgrm", times = length(low_Fgrm_off))), 
                 c(high_Fgrm_off, low_Fgrm_off))
colnames(dat) <- c("level", "offspring")

vioplot::vioplot(offspring~level, data = dat)
boxplot(offspring~level, data = dat)



#------------------------------------------------------------------------------------------------------------------------------

# modeling sexes separately 
#------------------------------------------------------------------------------------------------------------------------------
# set up binary offspring variable
data_scaled$offspring_binary <- NA
data_scaled[which(data_scaled$offspring > 0), "offspring_binary"] <- 1 
data_scaled[which(data_scaled$offspring == 0), "offspring_binary"] <- 0

plot(allEffects(mod1)) # mod1 shows a neg trend between offspring and Fgrm for males, but not for females
# but look at average number of offspring shows a greater impact on females? So what is happening? 

dim(subset(data_flt, sex ==1)) # 549
dim(subset(data_flt, sex ==1 & offspring == 0)) # 413, 75% have no offspring

dim(subset(data_flt, sex ==2)) # 428
dim(subset(data_flt, sex ==2 & offspring == 0)) # 339, 79% have no offspring 

# relationship in females is weaker than males, females are washing out signal w/o interaction? 

data_female <- subset(data_scaled, sex == 1)
plot(offspring~Fgrm, data = data_female, col = alpha("black", alpha = 0.25), pch = 19, main = "females")
data_male <- subset(data_scaled, sex == 2)
plot(offspring~Fgrm, data = data_male, col = alpha("black", alpha = 0.25), pch = 19, main = "males")

mod2_female <- glmmTMB(offspring ~ Fgrm + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                ziformula = ~Fgrm + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                data = data_female, 
                family = nbinom2)
drop1(mod2_female, test = "Chisq") # site, PC6, and years sig.
summary(mod2_female)

mod3_female <- glmmTMB(offspring ~ site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                       ziformula = ~site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                       data = data_female, 
                       family = nbinom2)

mod_list <- list(mod2_female, mod3_female)
AICcmodavg::aictab(mod_list) # mod2_female, mod3_female

# models don't converge... 
mod2_male <- glmmTMB(offspring ~ Fgrm + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                       ziformula = ~Fgrm + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                       data = data_male, 
                       family = nbinom2)
drop1(mod2_male, test = "Chisq") # PC2 and years sig.
summary(mod2_male)

mod3_male <- glmmTMB(offspring ~ site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                     ziformula = ~site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                     data = data_male, 
                     family = nbinom2)
mod_list <- list(mod2_male, mod3_male)
AICcmodavg::aictab(mod_list) # mod3_male, mod2_male, very similar AICc

mean(subset(data_female, offspring > 0, select = offspring)$offspring) # 3.338235
mean(subset(data_male, offspring > 0, select = offspring)$offspring) # 3.321839

data_female_noZero <- subset(data_female, offspring > 0)
data_male_noZero <- subset(data_male, offspring > 0)

plot(offspring~Fgrm, data = data_female_noZero, pch = 19, 
     col = alpha("black", alpha = 0.25), 
     main = "female")
plot(offspring~Fgrm, data = data_male_noZero, pch = 19, 
     col = alpha("black", alpha = 0.25), 
     main = "male")

# clearer trend in offspring counts in males than females (what we are seeing in plot(allEffects))

# all
plot(offspring_binary~Fgrm, data = data_scaled)
mod_logit <- glmmTMB(offspring_binary~Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, data = data_scaled, family = "binomial")
plot(DHARMa::simulateResiduals(mod_logit)) # look good
diagnose(mod_logit) # Unusually large Z-statistics

summary(mod_logit)
drop1(mod_logit, test = "Chisq") # Fgrm is marginally significant (0.055)
plot(allEffects(mod_logit), type = "response")

summary(mod_logit)
1-(exp(-0.62755 / sd(data_flt$Fgrm))^0.1) # females, a 0.1 increase in Fgrm is associated with a 78.74% decrease in the odds of having offspring
1-(exp((-0.62755 + 0.45472) / sd(data_flt$Fgrm))^0.1) # males, a 0.1 increase in Fgrm is associated with a 34.72% decrease in the odds of having offspring

mod2_logit <- glmmTMB(offspring_binary~Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, data = data_scaled, family = "binomial")
plot(allEffects(mod2_logit), type = "response")
summary(mod2_logit)
drop1(mod2_logit, test = "Chisq") # Fgrm is significant

1-(exp(-0.41065 / sd(data_flt$Fgrm))^0.1) # a 0.1 increase in Fgrm is associated with a 64.0% decrease in the odds of having offspring 

# females
plot(offspring_binary~Fgrm, data = data_female)
mod_logit_females <- glmmTMB(offspring_binary~Fgrm + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, data = data_female, family = "binomial")
plot(DHARMa::simulateResiduals(mod_logit_females)) # look okay
diagnose(mod_logit_females) # Unusually large Z-statistics

summary(mod_logit_females)
drop1(mod_logit_females, test = "Chisq") # Fgrm is significant 
plot(allEffects(mod_logit_females), type = "response")

# males
plot(offspring_binary~Fgrm, data = data_male)
mod_logit_males <- glmmTMB(offspring_binary~Fgrm + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, data = data_male, family = "binomial")
plot(DHARMa::simulateResiduals(mod_logit_males)) # look okay... not great
diagnose(mod_logit_males) # Unusually large Z-statistics

summary(mod_logit_males)
drop1(mod_logit_males, test = "Chisq") 
plot(allEffects(mod_logit_males), type = "response")

# My interpretation is that inbred females have trouble just having offspring, 
# but those who do are able to have approx. as many as non-inbred females. 
# Inbred males tend to have offspring at similar rates to noninbred males, 
# but do not have as many offspring
#------------------------------------------------------------------------------------------------------------------------------


# Effect sizes for most supported model 
#------------------------------------------------------------------------------------------------------------------------------

mod2_fixef <- fixef(mod2) # cond.Fgrm: -0.101503, zi: 0.34301

# convert to normal scale
mod2_fixef$cond["Fgrm"] / sd(data_flt$Fgrm) # -2.53204
mod2_fixef$zi["Fgrm"] / sd(data_flt$Fgrm) # 8.537613

emmeans(mod2, ~ yearsContribute2Ped, component = "cond", at = list(yearsContribute2Ped = range(data_scaled$yearsContribute2Ped)), type = "response") # 1.43
emmeans(mod2, ~ Fgrm, component = "cond", at = list(Fgrm = range(data_scaled$Fgrm)[1]), type = "response") # 1.43
emmeans(mod2, ~ Fgrm, component = "cond", at = list(Fgrm = range(data_scaled$Fgrm)[2]), type = "response") # 0.508

range(data_scaled$Fgrm)[1] * sd(data_flt$Fgrm) #  -0.0946
range(data_scaled$Fgrm)[2] * sd(data_flt$Fgrm) #  0.315

# As Fgrm changes from -0.09459466 to 0.3148291, the predicted number of offspring for an average individual changes from 1.43 to 0.508

mod2_fixef$zi["Fgrm"] / sd(data_flt$Fgrm)
exp(mod2_fixef$zi["Fgrm"] / sd(data_flt$Fgrm)) # odds ratio (quite large) for 1 unit increase
exp(mod2_fixef$zi["Fgrm"] / sd(data_flt$Fgrm))^0.01 # odds ratio (quite large) for 0.01 unit increase
((exp(mod2_fixef$zi["Fgrm"] / sd(data_flt$Fgrm))^0.01)-1) * 100 # an 0.01 increase in Fgrm is associated with an 8% decrease in the odds of having offspring
((exp(mod2_fixef$zi["Fgrm"] / sd(data_flt$Fgrm))^0.1)-1) * 100 # an 0.1 increase in Fgrm is associated with an 134.8% decrease in the odds of having offspring
((exp(mod2_fixef$zi["Fgrm"] / sd(data_flt$Fgrm))^0.05)-1) * 100 # an 0.05 increase in Fgrm is associated with an 53.25% decrease in the odds of having offspring

((exp(mod2_fixef$zi["Fgrm"] / sd(data_flt$Fgrm))^range(data_flt$Fgrm)[2] - range(data_flt$Fgrm)[1])-1) * 100 # an 0.05 increase in Fgrm is associated with an 53.25% decrease in the odds of having offspring


emmeans(mod2, ~ Fgrm, component = "zi", at = list(Fgrm = range(data_scaled$Fgrm)))

# plot 

# make confidence estimates: 

mod2_conf <- as.data.frame(confint(mod2))
rownames(mod2_conf) <- c("Cond(Intercept)", "Cond(Fgrm)", "Cond(sex)", "Cond(site)", "Cond(PC2)", "Cond(PC3)", "Cond(PC4)", "Cond(PC5)", "Cond(PC6)", "Cond(years)", "Zi(Intercept)", "Zi(Fgrm)", "Zi(sex)", "Zi(site)", "Zi(PC2)", "Zi(PC3)", "Zi(PC4)", "Zi(PC5)", "Zi(PC6)", "Zi(years)")
mod2_conf_85 <- as.data.frame(confint(mod2, level = 0.85))
mod2_conf <- cbind.data.frame(rownames(mod2_conf), mod2_conf, mod2_conf_85$`7.5 %`, mod2_conf_85$`92.5 %`)
colnames(mod2_conf) <- c("factor", "LCI_95", "UCI_95", "estimate", "LCI_85", "UCI_85")
mod2_conf <- mod2_conf[-c(1, 11),] # remove estimates for intercepts 

# Betas figure
mod2_betas <-  ggplot(data =  mod2_conf) + geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 9.5, color = "black", linetype = 2) +
  geom_linerange(aes(x = factor, 
                     y = estimate, ymin = LCI_85, ymax = UCI_85),linewidth=1.25) + 
  geom_pointrange(aes(x = factor, 
                      y = estimate, ymin = LCI_95, ymax = UCI_95),linewidth=.5) + 
  coord_flip() +
  labs(x = NULL, y = NULL)+ theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(color = "black", size=14)) +
  theme(axis.text.y = element_text(color = "black", size=14)) +
  scale_y_continuous(breaks=seq(-4, 4, by=1), limits=c(-4, 4))+
  scale_x_discrete(labels=factor)

mod2_betas

# mod 1
mod1_conf <- as.data.frame(confint(mod1))
rownames(mod1_conf) <- c("Cond(Intercept)", "Cond(Fgrm)", "Cond(sex)", "Cond(site)", "Cond(PC2)", "Cond(PC3)", "Cond(PC4)", "Cond(PC5)", "Cond(PC6)", "Cond(years)", "Cond(Fgrm*sex)", "Zi(Intercept)", "Zi(Fgrm)", "Zi(sex)", "Zi(site)", "Zi(PC2)", "Zi(PC3)", "Zi(PC4)", "Zi(PC5)", "Zi(PC6)", "Zi(years)", "Zi(Fgrm*sex)")
mod1_conf_85 <- as.data.frame(confint(mod1, level = 0.85))
mod1_conf <- cbind.data.frame(rownames(mod1_conf), mod1_conf, mod1_conf_85$`7.5 %`, mod1_conf_85$`92.5 %`)
colnames(mod1_conf) <- c("factor", "LCI_95", "UCI_95", "estimate", "LCI_85", "UCI_85")
mod1_conf <- mod1_conf[-c(1, 12),] # remove estimates for intercepts 

# Betas figure
mod1_betas <-  ggplot(data =  mod1_conf) + geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 9.5, color = "black", linetype = 2) +
  geom_linerange(aes(x = factor, 
                     y = estimate, ymin = LCI_85, ymax = UCI_85),linewidth=1.25) + 
  geom_pointrange(aes(x = factor, 
                      y = estimate, ymin = LCI_95, ymax = UCI_95),linewidth=.5) + 
  coord_flip() +
  labs(x = NULL, y = NULL)+ theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(color = "black", size=14)) +
  theme(axis.text.y = element_text(color = "black", size=14)) +
  scale_y_continuous(breaks=seq(-4, 4, by=1), limits=c(-3.5, 3.5))+
  scale_x_discrete(labels=factor)

mod1_betas


tiff(file="../figures/SupplementalFig12.tiff", 
     width=7, height=5, units="in", res=300)
mod2_betas
dev.off()

tiff(file="../figures/SupplementalFig13.tiff", 
     width=7, height=5, units="in", res=300)
mod1_betas
dev.off()


# effect size figure 

# one way
mod2_zprob <- predict(mod2, type = "zprob")
plot(mod2_zprob~Fgrm, data = data_flt)

# Effect plot for the zero-inflated portion
eff <- effect("Fgrm", mod2, which = "zero")
e <- allEffects(mod = mod2, xlevels = 200, which = "zero", type = "response") %>% as.data.frame() 

plot(allEffects(mod = mod2, xlevels = 200, which = "zero")) # what is the y axis here? 

Fgrm_effects_zi <- e[["Fgrm"]]
Fgrm_effects_zi$Fgrm <- Fgrm_effects_zi$Fgrm * sd(data_flt$Fgrm)

plot(fit~Fgrm, data = Fgrm_effects_zi, type = "l", col = "black", 
     xlab = "Fgrm", ylab = "Offspring", 
     cex.axis = 1.2, 
     cex.lab = 1.75)
polygon(x=c(Fgrm_effects_zi$Fgrm, rev(Fgrm_effects_zi$Fgrm)), y = c(Fgrm_effects_zi$upper, rev(Fgrm_effects_zi$lower)), 
        col = alpha("gray", alpha = 0.5), border = NA)

plot(prob~Fgrm, data = Fgrm_effects_zi, type = "l", col = "black", 
     xlab = "Fgrm", ylab = "Probability", 
     cex.axis = 1.2, 
     cex.lab = 1.75, 
     xlim = range(data_flt$Fgrm), 
     ylim = c(0,1))
polygon(x=c(Fgrm_effects_zi$Fgrm, rev(Fgrm_effects_zi$Fgrm)), y = c(Fgrm_effects_zi$upper_prob, rev(Fgrm_effects_zi$lower_prob)), 
        col = alpha("gray", alpha = 0.5), border = NA)

emm_results_zi <- emmeans(mod2, ~ Fgrm, at = list(Fgrm = round(seq(from = range(data_scaled$Fgrm)[1], to =  range(data_scaled$Fgrm)[2], length.out = 200), digits = 2)), component = "zi", type = "response")
plot(emm_results_zi, horizontal = FALSE, 
     xlab = "emmean (predicted offspring on the log scale), conditional model", 
     ylab = "scaled Fgrm", 
     comparisons = FALSE)
emm_results_zi_df <- as.data.frame(emm_results_zi)
emm_results_zi_df$orig_Fgrm <- emm_results_zi_df$Fgrm * sd(data_flt$Fgrm)
# create offspring binary
offspring_binary <- data_flt$offspring
offspring_binary[which(offspring_binary == 0)] <- NA # 0s become NAs (temp)
offspring_binary[which(offspring_binary > 0)] <- 0 # had offspring <- 0
offspring_binary[is.na(offspring_binary)] <- 1 # had offspring <- 1

plot(response~orig_Fgrm, data = emm_results_zi_df, type = "l", col = "black", 
     xlab = "Fgrm", ylab = "Probability of excess zero", 
     cex.axis = 1.2, 
     cex.lab = 1.75, 
     xlim = range(emm_results_zi_df$orig_Fgrm), 
     ylim = c(0,1))
points(offspring_binary~data_flt$Fgrm, pch = 19, col = alpha("black", alpha = 0.1))
polygon(x=c(emm_results_zi_df$orig_Fgrm, rev(emm_results_zi_df$orig_Fgrm)), 
        y = c(emm_results_zi_df$asymp.UCL, rev(emm_results_zi_df$asymp.LCL)), 
        col = alpha("gray", alpha = 0.5), border = NA)

# I like this plot
# export data to remake in make_figures.R
zi_plot <- list(emm_results_zi_df, offspring_binary, data_flt, emm_results_zi_obs_df)
save(zi_plot, file = "../inbreeding_models/zi_Fgrm_plot_data.Robj")


# get predicted probability of excess zero from observed Fgrm values 

emm_results_zi_obs <- emmeans(mod2, ~ Fgrm, 
                          at = list(Fgrm = data_scaled$Fgrm), 
                          component = "zi", type = "response")
emm_results_zi_obs_df <- as.data.frame(emm_results_zi_obs)
emm_results_zi_obs_df$orig_Fgrm <- emm_results_zi_obs_df$Fgrm * sd(data_flt$Fgrm)

points(response~orig_Fgrm, data = emm_results_zi_obs_df)



# another way
# Plot effect

hist(predict(mod2, type = "zprob"))

newdata = data_scaled[,2:10]
newdata$prob_zi <- predict(mod2, newdata, type = "zprob")
newdata$spacedFgrm <- seq(from = range(data_scaled$Fgrm)[1], to = range(data_scaled$Fgrm)[2], length.out = nrow(newdata))
plot(prob_zi~spacedFgrm, data = newdata, pch = 19, col = "black", 
     xlab = "Fgrm", ylab = "Probability", 
     cex.axis = 1.2, 
     cex.lab = 1.75, 
     xlim = range(data_scaled$Fgrm), 
     ylim = c(0,1))
j <- order(newdata$spacedFgrm)
lines(newdata$spacedFgrm, loess(prob_zi~spacedFgrm, data = newdata)$fitted[j])
# eeek, this is not the visualization that I want! 

emmeans(mod1_noOut, ~ Fgrm:sex, component = "cond")
emmeans(mod1, ~ Fgrm:sex, component = "cond", at = list(Fgrm = 0))
emmeans(mod1, ~ Fgrm:sex, component = "cond", at = list(Fgrm = 0.3))

emmeans(mod1, ~ Fgrm:sex, component = "zi")
emmeans(mod1, ~ Fgrm:sex, component = "zi", at = list(Fgrm = 0))
emmeans(mod1, ~ Fgrm:sex, component = "zi", at = list(Fgrm = 0.1))
emmeans(mod1, ~ Fgrm:sex, component = "zi", at = list(Fgrm = 0.3))
# males have more offspring than females when Fgrm is 0
# males have more offspring than females when Fgrm is 0.1
# males only have slightly more offspring than females when Fgrm is 0.3
contrast(emmeans(mod1, ~ Fgrm:sex, component = "zi", at = list(Fgrm = 0)), method = "pairwise")

emmeans(mod1, ~ Fgrm, component = "zi")

pdf(file="../inbreeding_models/emmeans_mod2.pdf", height = 6, width = 10)
emm_results_cond <- emmeans(mod2, ~ Fgrm, at = list(Fgrm = round(seq(-2.4, 6, length.out = 15), digits = 2)), component = "cond")
plot(emm_results_cond, horizontal = FALSE, 
     xlab = "emmean (predicted offspring on the log scale), conditional model", 
     ylab = "scaled Fgrm", 
     comparisons = FALSE)
emm_results_zi <- emmeans(mod2, ~ Fgrm, at = list(Fgrm = round(seq(-2.4, 6, length.out = 15), digits = 2)), component = "zi")
plot(emm_results_zi, horizontal = FALSE, 
     xlab = "log-odds of observation being an excess zero, zi model", 
     ylab = "scaled Fgrm", 
     comparisons = FALSE)
dev.off()

# model trend figure

e <- allEffects(mod = mod2, xlevels = 200) %>% as.data.frame() 

Fgrm_effects <- e$`Fgrm`
Fgrm_effects$Fgrm <- Fgrm_effects$Fgrm * sd(data_flt$Fgrm)

palette = "Hiroshige" #"OKeeffe2"
colors <- met.brewer(palette, n = 10)

p1 = ggplot(data = Fgrm_effects, aes(x=Fgrm, y=fit)) +
  geom_point(data = data_flt, aes(x=Fgrm, y=offspring, fill = sex), shape = 21, stroke = 0, color = "black", size = 2.5, alpha = 0.5) +
  geom_ribbon(aes(ymax=lower, ymin=upper), alpha=0.25) +
  geom_line() + 
  ylab("Offspring") + xlab("Fgrm") +
  scale_x_continuous(breaks=seq(-0.10, 0.3, by=0.05), limits=c(-0.0946,0.315))+
  scale_y_continuous(breaks=seq(0.0, 20, by=5), limits=c(0.0, 20)) +
  theme_classic() + theme(text = element_text(size=14),axis.text = element_text(colour='black'),
                          #legend.position.inside = c(0.9, 0.9),
                          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  labs(x = "Fgrm", y = "Number of offspring") + 
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) + theme(legend.position="none") +
  scale_fill_manual(values = c("1" = colors[2] ,"2" = colors[4])) 
#guides(fill = guide_legend(title = "sex", position = "inside"), linetype = guide_legend(title = "Sex", position = "inside"))

tiff(file="../inbreeding_models/mod2_Fgrm.tiff", 
     width=7, height=5, units="in", res=600)
p1
dev.off()

plot_list <- list(Fgrm_effects, data_flt)
save(plot_list, file = "../inbreeding_models/Fgrm_plot.Robj")
# years 
year_effects <- e$`yearsContribute2Ped`
year_effects$yearsContribute2Ped <- year_effects$yearsContribute2Ped * sd(data_flt$yearsContribute2Ped)

fig_years = ggplot(data = year_effects, aes(x=yearsContribute2Ped, y=fit)) +
  geom_point(data = data_flt, aes(x=yearsContribute2Ped, y=offspring, fill = sex), shape = 21, stroke = 0, color = "black", size = 2) +
  geom_ribbon(aes(ymax=lower, ymin=upper), alpha=0.5) +
  geom_line() + 
  ylab("Offspring") + xlab("yearsContribute2Ped") +
  scale_x_continuous(breaks=seq(0, 25, by=5), limits=c(0,26))+
  scale_y_continuous(breaks=seq(0.0, 20, by=5), limits=c(0.0, 20)) +
  theme_classic() + theme(text = element_text(size=14),axis.text = element_text(colour='black'),
                          #legend.position.inside = c(0.9, 0.9),
                          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  labs(x = "yearsContribute2Ped", y = "Number of offspring") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) + theme(legend.position="none") +
  scale_fill_manual(values = c("1" = "turquoise3","2" = "gray")) 
#guides(fill = guide_legend(title = "sex", position = "inside"), linetype = guide_legend(title = "Sex", position = "inside"))

fig_years

#-----------------------------------------------------------------------------------------------------------------------------

# Effect sizes for mod1
#------------------------------------------------------------------------------------------------------------------------------
# https://www.biorxiv.org/content/biorxiv/suppl/2017/05/01/132753.DC1/132753-2.pdf
library(plyr)
simple_mod <- glmmTMB(offspring~sex*Fgrm, zi = ~sex*Fgrm, data = data_scaled, family = nbinom2)
newdata0 = newdata = unique(data_scaled[,c("Fgrm","sex")])
temp = predict(simple_mod, newdata, se.fit=TRUE, type="response")
newdata$predFE = temp$fit
newdata$predFE.min = temp$fit-1.98*temp$se.fit
newdata$predFE.max = temp$fit+1.98*temp$se.fit
real=ddply(data_scaled, ~site+Fgrm+sex, summarize, m=mean(offspring))
ggplot(newdata, aes(Fgrm, predFE, colour=Fgrm))+geom_point()+
  geom_errorbar(aes(ymin=predFE.min, ymax=predFE.max))+
  geom_point(data=real, aes(x=Fgrm, y=m) )+
  ylab("Average Offspring")+
  xlab("Sex")

newdata0 = newdata = unique(data_scaled[,c("Fgrm","sex","site","PC2","PC3","PC4","PC5","PC6","yearsContribute2Ped")])
X.cond = model.matrix(lme4::nobars(formula(mod1)[-2]), newdata0)
beta.cond = fixef(mod1)$cond
pred.cond = X.cond %*% beta.cond
ziformula = mod1$modelInfo$allForm$ziformula
X.zi = model.matrix(lme4::nobars(ziformula), newdata0)
beta.zi = fixef(mod1)$zi
pred.zi = X.zi %*% beta.zi

pred.ucount = exp(pred.cond)*(1-plogis(pred.zi))


library(MASS)
set.seed(101)
pred.condpar.psim = mvrnorm(1000,mu=beta.cond,Sigma=vcov(mod1)$cond)
pred.cond.psim = X.cond %*% t(pred.condpar.psim)
pred.zipar.psim = mvrnorm(1000,mu=beta.zi,Sigma=vcov(mod1)$zi)
pred.zi.psim = X.zi %*% t(pred.zipar.psim)
pred.ucount.psim = exp(pred.cond.psim)*(1-plogis(pred.zi.psim))
ci.ucount = t(apply(pred.ucount.psim,1,quantile,c(0.025,0.975)))
ci.ucount = data.frame(ci.ucount)
names(ci.ucount) = c("ucount.low","ucount.high")
pred.ucount = data.frame(newdata0, pred.ucount, ci.ucount)

real.count = ddply(data_scaled, ~Fgrm+sex, summarize, m=median(offspring), mu=mean(offspring))
# ggplot(pred.ucount, aes(x=Fgrm, y=pred.ucount, colour=sex))+ # +geom_point(shape=1, size=2)
#   geom_line(aes(linetype=sex)) + 
#   geom_ribbon(aes(fill=sex, ymin=ucount.low, ymax=ucount.high), alpha = 0.5)+
#   geom_point(data=real.count, aes(x=Fgrm, y=m, colour=sex), shape=0, size=2)+
#   geom_point(data=real.count, aes(x=Fgrm, y=mu, colour=sex), shape=5, size=2)+
#   ylab("Offspring")+
#   xlab("Fgrm")

plot(pred.ucount~Fgrm, data = pred.ucount)
points(offspring~Fgrm, data = data_scaled, col = "red")
#------------------------------------------------------------------------------------------------------------------------------
fixef(mod1) # cond.Fgrm: -0.01007, zi: 0.603398
fixef(mod2) # cond.Fgrm: -0.120858, zi: 0.24137

mod1_sum <- summary(mod1)

mod1_sum$coefficients$cond[2] # scaled effect size of Fgrm for females
mod1_sum$coefficients$cond[2] + mod1_sum$coefficients$cond[11] # scaled effect size of Fgrm for males

mod1_sum$coefficients$zi[2] # scaled effect size of Fgrm for females
mod1_sum$coefficients$zi[2] + mod1_sum$coefficients$zi[11] # scaled effect size of Fgrm for males


# divide effect size by standard deviation of original data 
sd(data_flt$Fgrm)

# cond
mod1_sum$coefficients$cond["Fgrm", "Estimate"] / sd(data_flt$Fgrm) # effect size of Fgrm for females
(mod1_sum$coefficients$cond["Fgrm", "Estimate"] + mod1_sum$coefficients$cond["Fgrm:sex2", "Estimate"]) / sd(data_flt$Fgrm) # effect size of Fgrm for males

# zi
Beta_F_female <- mod1_sum$coefficients$zi["Fgrm", "Estimate"] / sd(data_flt$Fgrm) # effect size of Fgrm for females 
Beta_F_male <- (mod1_sum$coefficients$zi["Fgrm", "Estimate"] + mod1_sum$coefficients$zi["Fgrm:sex2", "Estimate"]) / sd(data_flt$Fgrm) # effect size of Fgrm for males

exp(Beta_F_female)^0.1
exp(Beta_F_male)^0.1

(1- exp(Beta_F_female)^0.1)*100 # 0.1 change: 373.328% less likely to have offspring
(1- exp(Beta_F_male)^0.1)*100 # 0.1 change: 30% more likely to have offspring


mod2_sum <- summary(mod2)
mod2_sum$coefficients$cond["Fgrm", "Estimate"] / sd(data_flt$Fgrm)
mod2_sum$coefficients$zi["Fgrm", "Estimate"] / sd(data_flt$Fgrm)

exp(mod2_sum$coefficients$zi["Fgrm", "Estimate"] / sd(data_flt$Fgrm))
(1- exp(mod2_sum$coefficients$zi["Fgrm", "Estimate"] / sd(data_flt$Fgrm))^0.1)*100 # 0.1 change: 86% less likely to have offspring

# odds ratios are excessively large... 
#------------------------------------------------------------------------------------------------------------------------------

data_noOut <- subset(data_scaled, offspring < 15)
plot(offspring~Fgrm, data = data_noOut)


# make models

mod1_noOut <- glmmTMB(offspring ~ Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                ziformula = ~., 
                data = data_noOut, 
                family = nbinom2)
# check model output
plot(DHARMa::simulateResiduals(mod1_noOut))
diagnose(mod1_noOut) # Unusually large Z-statistics
DHARMa::testDispersion(mod1_noOut)


# look at important effects
drop1(mod1_noOut, test = "Chisq") # p-value tells you term is explaining significant variation in the model
summary(mod1_noOut) # p-values not super informative, but effect sizes are (but they are scaled)
plot(allEffects(mod1_noOut)) 

mod2_noOut <- glmmTMB(offspring ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                ziformula = ~., 
                data = data_noOut, 
                family = nbinom2)

# check model output
plot(DHARMa::simulateResiduals(mod2_noOut))
diagnose(mod2_noOut) # Unusually large Z-statistics
DHARMa::testDispersion(mod2_noOut)

# look at important effects
drop1(mod2_noOut, test = "Chisq") # p-value tells you term is explaining significant variation in the model
summary(mod2_noOut) # p-values not super informative, but effect sizes are (but they are scaled)
plot(allEffects(mod2_noOut)) 

mod3_noOut <- glmmTMB(offspring ~ sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                ziformula = ~., 
                data = data_noOut, 
                family = nbinom2)

# check model output
plot(DHARMa::simulateResiduals(mod3_noOut))
diagnose(mod3_noOut) # Unusually large Z-statistics
DHARMa::testDispersion(mod3_noOut)

# look at important effects
drop1(mod3_noOut, test = "Chisq") # p-value tells you term is explaining significant variation in the model
summary(mod3_noOut) # p-values not super informative, but effect sizes are (but they are scaled)
plot(allEffects(mod3_noOut)) 

# use AICc to evaluate the best supported model: 
mod_list_noOut <- list(mod1_noOut, mod2_noOut, mod3_noOut)
AICcmodavg::aictab(mod_list_noOut) # mod1, mod2, mod3

# are odds ratios more reasonable? 
# cond
mod1_noOut_sum <- summary(mod1_noOut)
mod1_noOut_sum$coefficients$cond["Fgrm", "Estimate"] / sd(data_flt$Fgrm) # effect size of Fgrm for females
(mod1_noOut_sum$coefficients$cond["Fgrm", "Estimate"] + mod1_noOut_sum$coefficients$cond["Fgrm:sex2", "Estimate"]) / sd(data_flt$Fgrm) # effect size of Fgrm for males

# zi
Beta_F_female <- mod1_noOut_sum$coefficients$zi["Fgrm", "Estimate"] / sd(data_noOut$Fgrm) # effect size of Fgrm for females 
Beta_F_male <- (mod1_noOut_sum$coefficients$zi["Fgrm", "Estimate"] + mod1_noOut_sum$coefficients$zi["Fgrm:sex2", "Estimate"]) / sd(data_noOut$Fgrm) # effect size of Fgrm for males

exp(Beta_F_female)^0.1
exp(Beta_F_male)^0.1

(1- exp(Beta_F_female)^0.1)*100 # 0.1 change: 6.426641% less likely to have offspring
(1- exp(Beta_F_male)^0.1)*100 # 0.1 change: 0.2724786% more likely to have offspring

# mod 2 
mod2_noOut_sum <- summary(mod2_noOut)
mod2_noOut_sum$coefficients$cond["Fgrm", "Estimate"] / sd(data_scaled$Fgrm) 

# zi
mod2_noOut_sum$coefficients$zi["Fgrm", "Estimate"] / sd(data_scaled$Fgrm) 
exp(mod2_noOut_sum$coefficients$zi["Fgrm", "Estimate"] / sd(data_scaled$Fgrm))^0.1

# that looks reasonable

# outliers in data appear to be responsible for odds ratios...  results are consistent without outliers, with is reassuring... 

# potentially good info here: https://stats.stackexchange.com/questions/38493/how-to-deal-with-quasi-complete-separation-in-a-logistic-glmm


# mod2 
mod2_sum <- summary(mod2)
mod2_sum$coefficients$zi["Fgrm", "Estimate"] / sd(data_flt$Fgrm)
exp(mod2_sum$coefficients$zi["Fgrm", "Estimate"]/ sd(data_flt$Fgrm)) / 100
# for every 1.0 increase in Fgrm, the probability of being an excess 0 increases by 5


# confidence estimates: 
# mod1
mod1_conf <- as.data.frame(confint(mod1))
mod1_conf_85 <- as.data.frame(confint(mod1, level = 0.85))
mod1_conf <- cbind.data.frame(rownames(mod1_conf), mod1_conf, mod1_conf_85$`7.5 %`, mod1_conf_85$`92.5 %`)
colnames(mod1_conf) <- c("factor", "LCI_95", "UCI_95", "estimate", "LCI_85", "UCI_85")
mod1_conf <- mod1_conf[-c(1, 12),] # remove estimates for intercepts 

# mod1_noOut
mod1_noOut_conf <- as.data.frame(confint(mod1_noOut))
mod1_noOut_conf_85 <- as.data.frame(confint(mod1_noOut, level = 0.85))
mod1_noOut_conf <- cbind.data.frame(rownames(mod1_noOut_conf), mod1_noOut_conf, mod1_noOut_conf_85$`7.5 %`, mod1_noOut_conf_85$`92.5 %`)
colnames(mod1_noOut_conf) <- c("factor", "LCI_95", "UCI_95", "estimate", "LCI_85", "UCI_85")
mod1_noOut_conf <- mod1_noOut_conf[-c(1, 12),] # remove estimates for intercepts 


# Betas figure

mod1_betas <-  ggplot(data =  mod1_conf) + geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 10.5, color = "black", linetype = 2) +
  geom_linerange(aes(x = factor, 
                     y = estimate, ymin = LCI_85, ymax = UCI_85),size=1.25) + 
  geom_pointrange(aes(x = factor, 
                      y = estimate, ymin = LCI_95, ymax = UCI_95),size=.5) + 
  coord_flip() +
  labs(x = NULL, y = NULL)+ theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(color = "black", size=14)) +
  theme(axis.text.y = element_text(color = "black", size=14)) +
  scale_y_continuous(breaks=seq(-4, 4, by=1), limits=c(-4, 4))+
  scale_x_discrete(labels=factor)

mod1_betas

mod1_noOut_betas <-  ggplot(data =  mod1_noOut_conf) + geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 10.5, color = "black", linetype = 2) +
  geom_linerange(aes(x = factor, 
                     y = estimate, ymin = LCI_85, ymax = UCI_85),size=1.25) + 
  geom_pointrange(aes(x = factor, 
                      y = estimate, ymin = LCI_95, ymax = UCI_95),size=.5) + 
  coord_flip() +
  labs(x = NULL, y = NULL)+ theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(color = "black", size=14)) +
  theme(axis.text.y = element_text(color = "black", size=14)) +
  scale_y_continuous(breaks=seq(-4, 4, by=1), limits=c(-4, 4))+
  scale_x_discrete(labels=factor)

mod1_noOut_betas

tiff(file="../figures/SupplementalFig11.tiff", 
     width=7, height=5, units="in", res=300)

p2
dev.off()

# effect size figure 

emmeans(mod1_noOut, ~ Fgrm:sex, component = "cond")
emmeans(mod1, ~ Fgrm:sex, component = "cond", at = list(Fgrm = 0))
emmeans(mod1, ~ Fgrm:sex, component = "cond", at = list(Fgrm = 0.3))

emmeans(mod1, ~ Fgrm:sex, component = "zi")
emmeans(mod1, ~ Fgrm:sex, component = "zi", at = list(Fgrm = 0))
emmeans(mod1, ~ Fgrm:sex, component = "zi", at = list(Fgrm = 0.1))
emmeans(mod1, ~ Fgrm:sex, component = "zi", at = list(Fgrm = 0.3))
# males have more offspring than females when Fgrm is 0
# males have more offspring than females when Fgrm is 0.1
# males only have slightly more offspring than females when Fgrm is 0.3
contrast(emmeans(mod1, ~ Fgrm:sex, component = "zi", at = list(Fgrm = 0)), method = "pairwise")

emmeans(mod1, ~ Fgrm, component = "zi")

pdf(file="../inbreeding_models/emmeans_mod1_nodOut.pdf", height = 6, width = 10)
emm_results_cond <- emmeans(mod1_noOut, ~ Fgrm*sex, at = list(Fgrm = round(seq(-2.4, 6, length.out = 10), digits = 2)), component = "cond")
plot(emm_results_cond, horizontal = FALSE, 
     xlab = "emmean (predicted offspring on the log scale), conditional model", 
     ylab = "scaled Fgrm", 
     comparisons = FALSE)
emm_results_zi <- emmeans(mod1_noOut, ~ Fgrm*sex, at = list(Fgrm = round(seq(-2.4, 6, length.out = 10), digits = 2)), component = "zi")
plot(emm_results_zi, horizontal = FALSE, 
     xlab = "log-odds of observation being an excess zero, zi model", 
     ylab = "scaled Fgrm", 
     comparisons = FALSE)

emm_results_cond_onlyF <- emmeans(mod1_noOut, ~ Fgrm, at = list(Fgrm = round(seq(-2.4, 6, length.out = 10), digits = 2)), component = "cond")
plot(emm_results_cond_onlyF, horizontal = FALSE, 
     xlab = "emmean (predicted offspring on the log scale), conditional model", 
     ylab = "scaled Fgrm", 
     comparisons = FALSE)
emm_results_zi_onlyF <- emmeans(mod1_noOut, ~ Fgrm, at = list(Fgrm = round(seq(-2.4, 6, length.out = 10), digits = 2)), component = "zi")
plot(emm_results_zi_onlyF, horizontal = FALSE, 
     xlab = "log-odds of observation being an excess zero, zi model", 
     ylab = "scaled Fgrm", 
     comparisons = FALSE)
dev.off()
emm_results[1][1]

summary(emmeans(mod1_noOut, "sex"), infer = TRUE, type = "response")


plot(emmeans(mod1_NB, ~ Fgrm:sex))

mod1_means = emmeans(mod1, ~ Fgrm:sex)
contrast(mod1_means, method = "pairwise")

plot(allEffects(mod2_NB))
e <- allEffects(mod = mod1_noOut, xlevels = 200) %>% as.data.frame() 
effect("Fgrm", mod1_noOut, xlevels = 200)

Fgrm_effects <- e$`Fgrm:sex`
Fgrm_effects$Fgrm <- Fgrm_effects$Fgrm * sd(data_flt$Fgrm)

p1 = ggplot(data = Fgrm_effects, aes(x=Fgrm, y=fit)) +
  geom_point(data = data_flt, aes(x=Fgrm, y=offspring, fill = sex), shape = 21, stroke = 0, color = "black", size = 2) + 
  geom_ribbon(aes(fill=sex, ymax=lower, ymin=upper), alpha=0.5) +
  geom_line(aes(linetype=sex)) + 
  ylab("Offspring") + xlab("Fgrm") +
  scale_x_continuous(breaks=seq(-0.10, 0.3, by=0.05), limits=c(-0.09459466,0.25741486))+
  scale_y_continuous(breaks=seq(0.0, 20, by=5), limits=c(0.0, 20)) +
  theme_classic() + theme(text = element_text(size=14),axis.text = element_text(colour='black'),
                          #legend.position.inside = c(0.9, 0.9),
                          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  labs(x = "Fgrm", y = "Number of offspring") + 
  theme(axis.text.x=element_text(angle=45,hjust=1)) + theme(legend.position="none") +
  scale_fill_manual(values = c("1" = "turquoise3","2" = "gray")) 
  #guides(fill = guide_legend(title = "sex", position = "inside"), linetype = guide_legend(title = "Sex", position = "inside"))

p1

#------------------------------------------------------------------------------------------------------------------------------
# I would like to try to visualize these results somehow, but I'm not sure how to proceed with a 
# zero-inflation model. A sentence like "With every increase of 0.1 in Fgrm, the probability of having 
# offspring goes down by X. 

# I can make a plot using ggpredit, but I think this is visualizing the slope for
# the count model. 
plot(ggpredict(mod1, terms = c("Fgrm", "site"), type = "fixed"))

# Summary/Questions:
# (1) I'm curious if you have any thoughts or opinions on our modeling approach, and if you would recommend 
#     changing anything. 
# (2) I don't know how to interpret the "Unusually large coefficients" and "Unusually large Z-statistics" warnings
#     from the diagnose() function. How concerned should I be about these, and is there anything I can do about it? 
# (1) It's unclear to me how to proceed with model selection given the similar deltaAICc values. I've run 
#     some non-parametric permutation tests, and the effect of Fgrm on whether or not an individual 
#     has offspring seems real, but it doesn't seem like this is the conclusion from the parametric tests. 
# (2) how can I get a summary statement about the effect of Fgrm on offspring in the zero-inflation model?
# (3) although I don't present these data here (still in the works!), I eventually would like to incorporate 
#     individual longevity as either a fixed effect or an offset in the model. We expect that longevity will 
#     impact the number of offspring an individual has (older individuals will reproduce more), and whether 
#     an individual will have offspring (some individuals may have died before reaching sexual maturity). I
#     am not sure how glmmTMB interprets offset parameters in the zero inflation model, so any thoughts
#     about the best way to approach this would be helpful! 


#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------

# MISC
# convo with rachel 
# time to reproductive size ~ Fgrm 
# are direction of effects consistent between models? 
# plot raw data... does that make sense? 
# interaction term: 
# 
plot(mod1)

# (1) use drop1 for p-values
drop1(mod1, test = "Chisq") # not sure what this output indicates, p-value tells you term is explaining significant variation in the model
drop1(mod2, test = "Chisq") # fitting model with everything, and then everything except for one factor and doing likelihood ratios between those two

# (2) summary to look at direction and magnitude of effects NOT p-values, p-values are in relation to intercept or something? 

# plot()
plot(allEffects(mod1)) 


data_noZero <- subset(data_flt, offspring > 0)
plot(offspring~Fgrm, data_noZero)

mod1 <- glmmTMB(offspring ~ Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                data = data_noZero, 
                family = nbinom2)
plot(allEffects(mod1)) 

plot(offspring~Fgrm, data = subset(data_noZero, sex == 1), main = "females")
plot(offspring~Fgrm, data = subset(data_noZero, sex == 2), main = "males")


# remove outliers and see what happens
data_noZero_noOut <- subset(data_noZero, sex == 2 && offspring < 12)
plot(offspring~Fgrm, data = data_noZero_noOut, main = "males")

mod1 <- glmmTMB(offspring ~ Fgrm + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                data = data_noZero_noOut, 
                family = nbinom2)

summary(mod1)

# remove outliers but keep zeros, just males
data_noOut <- subset(data_flt, sex == 2 & offspring < 12)

mod1 <- glmmTMB(offspring ~ Fgrm + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                ziformula = ~., 
                data = data_noOut, 
                family = nbinom2)
summary(mod1)

subset(data_flt, offspring > 12)

# remove outliers but keep zeros, males and females
data_noOut <- subset(data_flt, offspring < 12)
plot(offspring~Fgrm, data = data_noOut)
mod1 <- glmmTMB(offspring ~ Fgrm * sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + yearsContribute2Ped, 
                ziformula = ~., 
                data = data_noOut, 
                family = nbinom2)
summary(mod1)
drop1(mod1, test = "Chisq")

plot(allEffects(mod1))
# effect size, directions, rank order, raw data, relationships vs model outputs
# pull out predicted values, values that the model would predict for dataset, and plot, compare model1 and model 2, 
# are they correlated? plot zero inflated vs non-zero inflated 
predict(mod1, type = "response", re.form = NA)



# plot(simulateresiduals)

plot(offspring~Fgrm, data = data_flt)
abline(lm(offspring~Fgrm, data = data_flt))

plot(offspring~Fgrm, data = subset(data_flt, sex == 2))
abline(lm(offspring~Fgrm, data = subset(data_flt, sex == 2)))

plot(offspring~Fgrm, data = subset(data_flt, sex == 1))
abline(lm(offspring~Fgrm, data = subset(data_flt, sex == 1)))

# look into how to interpret effect sizes of interactions! need to add or subtract interaction term from one of the main predictors! 


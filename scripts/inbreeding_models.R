### File set up  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
library(DHARMa)
library(glmmTMB)
library(GGally)
library(stringr)

# date <- "05262023"
# date <- "04132024"
date <- "05072024"

load(file = "../vcf_filtering/pop_gen_stats_02082024.Robj", verbose = T)
pop_gen_stats$id <- str_remove(pop_gen_stats$id, "_P[0-9]+")

load(paste0("../pedigree_reconstruction/ELF_pedigree_results_", date, ".Robj"), verbose = T)
load(paste0("../pedigree_reconstruction/PCC_pedigree_results_", date, ".Robj"), verbose = T)

# metadata from pedigree reconstruction
ELF_meta <- ELF_results[[2]]
PCC_meta <- PCC_results[[2]]
PCC_meta$ID <- str_remove(PCC_meta$ID, "_P[0-9]+")
ELF_meta$ID <- str_remove(ELF_meta$ID, "_P[0-9]+")
#ELF_meta <- ELF_meta[-which(ELF_meta$ID =="ELF_504"),]
#ELF_meta <- ELF_meta[-which(ELF_meta$ID =="ELF_5440"),]

# pedigrees
PCC_pedigree <- PCC_results[[4]]$Pedigree
load(paste0("../pedigree_exploration/EMR_ped_flt_", date, ".Robj"), verbose = TRUE)
ELF_pedigree <- ELF_pedigree_flt

# fix tech rep ids 
PCC_pedigree$id <- str_remove(PCC_pedigree$id, "_P[0-9]+")
PCC_pedigree$dam <- str_remove(PCC_pedigree$dam, "_P[0-9]+")
PCC_pedigree$sire <- str_remove(PCC_pedigree$sire, "_P[0-9]+")

ELF_pedigree$id <- str_remove(ELF_pedigree$id, "_P[0-9]+")
ELF_pedigree$dam <- str_remove(ELF_pedigree$dam, "_P[0-9]+")
ELF_pedigree$sire <- str_remove(ELF_pedigree$sire, "_P[0-9]+")

# set up vectors of individuals
PCC_inds <- pop_gen_stats$id[which(pop_gen_stats$site == "PCC")]
ELF_inds <- pop_gen_stats$id[which(pop_gen_stats$site == "ELF")]

load(paste0("../pedigree_exploration/ELF_expanded_metaData.Robj"), verbose = T) # ELF_exp_metaData
load(paste0("../pedigree_exploration/PCC_expanded_metaData.Robj"), verbose = T) # PCC_exp_metaData

PCC_dupe_inds <- read.csv("../pedigree_reconstruction/PCC_LifeHistData_dupes.csv", header=FALSE)


# ------------------------------------------------------------------------------------------------------------

## Define custom functions ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
source("./pedigree_funcs.R") 

# ------------------------------------------------------------------------------------------------------------

### Get number of offspring  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
PCC_off <- offspring_per_ind(PCC_pedigree, plot = F)
ELF_off <- offspring_per_ind(ELF_pedigree, plot = F)
# ------------------------------------------------------------------------------------------------------------

### Compile data  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# combine pop gen data
data <- cbind.data.frame(pop_gen_stats$id, pop_gen_stats$site, pop_gen_stats$Fgrm)
colnames(data) <- c("id", "site", "Fgrm")

# format offspring data and combine
moms <- c(ELF_off[[1]], PCC_off[[1]])
moms <- moms[which(grepl("_", names(moms)))]
dads <- c(ELF_off[[2]], PCC_off[[2]])
dads <- dads[which(grepl("_", names(dads)))]
offspring_data <- cbind.data.frame(c(names(moms), names(dads)), c(moms, dads))
colnames(offspring_data) <- c("id", "offspring")

data <- merge(x = data, y = offspring_data, by = "id", all.x = TRUE)

# add sex and known age information
meta_mega <- rbind.data.frame(ELF_meta, PCC_meta)
data <- merge(x = data, y = meta_mega, by.x = "id", by.y = "ID", all.x = TRUE)

# covert -9 to NA in sex and BirthYear columns and make NA for offspring 0
data[which(is.na(data$offspring)),"offspring"] <- 0

head(data)
plot(data$offspring~data$Fgrm)

# remove individual with high Fgrm (likely migrant?)

data <- subset(data, Fgrm < 0.8)
# ------------------------------------------------------------------------------------------------------------
 
### Initial strategy: Age estimates  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
ELF.by.probs <- CalcBYprobs(Pedigree = ELF_pedigree, LifeHistData = ELF_meta)
PCC.by.probs <- CalcBYprobs(Pedigree = PCC_pedigree, LifeHistData = PCC_meta)

ELF.by.est <- vector(length = nrow(ELF.by.probs))
for(i in 1:nrow(ELF.by.probs)){
  ind <- rownames(ELF.by.probs)[i]
  if(sum(ELF.by.probs[i,]) > 0){
    if(max(ELF.by.probs[i,])  > 0.8){ # if there is a probable birth year greater than 80%, record as birth year
      ELF.by.est[i] <- as.numeric(colnames(ELF.by.probs)[which(ELF.by.probs[i,]==max(ELF.by.probs[i,]))])
    }else if(sum(ELF.by.probs[i,] > 0) <= 7){ # if there are less than or equal to 5 birth years with probabilities... 
      ELF.by.est[i] <- sum(as.numeric(names(ELF.by.probs[i,]))*ELF.by.probs[i,])/sum(ELF.by.probs[i,]) # otherwise calculate a weighted mean to use as birth year
    }else{ # if the probability distribution is too flat, throw it out! 
      ELF.by.est[i] <- NA
    }
  }else if(sum(ELF.by.probs[i,]) == 0){
    ELF.by.est[i] <- NA
  }
}
ELF.by.est <- cbind.data.frame(ELF.by.est, rownames(ELF.by.probs))
colnames(ELF.by.est) <- c("by.est", "id")

sum(!is.na(ELF.by.est$by.est)) # 198 to 171 with more stringent filters, 526 run: 326???

PCC.by.est <- vector(length = nrow(PCC.by.probs))
for(i in 1:nrow(PCC.by.probs)){
  ind <- rownames(PCC.by.probs)[i]
  if(sum(PCC.by.probs[i,]) > 0){
    if(max(PCC.by.probs[i,])  > 0.8){ # if there is a probable birth year greater than 80%, record as birth year
      PCC.by.est[i] <- as.numeric(colnames(PCC.by.probs)[which(PCC.by.probs[i,]==max(PCC.by.probs[i,]))])
    }else if(sum(PCC.by.probs[i,] > 0) <= 7){ # if there are less than or equal to 5 birth years with probabilities... 
      PCC.by.est[i] <- sum(as.numeric(names(PCC.by.probs[i,]))*PCC.by.probs[i,])/sum(PCC.by.probs[i,]) # otherwise calculate a weighted mean to use as birth year
    }else{ # if the probability distribution is too flat, throw it out! 
      PCC.by.est[i] <- NA
    }
  }else if(sum(PCC.by.probs[i,]) == 0){
    PCC.by.est[i] <- NA
  }
}
PCC.by.est <- cbind.data.frame(PCC.by.est, rownames(PCC.by.probs))
colnames(PCC.by.est) <- c("by.est", "id")

sum(!is.na(PCC.by.est$by.est)) # 173 to 35 with more stringent filters, 526: 72


by.est <- rbind.data.frame(PCC.by.est, PCC.by.est)

# double checking all inds with estimated by don't have a confirmed by 
# sum(is.na(data[match(by.est$id, data$id),"BirthYear"])) == length(by.est$id) # confirmed! 

by.est <- by.est[which(grepl("_", by.est$id)),] # filter out dummy individuals

# add estimated birth years to data object
data$est.BY <- data$BirthYear
data[match(by.est$id, data$id),"est.BY"] <- by.est[na.omit(match(data$id, by.est$id)),"by.est"] 
sum(!is.na(data$BirthYear)) # 558 with known BY
sum(!is.na(data$est.BY)) # 628 with estimated birth years



# Age estimates
# for individuals with no age estimates but mult recaps, use minimum age based on last capture and age class at first capture
# for individuals at 10 yrs, inform using recap data, use year after last capture?
# compare hist of ages --> is it more normal/left skewed?
# or for everyone use year after last recap for all individuals, and estimated age a first capture


MakeAgePrior(Pedigree = ELF_pedigree, LifeHistData = ELF_meta) # MaxAgeParent = 8,9
MakeAgePrior(Pedigree = PCC_pedigree, LifeHistData = PCC_meta) # MaxAgeParent = 13,13

# going with a max age of 10? 
data_flt <- data[which(!is.na(data$est.BY)),] # not sure this is correct here


# last sampling year is 2019
data_flt$Age <- NA
data_flt[,"Age"] <- 2019 - as.numeric(data_flt[,"BirthYear"])
hist(as.numeric(data_flt[,"BirthYear"]))
hist(data_flt$Age)

# the above distribution now actually looks okay... maybe just use that? 

# if an individual is over ten, we are using their last year of capture or the last recorded birthyear of their offspring as their "last year"
over_10 <- data_flt[which(data_flt$Age > 10),] 

last_evidence <- vector(length = nrow(over_10))
for(i in 1:nrow(over_10)){
  if(over_10$site[i] == "ELF"){
    capture_data <- ELF_exp_metaData[[over_10$id[i]]] # capture 1 time 2011
    last_cap <- as.numeric(max(capture_data$year, na.rm = T))
    
    # what about offspring birth years?
    offspring <- tryCatch({GetDescendants(over_10$id[i], ELF_pedigree)$offspring}, error = function(e){}) 
  }else if(over_10$site[i] == "PCC"){
    capture_data <- PCC_exp_metaData[[over_10$id[i]]] # capture 1 time 2011
    last_cap <- as.numeric(max(capture_data$year, na.rm = T))
    
    # what about offspring birth years?
    offspring <- tryCatch({GetDescendants(over_10$id[i], PCC_pedigree)$offspring}, error = function(e){}) 
  }
  if(!is.null(offspring)){
    off.by <- vector(length=length(offspring))
    for(o in 1:length(offspring)){
      if(length(data[which(data$id == offspring[o]),"BirthYear"]) == 1){
        off.by[o] <- data[which(data$id == offspring[o]),"BirthYear"] 
      }else{off.by[o] <- NA}
    }
    last_off <- max(as.numeric(off.by), na.rm = T)
  }else{last_off <- NA}
  
  last_evidence[i] <- max(last_off, last_cap, na.rm=T)
}
over_10$last.year <- last_evidence

(1+over_10$last.year) - as.numeric(over_10$BirthYear)

data_flt[which(data_flt$Age > 10),"Age"] <- (1+over_10$last.year) - as.numeric(over_10$BirthYear)

hist(2019 - as.numeric(data_flt[,"BirthYear"]), main = NULL, xlab = "age")
hist(data_flt$Age, add = TRUE, col = "lavender")
legend("topright", legend = c("2019-BY", "inds over 10 re-assessed with last inferred year"), col = c("gray", "lavender"), pch = 15)

# make age and birth year numeric
class(data_flt$Age)
data_flt$BirthYear <- as.numeric(data_flt$BirthYear)

# filter out individulas without sex data
# remove those without known sex information 
data_flt <- data_flt[-which(is.na(data_flt$sex)),]

dim(data_flt)

#------------------------------------------------------------------------------
# calculate all death years via last observation? 

# if an individual is over ten, we are using their last year of capture or the last recorded birthyear of their offspring as their "last year"
over_10 <- data_flt

last_evidence <- vector(length = nrow(over_10))
for(i in 1:nrow(over_10)){
  if(over_10$site[i] == "ELF"){
    capture_data <- ELF_exp_metaData[[over_10$id[i]]] # capture 1 time 2011
    last_cap <- as.numeric(max(capture_data$year, na.rm = T))
    
    # what about offspring birth years?
    offspring <- tryCatch({GetDescendants(over_10$id[i], ELF_pedigree)$offspring}, error = function(e){}) 
  }else if(over_10$site[i] == "PCC"){
    capture_data <- PCC_exp_metaData[[over_10$id[i]]] # capture 1 time 2011
    last_cap <- as.numeric(max(capture_data$year, na.rm = T))
    
    # what about offspring birth years?
    offspring <- tryCatch({GetDescendants(over_10$id[i], PCC_pedigree)$offspring}, error = function(e){}) 
  }
  if(!is.null(offspring)){
    off.by <- vector(length=length(offspring))
    for(o in 1:length(offspring)){
      if(length(data[which(data$id == offspring[o]),"BirthYear"]) == 1){
        off.by[o] <- data[which(data$id == offspring[o]),"BirthYear"] 
      }else{off.by[o] <- NA}
    }
    last_off <- max(as.numeric(off.by), na.rm = T)
  }else{last_off <- NA}
  
  last_evidence[i] <- max(last_off, last_cap, na.rm=T)
}
over_10$last.year <- last_evidence

over_10$age_recap <- over_10$last.year - as.numeric(over_10$BirthYear)

hist(2019 - as.numeric(data_flt[,"BirthYear"]), main = NULL, xlab = "age")
hist(over_10$age_recap, add = TRUE, col = "lavender")
legend("topright", legend = c("2019-BY", "inds over 10 re-assessed with last inferred year"), col = c("gray", "lavender"), pch = 15)


# ------------------------------------------------------------------------------------------------------------

### Age estimates from SVL  --------------------------------------------------------------------

# get SVL data 
# make data frame of asymptotic size values
asymptotic_L <- data.frame(matrix(NA, nrow = 4, ncol = 5))
colnames(asymptotic_L) <- c("site", "sex", "lower_CI", "upper_CI", "L_est")
asymptotic_L[1,] <- c("cass", "female", 57.5, 59.8, 58.7) # 58.7 cm (95% CI = 57.5 – 59.8) 
asymptotic_L[2,] <- c("cass", "male", 60.3, 63.3, 61.8) # 61.8 cm (95% CI = 60.3 – 63.3)
asymptotic_L[3,] <- c("barry", "female", 58.9, 61.6, 60.3) # 60.3 cm (95% CI = 58.9 – 61.6)
asymptotic_L[4,] <- c("barry", "male", 61.8, 65.1, 63.4) # 63.4 cm (95% CI = 61.8 – 65.1) 

asymptotic_L$lower_CI <- as.numeric(asymptotic_L$lower_CI)
asymptotic_L$upper_CI <- as.numeric(asymptotic_L$upper_CI)
asymptotic_L$L_est <- as.numeric(asymptotic_L$L_est)

# generate SVL data frames

# PCC
PCC_pedigree <- PCC_pedigree[grepl("PCC", PCC_pedigree$id),]
PCC_SVL_data <- data.frame(matrix(NA, nrow = length(PCC_pedigree$id), ncol = 4))
colnames(PCC_SVL_data) <- c("ID", "sex", "year", "SVL")
# get smallest size! 
for(i in 1:length(PCC_pedigree$id)){ # for each individual in the pedigree
  
  # get id of focal individual 
  focal_ind_search <- PCC_pedigree$id[i]
  #focal_ind_search <- str_remove(PCC_pedigree$id[i], "_P[0-9]")
  PCC_SVL_data[i, "ID"] <- focal_ind_search
  if(grepl("_", focal_ind_search)){
    # identify minimum SVL 
    focal_SVL <- min(as.numeric(PCC_exp_metaData[[focal_ind_search]]$SVL))
    if(!is.na(focal_SVL)){
      # find sex of focal individual
      PCC_SVL_data$sex[i] <- PCC_meta[which(PCC_meta$ID == PCC_pedigree$id[i]),"Sex"]
      PCC_SVL_data$year[i] <- PCC_exp_metaData[[focal_ind_search]]$year[which(as.numeric(PCC_exp_metaData[[focal_ind_search]]$SVL) == focal_SVL)] 
      PCC_SVL_data$SVL[i] <- focal_SVL
    }
  }
}

PCC_SVL_data$SVL <- as.numeric(PCC_SVL_data$SVL)
PCC_SVL_data$year <- as.numeric(PCC_SVL_data$year)

#ELF
ELF_SVL_data <- data.frame(matrix(NA, nrow = length(ELF_pedigree_flt$id), ncol = 4))
colnames(ELF_SVL_data) <- c("ID", "sex", "year", "SVL")

# get smallest size! 
for(i in 1:length(ELF_pedigree_flt$id)){ # for each individual in the pedigree
  
  # get id of focal individual 
  focal_ind <- ELF_pedigree_flt$id[i]
  focal_ind_search <- str_remove(ELF_pedigree_flt$id[i], "_P[0-9]")
  ELF_SVL_data[i, "ID"] <- focal_ind_search
  if(grepl("_", focal_ind)){
    # identify minimum SVL 
    focal_SVL <- min(as.numeric(ELF_exp_metaData[[focal_ind_search]]$SVL), na.rm = TRUE)
    if (min(as.numeric(ELF_exp_metaData[[focal_ind_search]]$SVL), na.rm = TRUE) == Inf){
      focal_SVL <- NA
    }
    if(!is.na(focal_SVL)){
      # find sex of focal individual
      ELF_SVL_data$sex[i] <- ELF_meta[which(ELF_meta$ID == ELF_pedigree_flt$id[i]),"Sex"]
      ELF_SVL_data$year[i] <- ELF_exp_metaData[[focal_ind_search]]$year[which(as.numeric(ELF_exp_metaData[[focal_ind_search]]$SVL) == focal_SVL)] 
      ELF_SVL_data$SVL[i] <- focal_SVL
    }
  }
}

ELF_SVL_data$SVL <- as.numeric(ELF_SVL_data$SVL)
ELF_SVL_data$year <- as.numeric(ELF_SVL_data$year)

# infer birth years
head(data)
data$SVL_by <- NA
for(i in 1:nrow(data)){
  ind <- data[i,"id"]
  if(data[i,"site"] == "ELF"){
    if(!length(subset(ELF_SVL_data, ID == ind)$SVL) == 0){ # does SVL data exist for individual? 
      if(!is.na(subset(ELF_SVL_data, ID == ind)$SVL)){
        if(!subset(ELF_SVL_data, ID == ind)$sex == -9){ # is sex -9? 
          data$SVL_by[i] <- infer_birth_year(ind, "cass", ELF_SVL_data)[1]
        }
      }
    }else(data$SVL_by[i] <- NA)
  }else if(data[i,"site"] == "PCC"){
    if(!length(subset(PCC_SVL_data, ID == ind)$SVL) == 0){
      if(!is.na(subset(PCC_SVL_data, ID == ind)$SVL)){
        if(!subset(PCC_SVL_data, ID == ind)$sex == -9){
          data$SVL_by[i] <- infer_birth_year(ind, "barry", PCC_SVL_data)[1]
        }
      }
    }else(data$SVL_by[i] <- NA)
  }
}

plot(x = data$BirthYear, y = data$SVL_by, col = as.factor(data$site))
abline(a = 0, b = 1, lty = 2 )

data$SVL_by

plot(x = data$BirthYear, y = data$Fgrm)
plot(x = data$SVL_by, y = data$Fgrm)


### Build models  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

data_flt <- subset(data, !is.na(SVL_by))

# center SVL_by? 
data_flt$SVL_by_centered <- data_flt$SVL_by - mean(data_flt$SVL_by)

# save 
save(data_flt, file = "../inbreeding_models/data_flt_04132024.RObj")
# data_flt <- over_10
# data_flt$Age <- data_flt$age_recap
# check data
plot(data_flt$offspring~data_flt$Fgrm, col = alpha("black", alpha = 0.3), pch = 19, xlab = "Fgrm", ylab = "number of offspring")
plot(data_flt$offspring~data_flt$SVL_by, col = alpha("black", alpha = 0.3), pch = 19, xlab = "birth year", ylab = "number of offspring")

data_flt_old <- subset(data_flt, SVL_by < 2011)
plot(data_flt_old$offspring~data_flt_old$Fgrm, col = alpha("black", alpha = 0.3), pch = 19, xlab = "Fgrm", ylab = "number of offspring")

dim(data_flt_old)
dim(data_flt)


ggpairs(data_flt[,-1])

# make models

# example of zero inflation model from Isabela's code
# ts1 <- glmmTMB(total_seed_number ~ height_at_flowering + nodule_number * ped_inbreeding + I(ped_inbreeding^2) + rhizobia + prop.E +
#                  (1 | cross_ID) + (1|planting_date_corrected_daymonthyear) + (1|inoculation_date_corrected_daymonthyear),
#                ziformula = ~ nodule_number + prop.E,
#                data = d, family=nbinom1)
# summary(ts1) # interaction is sig


mod_lr <- glm(formula = offspring ~ Fgrm + sex + site + SVL_by_centered, data = data_flt_off, family = 'gaussian')

# build poisson regression 

mod1 <- glm(formula = offspring ~ Fgrm + sex + site + SVL_by_centered, data = data_flt, family = 'poisson')
testZeroInflation(mod1)

# mod2: negative binomial,
mod2 <- glmmTMB(offspring ~ Fgrm + sex + site + (1|SVL_by_centered),
                ziformula = ~ Fgrm + sex + site + (1|SVL_by_centered),
                data = data_flt, 
                family = nbinom1)

# mod3: Conway-Maxwell Poisson with Fgrm
mod3 <- glmmTMB(offspring ~ Fgrm + Sex + site,
                ziformula = ~1,
                data = data_flt_old, 
                family = compois)
summary(mod3)
diagnose(mod3)
testDispersion(mod3)
AIC(mod3) # 1316.707

# mod3: Conway-Maxwell Poisson  Fgrm
mod3 <- glmmTMB(offspring ~ Fgrm + Sex + site + Sex*Fgrm,
                ziformula = ~ Fgrm + site + Sex + Sex*Fgrm,
                data = data_flt_old, 
                family = compois)
summary(mod3)
diagnose(mod3)
testDispersion(mod3)
AIC(mod3) # 1311.585

mod4 <- glmmTMB(offspring ~ Fgrm + sex + site,
                ziformula = ~ site + sex,
                data = data_flt_old, 
                family = compois)
summary(mod4)
AIC(mod4) # 1313.806

# poisson with birth year as offset
mod4 <- glmmTMB(offspring ~ Fgrm + sex + site,
                ziformula = ~ Fgrm + site, 
                offset = SVL_by_centered,
                data = data_flt, 
                family = poisson)
diagnose(mod4)
testDispersion(mod4)
plot(simulateResiduals(mod4))

# negative binomial with age as offset
mod5 <- glmmTMB(offspring ~ Fgrm + sex + site,
                ziformula = ~ Fgrm + site, 
                offset = SVL_by_centered,
                data = data_flt, 
                family = nbinom1)
diagnose(mod5)
testDispersion(mod5)
plot(simulateResiduals(mod5))

# negative binomial with age as offset, changing zero inflated formula
mod6 <- glmmTMB(offspring ~ Fgrm + sex + site,
                ziformula = ~site, # site in zero inflated b/c different sampling efforts, anything biologically different between sites
                offset = SVL_by_centered,
                data = data_flt, 
                family = nbinom1)
#diagnose(mod6)
testDispersion(mod6)
plot(simulateResiduals(mod6))

# negative binomial with age as offset, changing zero inflated formula
mod7 <- glmmTMB(offspring ~ Fgrm + sex + site,
                ziformula = ~., # same as full conditional formula
                offset = SVL_by_centered,
                data = data_flt, 
                family = nbinom1)
#diagnose(mod7)
testDispersion(mod7)
plot(simulateResiduals(mod7))

# negative binomial with age as offset, including birth year as proxy for search effort over time
mod8 <- glmmTMB(offspring ~ Fgrm + sex + site*SVL_by_centered,
                     ziformula = ~site*SVL_by_centered, 
                     offset = SVL_by_centered,
                     data = data_flt, 
                     family = nbinom1)
#diagnose(mod8)
testDispersion(mod8)
plot(simulateResiduals(mod8))

# negative binomial with age as offset, including birth year as proxy for search effort over time, NO FGRM
mod9 <- glmmTMB(offspring ~ sex + site*SVL_by_centered,
                ziformula = ~site*SVL_by_centered, 
                offset = SVL_by_centered,
                data = data_flt, 
                family = nbinom1)
#diagnose(mod8)
testDispersion(mod8)
plot(simulateResiduals(mod8))


mod.AIC = AIC(mod1,mod2,mod3,mod4,mod5, mod6, mod7, mod8)

mod.AIC[order(mod.AIC$AIC),]

mod.AIC = AIC(mod1,mod2,mod3,mod4,mod5, mod6, mod7, mod8, mod9)

mod.AIC[order(mod.AIC$AIC),"AIC"][1] - mod.AIC[order(mod.AIC$AIC),"AIC"]



### plots  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
site_palette <- met.brewer("Archambault", n = 6)

plot(data_flt$offspring~data_flt$Fgrm, col = alpha("black", alpha = 0.3), pch = 19, xlab = "Fgrm", ylab = "number of offspring")


data_flt <- data_flt[data_flt$Fgrm < 0.4,]
pdf(file = "../inbreeding_models/off_v_fgrm.pdf", width = 10, height = 5)
par(mfrow=c(1,2))
plot(data_flt[which(data_flt$site == "ELF"), "offspring"]~data_flt[which(data_flt$site == "ELF"), "Fgrm"], col = alpha(site_palette[2], alpha = 0.4), pch = 19, 
     xlab = "Fgrm", ylab = "number of offspring", main = "Cass")
plot(data_flt[which(data_flt$site == "PCC"), "offspring"]~data_flt[which(data_flt$site == "PCC"), "Fgrm"], col = alpha(site_palette[1], alpha = 0.4), pch = 19, 
     xlab = "Fgrm", ylab = "number of offspring", main = "Barry")
dev.off()


pdf(file = "../inbreeding_models/off_v_fgrm_poster.pdf", width = 10, height = 5)
par(mfrow=c(1,2))
cex_val = 2
plot(data_flt[which(data_flt$site == "ELF"), "offspring"]~data_flt[which(data_flt$site == "ELF"), "Fgrm"], col = alpha(site_palette[2], alpha = 0.4), pch = 19, 
     xlab = "Fgrm", ylab = "number of offspring", main = "Cass", cex = cex_val, cex.lab = cex_val, cex.axis = cex_val, cex.main =cex_val)
plot(data_flt[which(data_flt$site == "PCC"), "offspring"]~data_flt[which(data_flt$site == "PCC"), "Fgrm"], col = alpha(site_palette[1], alpha = 0.4), pch = 19, 
     xlab = "Fgrm", ylab = "number of offspring", main = "Barry", cex = cex_val, cex.lab = cex_val, cex.axis = cex_val, cex.main =cex_val)
dev.off()

### Inbreeding over time  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

plot(data_flt[which(data_flt$site == "ELF"), "Fgrm"]~data_flt[which(data_flt$site == "ELF"), "SVL_by"])
plot(data_flt[which(data_flt$site == "PCC"), "Fgrm"]~data_flt[which(data_flt$site == "PCC"), "SVL_by"])

### Non-parametric/permutation approach  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# LR version
plot(data_flt$offspring~ data_flt$Fgrm)
abline(lm(data_flt$offspring~ data_flt$Fgrm))
real_mod <- lm(data_flt$offspring~ data_flt$Fgrm)
real_coeff <- real_mod$coefficients[2]

# permute Fgrm values and re run
permuts <- 10000
permut_slopes <- vector(length = permuts)
for(i in 1:permuts){
  F_shuffle <- sample(data_flt$Fgrm)
  shuffle_mod <- lm(data_flt$offspring ~F_shuffle)
  permut_slopes[i] <- shuffle_mod$coefficients[2]
}
quants <- quantile(permut_slopes, probs = c(0.05, 0.95))

# plot 
pdf(file = "../inbreeding_models/permutation_dist_LR.pdf", width = 10, height = 5)
par(mfrow=c(1,2))
plot(data_flt$offspring~ data_flt$Fgrm, pch = 19, col = alpha("gray19", alpha = 0.5), 
     xlab = "Fgrm", 
     ylab = "number of offspring")
abline(real_mod, col = site_palette[4], lwd = 1.5)
hist(permut_slopes, main = NULL, xlab = "slope of offspring~inbreeding", 
     xlim = c(-10, 15), breaks = 50, col = "gray80", border = "gray80")
mtext("10,000 permutations")
abline(v = quants, col = "black", lty = 2)
abline(v = real_coeff, col = site_palette[4], lwd =1.5)
dev.off()

# poisson regresson
mod1 <- glm(formula = offspring ~ Fgrm, data = data_flt_old, family = 'poisson')
real_coeff <- mod1$coefficients[2]

permut_slopes <- vector(length = permuts)
for(i in 1:permuts){
  data_flt_old$F_shuffle <- sample(data_flt_old$Fgrm)
  shuffle_mod <- glm(formula = offspring ~ F_shuffle, data = data_flt_old, family = 'poisson')
  permut_slopes[i] <- shuffle_mod$coefficients[2]
}
quants <- quantile(permut_slopes, probs = c(0.05, 0.95))

# predict POIS regression 
newdat <- data.frame(Fgrm = seq(-0.2, 0.4, 0.1))
newdat$predPOIS <- predict(mod1, newdata = newdat, type = "response")

# plot 
pdf(file = "../inbreeding_models/permutation_dist_Pois.pdf", width = 10, height = 5)
par(mfrow=c(1,2))
plot(data_flt$offspring~ data_flt$Fgrm, pch = 19, col = alpha("gray19", alpha = 0.5), 
     xlab = "Fgrm", 
     ylab = "number of offspring")
lines(predPOIS ~ Fgrm, newdat, col = site_palette[4], lwd = 2)

hist(permut_slopes, main = NULL, xlab = "slope of offspring~inbreeding", 
     xlim = c(-10, 15), breaks = 50, col = "gray80", border = "gray80")
mtext("10,000 permutations")
abline(v = quants, col = "black", lty = 2)
abline(v = real_coeff, col = site_palette[4], lwd =2)
dev.off()


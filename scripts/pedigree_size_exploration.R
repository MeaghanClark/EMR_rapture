
# notes
# double check how I am dealing with multiple captures
library(stringr)
library(scales)
library(sequoia)
library(GGally)

date = "02082024"

# define custom functions ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

vb_growth_curve <- function(t, Loo, k = 0.378, Lo = 19.24){
  # from http://elasmollet.org/VBGF/VBGF.html
  return(Loo - (Loo - Lo)*exp(-k*t))
}

est_age_from_curve <- function(L, Loo, k = 0.378, Lo = 19.24){

  if(L < Loo){
    return(-(log((Loo-L)/(Loo-Lo))) / (k))
  }else{
    return(NULL) 
  }
}

plot_curves <- function(){
  x <- seq(0,20,length.out=20)
  # cass males
  plot(x,vb_growth_curve(t = x, Loo = 61.8, k = 0.378, Lo = 19.24),type='l',xlab="age",ylab="SVL (cm)",main="", 
       xlim = c(0, 20), 
       ylim = c(15, 65), 
       col = "black", 
       lty = 2)
  # cass females
  lines(x,vb_growth_curve(t = x, Loo = 58.7, k = 0.378, Lo = 19.24),type='l',xlab="",ylab="",main="", 
        col = "black")
  
  # barry males
  lines(x,vb_growth_curve(t = x, Loo = 63.4, k = 0.378, Lo = 19.24),type='l',xlab="",ylab="",main="", 
        col = "gray", 
        lty = 2)
  # barry females
  lines(x,vb_growth_curve(t = x, Loo = 60.3, k = 0.378, Lo = 19.24),type='l',xlab="",ylab="",main="", 
        col = "gray")
  legend("bottomright", legend = c("Cass Males", "Cass Females", "Barry Males", "Barry Females"), col = c("black", "black", "gray", "gray"), lty = c(2, 1, 2, 1))
}

plot_curves_with_CI <- function(asymptotic_L){
  x <- seq(0,20,length.out=20)
  # cass males
  plot(NULL,type='l',xlab="age",ylab="SVL (cm)",main="", 
       xlim = c(0, 20), 
       ylim = c(15, 65), 
       col = "black", 
       lty = 1)
  
  polygon(x = c(x, rev(x)), y = c(vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[2,"upper_CI"])), rev(vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[2,"lower_CI"])))), 
          border = NA, 
          col = alpha(alpha = 0.25, "thistle4"))
  lines(x,vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[2,"L_est"])))
  
  # cass females
  polygon(x = c(x, rev(x)), y = c(vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[1,"upper_CI"])), rev(vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[1,"lower_CI"])))), 
          border = NA, 
          col = alpha(alpha = 0.25, "thistle"))
  lines(x,vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[1,"L_est"])))
  
  # barry males
  polygon(x = c(x, rev(x)), y = c(vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[4,"upper_CI"])), rev(vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[4,"lower_CI"])))), 
          border = NA, 
          col = alpha(alpha = 0.25, "skyblue4"))
  lines(x,vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[4,"L_est"])))
  # barry females
  polygon(x = c(x, rev(x)), y = c(vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[3,"upper_CI"])), rev(vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[3,"lower_CI"])))), 
          border = NA, 
          col = alpha(alpha = 0.25, "skyblue"))
  lines(x,vb_growth_curve(t = x, Loo = as.numeric(asymptotic_L[3,"L_est"])))
  
  legend("bottomright", legend = c("Cass Males", "Cass Females", "Barry Males", "Barry Females"), col = c("thistle4", "thistle", "skyblue4", "skyblue"), pch = 15)
}

plot_pair_asym_L <- function(focal_ind, parent_ind, location, data){
  #focal_ind = ped.PCC.ped$Pedigree[i,"id"]
  #parent_ind = ped.PCC.ped$Pedigree[i,"dam"]
  
  par(mfrow = c(1,2))
  focal_SVL <- data[which(data$ID == focal_ind),"SVL"]
  parent_SVL <- data[which(data$ID == parent_ind),"SVL"]
  focal_year <- data[which(data$ID == focal_ind),"year"]
  parent_year <- data[which(data$ID == parent_ind),"year"]
  # size
  plot(x = c(focal_year, parent_year), 
       y = c(focal_SVL, parent_SVL), 
       xlab = "year", 
       ylab = "SVL", 
       xlim = c(2010, 2020), 
       #ylim = c(50, 55), 
       pch = 19, 
       col = c("red",  "blue"))
  
  if(data[which(data$ID == focal_ind),"sex"] == 1){
    focal_sex <- "female"
  }else if(data[which(data$ID == focal_ind),"sex"] == 2){
    focal_sex <- "male"
  }
  if(data[which(data$ID == parent_ind),"sex"] == 1){
    parent_sex <- "female"
  }else if(data[which(data$ID == parent_ind),"sex"] == 2){
    parent_sex <- "male"
  }
  
  focal_age <- est_age_from_curve(L = focal_SVL, Loo = subset(asymptotic_L, site == location & sex == focal_sex, select = L_est))
  parent_age <- est_age_from_curve(L = parent_SVL, Loo = subset(asymptotic_L, site == location & sex == parent_sex, select = L_est))
  
  focal_age_upper <- est_age_from_curve(L = focal_SVL, Loo = subset(asymptotic_L, site == location & sex == focal_sex, select = upper_CI))
  parent_age_upper <- est_age_from_curve(L = parent_SVL, Loo = subset(asymptotic_L, site == location & sex == parent_sex, select = upper_CI))
  focal_age_lower <- est_age_from_curve(L = focal_SVL, Loo = subset(asymptotic_L, site == location & sex == focal_sex, select = lower_CI))
  parent_age_lower <- est_age_from_curve(L = parent_SVL, Loo = subset(asymptotic_L, site == location & sex == parent_sex, select = lower_CI))

  if(is.null(parent_age)){
    parent_age = 20 # some very large number
  }
  if(is.null(parent_age_upper)){
    parent_age_upper = 20 # some very large number
  }
  if(is.null(parent_age_lower)){
    parent_age_lower = 20 # some very large number
  }
  
  if(is.null(focal_age)){
    focal_age = 20 # some very large number
  }
  if(is.null(focal_age_upper)){
    focal_age_upper = 20 # some very large number
  }
  if(is.null(focal_age_lower)){
    focal_age_lower = 20 # some very large number
  }
  
  plot(x = c(focal_year, parent_year), 
       y = c(focal_age, parent_age), 
       xlab = "year", 
       ylab = "age", 
       xlim = c(2010, 2022), 
       ylim = c(0, 20), 
       pch = 19, 
       col = c("red",  "blue"))
  
  lines(x = rep(focal_year, 2), y = c(focal_age_upper, focal_age_lower))
  lines(x = rep(parent_year, 2), y = c(parent_age_upper, parent_age_lower))
  
  latest_parent_age = parent_age_upper + (focal_year - parent_year)
  oldest_child_age = focal_age_lower
  
  points(x = focal_year, y = latest_parent_age, col = "blue")
  points(x = focal_year, y = oldest_child_age, col = "red")
  
  legend(x = "topright", legend = c("focal", "parent", "focal, upper est", "parent, proj"), pch = c(19, 19, 1, 1), col = c("red", "blue"))
  mtext(paste0("age difference in focal year ", latest_parent_age - oldest_child_age))
  
  # what do we want to return? 
    # if parent and offspring are inferred correctly (is the project parent above the focal upper est?)
    # latest_parent_age - oldest_child_age
  return(as.numeric(latest_parent_age - oldest_child_age))
}

eval_pair_SVL <- function(focal_ind, parent_ind, data){
  # focal_ind = focal_ind
  # parent_ind = dam
  # data = SVL_data
  # define site
  if(grepl("PCC", focal_ind)){
    site = "barry"
  }else if(grepl("ELF", focal_ind)){
    site = "cass"
  }
  if(!is.na(parent_ind)){
    age_dif <- plot_pair_asym_L(focal_ind, parent_ind, location = site, data) 
  }else{
    print("parent ind is NA")
    age_dif <- NA
    }
  
  return(age_dif)
}

infer_birth_year <- function(ind, location, data){
    ind_SVL <- data[which(data$ID == ind),"SVL"]
    ind_year <- data[which(data$ID == ind),"year"]

    if(data[which(data$ID == ind),"sex"] == 1){
      ind_sex <- "female"
    }else if(data[which(data$ID == ind),"sex"] == 2){
      ind_sex <- "male"
    }
    ind_age <- est_age_from_curve(L = ind_SVL, Loo = subset(asymptotic_L, site == location & sex == ind_sex, select = L_est))

    ind_age_upper <- est_age_from_curve(L = ind_SVL, Loo = subset(asymptotic_L, site == location & sex == ind_sex, select = upper_CI))
    ind_age_lower <- est_age_from_curve(L = ind_SVL, Loo = subset(asymptotic_L, site == location & sex == ind_sex, select = lower_CI))

    if(is.null(ind_age)){
      ind_age = NA # some very large number
    }
    if(is.null(ind_age_upper)){
      ind_age_upper = NA  # some very large number
    }
    if(is.null(ind_age_lower)){
      ind_age_lower = NA # some very large number
    }
    
    return(c(as.numeric(ind_year - ind_age), as.numeric(ind_year - ind_age_upper), as.numeric(ind_year - ind_age_lower)))
}



infer_birth_year_wSVL <- function(ind, location, ind_SVL, ind_year, ind_sex){
  # ind = PCC_known_by_inds[i]
  # location = "barry"
  # ind_SVL = as.numeric(PCC_exp_metaData[[PCC_known_by_inds[i]]][n,"SVL"])
  # ind_year =PCC_exp_metaData[[PCC_known_by_inds[i]]][n,"year"] 
  # ind_sex = 1
  
  ind_year <- as.numeric(ind_year)
  if(ind_sex == 1){
    ind_sex = "female"
  }else if(ind_sex == 2){
    ind_sex = "male"
  }
  ind_age <- est_age_from_curve(L = ind_SVL, Loo = subset(asymptotic_L, site == location & sex == ind_sex, select = L_est))
  ind_age_upper <- est_age_from_curve(L = ind_SVL, Loo = subset(asymptotic_L, site == location & sex == ind_sex, select = upper_CI))
  ind_age_lower <- est_age_from_curve(L = ind_SVL, Loo = subset(asymptotic_L, site == location & sex == ind_sex, select = lower_CI))
  
  if(is.null(ind_age)){
    ind_age = NA # some very large number
  }
  if(is.null(ind_age_upper)){
    ind_age_upper = NA  # some very large number
  }
  if(is.null(ind_age_lower)){
    ind_age_lower = NA # some very large number
  }
  
  return(c(as.numeric(ind_year - ind_age), as.numeric(ind_year - ind_age_upper), as.numeric(ind_year - ind_age_lower)))
}

# ------------------------------------------------------------------------------------------------------------

# load data ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

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

# "extra" metadata for analyses
load(paste0("../pedigree_exploration/ELF_expanded_metaData.Robj"), verbose = T) # ELF_exp_metaData
load(paste0("../pedigree_exploration/PCC_expanded_metaData.Robj"), verbose = T) # PCC_exp_metaData
PCC_exp_metaData[["PCC_42"]]$year[1] <- "2015"

load(paste0("../pedigree_reconstruction/PCC_pedigree_results_", "05072024", ".Robj"), verbose = T)
PCC_gt_for_ped_drop2 <- PCC_results[[1]]
PCC_meta <- PCC_results[[2]]
ped.PCC.ped <- PCC_results[[4]]
error_rate <- PCC_results[[5]]

load(paste0("../pedigree_reconstruction/ELF_pedigree_results_", "05072024", ".Robj"), verbose = T)
ELF_gt_for_ped_drop2 <- ELF_results[[1]]
ELF_meta <- ELF_results[[2]]
ped.ELF.ped <- ELF_results[[4]]


# generate SVL data frames

# PCC
PCC_SVL_data <- data.frame(matrix(NA, nrow = length(ped.PCC.ped$Pedigree$id), ncol = 4))
colnames(PCC_SVL_data) <- c("ID", "sex", "year", "SVL")

# get smallest size! 
for(i in 1:length(ped.PCC.ped$Pedigree$id)){ # for each individual in the pedigree
  
  # get id of focal individual 
  focal_ind <- ped.PCC.ped$Pedigree$id[i]
  focal_ind_search <- str_remove(ped.PCC.ped$Pedigree$id[i], "_P[0-9]")
  PCC_SVL_data[i, "ID"] <- focal_ind_search
  if(grepl("_", focal_ind)){
    # identify minimum SVL 
    focal_SVL <- min(as.numeric(PCC_exp_metaData[[focal_ind_search]]$SVL))
    if(!is.na(focal_SVL)){
      # find sex of focal individual
      PCC_SVL_data$sex[i] <- PCC_meta[which(PCC_meta$ID == ped.PCC.ped$Pedigree$id[i]),"Sex"]
      PCC_SVL_data$year[i] <- PCC_exp_metaData[[focal_ind_search]]$year[which(as.numeric(PCC_exp_metaData[[focal_ind_search]]$SVL) == focal_SVL)] 
      PCC_SVL_data$SVL[i] <- focal_SVL
    }
  }
}

PCC_SVL_data$SVL <- as.numeric(PCC_SVL_data$SVL)
PCC_SVL_data$year <- as.numeric(PCC_SVL_data$year)

#ELF
ELF_SVL_data <- data.frame(matrix(NA, nrow = length(ped.ELF.ped$Pedigree$id), ncol = 4))
colnames(ELF_SVL_data) <- c("ID", "sex", "year", "SVL")

# get smallest size! 
for(i in 1:length(ped.ELF.ped$Pedigree$id)){ # for each individual in the pedigree
  
  # get id of focal individual 
  focal_ind <- ped.ELF.ped$Pedigree$id[i]
  focal_ind_search <- str_remove(ped.ELF.ped$Pedigree$id[i], "_P[0-9]")
  ELF_SVL_data[i, "ID"] <- focal_ind_search
  if(grepl("_", focal_ind)){
    # identify minimum SVL 
    focal_SVL <- min(as.numeric(ELF_exp_metaData[[focal_ind_search]]$SVL), na.rm = TRUE)
    if (min(as.numeric(ELF_exp_metaData[[focal_ind_search]]$SVL), na.rm = TRUE) == Inf){
      focal_SVL <- NA
    }
    if(!is.na(focal_SVL)){
      # find sex of focal individual
      ELF_SVL_data$sex[i] <- ELF_meta[which(ELF_meta$ID == ped.ELF.ped$Pedigree$id[i]),"Sex"]
      ELF_SVL_data$year[i] <- ELF_exp_metaData[[focal_ind_search]]$year[which(as.numeric(ELF_exp_metaData[[focal_ind_search]]$SVL) == focal_SVL)] 
      ELF_SVL_data$SVL[i] <- focal_SVL
    }
  }
}

ELF_SVL_data$SVL <- as.numeric(ELF_SVL_data$SVL)
ELF_SVL_data$year <- as.numeric(ELF_SVL_data$year)

# ------------------------------------------------------------------------------------------------------------

# 1. test to see if asym_L test works. For each individual in the pedigree, compare size to their dam and sire, and check if we called the correct relationship
# ------------------------------------------------------------------------------------------------------------

# PCC # -----------------------------------------------------------------------------------------------------------
PCC_pedigree_eval <- ped.PCC.ped$Pedigree
PCC_pedigree_eval$dam_age_dif <- NA
PCC_pedigree_eval$sire_age_dif <- NA

# for each individual in the pedigree... 
for(i in 1:nrow(ped.PCC.ped$Pedigree)){
  focal_ind <- ped.PCC.ped$Pedigree[i,"id"]
  dam <- ped.PCC.ped$Pedigree[i,"dam"]
  sire <- ped.PCC.ped$Pedigree[i,"sire"]
  # who has metadata? 
  
  if(length(PCC_SVL_data[which(PCC_SVL_data$ID == focal_ind),"SVL"]) == 1 
     && length(PCC_SVL_data[which(PCC_SVL_data$ID == dam),"SVL"]) == 1 
     && !is.na(PCC_SVL_data[which(PCC_SVL_data$ID == dam),"SVL"]) 
     && !is.na(PCC_SVL_data[which(PCC_SVL_data$ID == focal_ind),"SVL"]) 
     && !(PCC_SVL_data[which(PCC_SVL_data$ID == focal_ind),"sex"] == -9)
     && !(PCC_SVL_data[which(PCC_SVL_data$ID == dam),"sex"] == -9)) { # if dam and focal ind have SVL data
    
    # EVAL FUNC
    PCC_pedigree_eval$dam_age_dif[i] <- eval_pair_SVL(focal_ind, dam, PCC_SVL_data)
    
  }else(PCC_pedigree_eval$dam_age_dif[i] <- NA)
  if(length(PCC_SVL_data[which(PCC_SVL_data$ID == focal_ind),"SVL"]) == 1 
     && length(PCC_SVL_data[which(PCC_SVL_data$ID == sire),"SVL"]) == 1 
     && !is.na(PCC_SVL_data[which(PCC_SVL_data$ID == sire),"SVL"]) 
     && !is.na(PCC_SVL_data[which(PCC_SVL_data$ID == focal_ind),"SVL"]) 
     && !(PCC_SVL_data[which(PCC_SVL_data$ID == focal_ind),"sex"] == -9)
     && !(PCC_SVL_data[which(PCC_SVL_data$ID == sire),"sex"] == -9)){ # if sire and focal individual have SVL data
    # EVAL FUNC
    PCC_pedigree_eval$sire_age_dif[i] <- eval_pair_SVL(focal_ind, sire, PCC_SVL_data)
  }else(PCC_pedigree_eval$sire_age_dif[i] <- NA)
}

hist(PCC_pedigree_eval$dam_age_dif, main = "PCC", xlab = "inferred age difference between focal and dam", breaks = 25)
hist(PCC_pedigree_eval$sire_age_dif, main = NULL, xlab = "inferred age difference between focal and sire", breaks = 25)


# ELF # ------------------------------------------------------------------------------------------------------------

ELF_pedigree_eval <- ped.ELF.ped$Pedigree
ELF_pedigree_eval$dam_age_dif <- NA
ELF_pedigree_eval$sire_age_dif <- NA

# for each individual in the pedigree... 
for(i in 1:nrow(ped.ELF.ped$Pedigree)){
  focal_ind <- ped.ELF.ped$Pedigree[i,"id"]
  dam <- ped.ELF.ped$Pedigree[i,"dam"]
  sire <- ped.ELF.ped$Pedigree[i,"sire"]
  # who has metadata? 
  
  if(length(ELF_SVL_data[which(ELF_SVL_data$ID == focal_ind),"SVL"]) == 1 && length(ELF_SVL_data[which(ELF_SVL_data$ID == dam),"SVL"]) == 1 && !is.na(ELF_SVL_data[which(ELF_SVL_data$ID == dam),"SVL"]) && !is.na(ELF_SVL_data[which(ELF_SVL_data$ID == focal_ind),"SVL"])) { # if dam and focal ind have SVL data
    # EVAL FUNC
    ELF_pedigree_eval$dam_age_dif[i] <- eval_pair_SVL(focal_ind, dam, ELF_SVL_data)
    
  }else(ELF_pedigree_eval$dam_age_dif[i] <- NA)
  if(length(ELF_SVL_data[which(ELF_SVL_data$ID == focal_ind),"SVL"]) == 1 && length(ELF_SVL_data[which(ELF_SVL_data$ID == sire),"SVL"]) == 1 && !is.na(ELF_SVL_data[which(ELF_SVL_data$ID == sire),"SVL"]) && !is.na(ELF_SVL_data[which(ELF_SVL_data$ID == focal_ind),"SVL"])){ # if sire and focal individual have SVL data
    # EVAL FUNC
    ELF_pedigree_eval$sire_age_dif[i] <- eval_pair_SVL(focal_ind, sire, ELF_SVL_data)
  }else(ELF_pedigree_eval$sire_age_dif[i] <- NA)
}


hist(ELF_pedigree_eval$dam_age_dif, main = "ELF", xlab = "inferred age difference between focal and dam", breaks = 25)
hist(ELF_pedigree_eval$sire_age_dif, main = NULL, xlab = "inferred age difference between focal and sire", breaks = 25)


# how do we feel about instances of negative age differences?

PCC_pedigree_eval[which(PCC_pedigree_eval$dam_age_dif < 0),]
# PCC_196 and PCC_182
PCC_exp_metaData[["PCC_196"]]
PCC_exp_metaData[["PCC_182"]] # error in SVL measurements, confident to attribute this neg. value to error in SVL by inference

# PCC_272 and PCC_151
PCC_exp_metaData[["PCC_272"]]
PCC_exp_metaData[["PCC_151"]] # both individuals are adults captured once as near asymptotic size 

ELF_pedigree_eval[which(ELF_pedigree_eval$dam_age_dif < 0),]
# ELF_788 and ELF_67
ELF_exp_metaData[["ELF_788"]]
ELF_exp_metaData[["ELF_67"]] # both captured for the first time at/near asymptotic size 

# ELF_846 and ELF_65
ELF_exp_metaData[["ELF_846"]] # focal
ELF_exp_metaData[["ELF_65"]] # dam, perhaps some error around SVL from multiple capture events? 
infer_birth_year_wSVL("ELF_65", "cass", 52.1, 2010, 1) # 2005
infer_birth_year_wSVL("ELF_65", "cass", 56.5, 2012, 1) # 2004
infer_birth_year_wSVL("ELF_846", "cass", 58.5, 2016, 1) # 2002
# both > 50 cm when first captured, appears to be some error around SVL from different capture events of ELF_65, 

# ELF_930 and ELF_187
ELF_exp_metaData[["ELF_930"]] # focal
ELF_exp_metaData[["ELF_187"]] # dam 
# both over 50 cm when captured

ELF_pedigree_eval[which(ELF_pedigree_eval$sire_age_dif < 0),]
# ELF_153 and ELF_205
ELF_exp_metaData[["ELF_153"]] # focal
ELF_exp_metaData[["ELF_205"]] # sire
# both over 50 SVL when captured

# ELF_321 and ELF_130
ELF_exp_metaData[["ELF_321"]] # focal
ELF_exp_metaData[["ELF_130"]] # sire
# both over 50 SVL when captured


# Plot ------------------------------------------------------------------------------------------------------------

par(mfrow=c(2,2))
hist(PCC_pedigree_eval$dam_age_dif, main = "PCC", xlab = "inferred age dif, focal and dam", breaks = 25)
abline(v = 0, col = "red", lty = 2)
hist(PCC_pedigree_eval$sire_age_dif, main = NULL, xlab = "inferred age dif, focal and sire", breaks = 25)
abline(v = 0, col = "red", lty = 2)

hist(ELF_pedigree_eval$dam_age_dif, main = "ELF", xlab = "inferred age dif, focal and dam", breaks = 25)
abline(v = 0, col = "red", lty = 2)
hist(ELF_pedigree_eval$sire_age_dif, main = NULL, xlab = "inferred age dif, focal and sire", breaks = 25)
abline(v = 0, col = "red", lty = 2)


# ------------------------------------------------------------------------------------------------------------


# Of individuals with known birth years... do we infer those years using SVL data?  ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# get all SVL recorded for individuals to assess relationship between size and accuracy ------------------------------------------------------------------------------------------------------------

PCC_known_by_inds <- unique(subset(x = PCC_meta, subset = !(BirthYear == -9), select = ID))$ID

by_test <- vector()
for(i in 1:length(PCC_known_by_inds)){
  # how many captures? 
  no_caps <- dim(PCC_exp_metaData[[PCC_known_by_inds[i]]])[1]
  known_by <- PCC_meta[which(PCC_meta$ID == PCC_known_by_inds[i]), "BirthYear"]
  sex <-  PCC_meta[which(PCC_meta$ID == PCC_known_by_inds[i]), "Sex"]
  tryCatch(
    {
      for(n in 1:no_caps){
        SVL <- as.numeric(PCC_exp_metaData[[PCC_known_by_inds[i]]][n,"SVL"])
        inferred_by <- infer_birth_year_wSVL(ind = PCC_known_by_inds[i], location = "barry", ind_SVL = as.numeric(PCC_exp_metaData[[PCC_known_by_inds[i]]][n,"SVL"]), ind_year =PCC_exp_metaData[[PCC_known_by_inds[i]]][n,"year"], ind_sex = sex)[1]
        inferred_by_upper <- infer_birth_year(PCC_known_by_inds[i], location = "barry", PCC_SVL_data)[2]
        inferred_by_lower <- infer_birth_year(PCC_known_by_inds[i], location = "barry", PCC_SVL_data)[3]
        data_vec <- c(PCC_known_by_inds[i], "PCC", sex, as.integer(known_by), as.numeric(SVL), as.numeric(inferred_by), as.numeric(inferred_by_upper), as.numeric(inferred_by_lower))
        by_test <- rbind.data.frame(by_test, data_vec)
      }
    },
    error = function(e){
      cat("Error occured at interation", i, ": ", conditionMessage(e), "\n")
    }
  )
}

colnames(by_test) <- c("ID", "Site", "Sex", "BirthYear", "SVL", "inferred_by", "inferred_by_upper", "inferred_by_lower")

# add ELF
ELF_known_by_inds <- unique(subset(x = ELF_meta, subset = !(BirthYear == -9), select = ID))$ID

for(i in 1:length(ELF_known_by_inds)){
  # how many captures? 
  no_caps <- dim(ELF_exp_metaData[[ELF_known_by_inds[i]]])[1]
  known_by <- ELF_meta[which(ELF_meta$ID == ELF_known_by_inds[i]), "BirthYear"]
  sex <- ELF_meta[which(ELF_meta$ID == ELF_known_by_inds[i]), "Sex"]
  tryCatch(
    {
      for(n in 1:no_caps){
        SVL <- as.numeric(ELF_exp_metaData[[ELF_known_by_inds[i]]][n,"SVL"])
        inferred_by <- infer_birth_year_wSVL(ind = ELF_known_by_inds[i], location = "cass", ind_SVL = as.numeric(ELF_exp_metaData[[ELF_known_by_inds[i]]][n,"SVL"]), ind_year =ELF_exp_metaData[[ELF_known_by_inds[i]]][n,"year"], ind_sex = sex)[1]
        inferred_by_upper <- infer_birth_year(ELF_known_by_inds[i], location = "cass", ELF_SVL_data)[2]
        inferred_by_lower <- infer_birth_year(ELF_known_by_inds[i], location = "cass", ELF_SVL_data)[3]
        data_vec <- c(ELF_known_by_inds[i], "ELF", sex, as.integer(known_by), as.numeric(SVL), as.numeric(inferred_by), as.numeric(inferred_by_upper), as.numeric(inferred_by_lower))
        by_test <- rbind.data.frame(by_test, data_vec)
      }
    },
    error = function(e){
      cat("Error occured at interation", i, ": ", conditionMessage(e), "\n")
    }
  )
}

by_test$BirthYear <- as.integer(by_test$BirthYear)
by_test$inferred_by <- as.numeric(by_test$inferred_by)
by_test$inferred_by_upper <- as.numeric(by_test$inferred_by_upper)
by_test$inferred_by_lower <- as.numeric(by_test$inferred_by_lower)


# error in birth year 
by_test$by_error <- by_test$BirthYear - by_test$inferred_by
by_test$by_error_upper <- by_test$BirthYear - by_test$inferred_by_upper
by_test$by_error_lower <- by_test$BirthYear - by_test$inferred_by_lower
by_test$by_error_round <- by_test$BirthYear - round(by_test$inferred_by)


colors <- by_test$Site
colors[which(by_test$Site == "PCC")] <- "green"
colors[which(by_test$Site == "ELF")] <- "blue"

plot(x = by_test$BirthYear, y = by_test$inferred_by, col = alpha(colors, alpha = 0.25), pch = 19, xlab = "known birth year", ylab = "inferred birht year")
abline(a = 0, b = 1)
abline(a = 1, b = 1, lty = 2)
abline(a = -1, b = 1, lty = 2)

SVL_palette = "OKeeffe2"
SVL_colors <- met.brewer(SVL_palette, n = length(unique(by_test$SVL)))
palette(SVL_colors)

plot(x = by_test$BirthYear, y = by_test$inferred_by, col = as.factor(by_test$SVL), pch = 19, xlab = "known birth year", ylab = "inferred birht year")
abline(a = 0, b = 1)
abline(a = 1, b = 1, lty = 2)
abline(a = -1, b = 1, lty = 2)


colors <- by_test[,c("Site", "Sex")]
colors[which(by_test$Site == "PCC" & by_test$Sex == 1),"color"] <- "palegreen3"
colors[which(by_test$Site == "PCC" & by_test$Sex == 2),"color"] <- "darkgreen"
colors[which(by_test$Site == "ELF" & by_test$Sex == 1),"color"] <- "tomato"
colors[which(by_test$Site == "ELF" & by_test$Sex == 2),"color"] <- "violetred3"

plot(x = by_test$SVL, y = by_test$by_error, xlab = "SVL at capture", ylab = "known birth year - inferred birth year", 
     col = alpha(colors$color, 0.25), 
     pch = 19)
abline(h = 0, col = "red", lty = 2)
abline(h = 2, col = "gray", lty = 2)
abline(h = -2, col = "gray", lty = 2)

abline(v = 45, col = "black", lty = 2)

legend("topleft", legend = c("PCC females", "PCC males", "ELF females", "ELF males"), pch = 19, col = c("palegreen3", "darkgreen", "tomato", "violetred3"))

# rounded
plot(x = by_test$SVL, y = by_test$by_error_round, xlab = "SVL at capture", ylab = "known birth year - rounded inferred birth year", 
     col = alpha(colors$color, 0.25), 
     pch = 19)
abline(h = 0, col = "red", lty = 2)
abline(h = 2, col = "gray", lty = 2)
abline(h = -2, col = "gray", lty = 2)

abline(v = 45, col = "black", lty = 2)

legend("topleft", legend = c("PCC females", "PCC males", "ELF females", "ELF males"), pch = 19, col = c("palegreen3", "darkgreen", "tomato", "violetred3"))

bound = 45

hist(by_test$by_error)
#hist(by_test$by_error[which(by_test$Sex == 1)], add = T, col = "tomato")
#hist(by_test$by_error[which(by_test$Sex == 2)], add = T, col = "violetred3")
# hist(by_test$by_error[which(by_test$SVL <= bound)], add = T, col = "turquoise")
hist(by_test$by_error[which(by_test$SVL > bound)], add = T, col = "skyblue")



by_test$ID


#### original: 

PCC_known_by <- subset(x = PCC_meta, subset = !(BirthYear == -9))
# deal with duplicate
for(i in 1:nrow(PCC_known_by)){
  tryCatch(
    {
      if(length(which(PCC_SVL_data$ID == PCC_known_by[i,"ID"])) == 0){
        PCC_known_by$SVL[i] <- NA
        PCC_known_by$inferred_by[i] <- NA
        PCC_known_by$inferred_by_upper[i] <- NA
        PCC_known_by$inferred_by_lower[i] <- NA
      }else{
        PCC_known_by$SVL[i] <- PCC_SVL_data[which(PCC_SVL_data$ID == PCC_known_by[i,"ID"]), "SVL"]
        PCC_known_by$inferred_by[i] <- infer_birth_year(PCC_known_by[i,"ID"], location = "barry", PCC_SVL_data)[1]
        PCC_known_by$inferred_by_upper[i] <- infer_birth_year(PCC_known_by[i,"ID"], location = "barry", PCC_SVL_data)[2]
        PCC_known_by$inferred_by_lower[i] <- infer_birth_year(PCC_known_by[i,"ID"], location = "barry", PCC_SVL_data)[3]
      }
    },
    error = function(e){
      cat("Error occured at interation", i, ": ", conditionMessage(e), "\n")
    }
  )
}

plot(NULL, xlim = c(2012,2020), 
     ylim = c(2012,2020), 
     xlab = "known birth year", 
     ylab = "inferred birth year")
points(x = PCC_known_by$BirthYear, y = round(PCC_known_by$inferred_by), pch = 19, col = alpha("black", alpha = 0.25))

abline(a = 0, b = 1, lty = 2 )


ELF_known_by <- subset(x = ELF_meta, subset = !(BirthYear == -9))
ELF_known_by <- ELF_known_by[!grepl("^\\w{3}_\\d{1,3}(_\\w{2,3})$", ELF_known_by$ID), ]
for(i in 1:nrow(ELF_known_by)){
  tryCatch(
    {
      
      if(is.na(ELF_SVL_data[which(ELF_SVL_data$ID == ELF_known_by[i,"ID"]), "SVL"]) | (ELF_SVL_data[which(ELF_SVL_data$ID == ELF_known_by[i,"ID"]),"sex"] == -9)){
        ELF_known_by$SVL[i] <- NA
        ELF_known_by$inferred_by[i] <- NA
        ELF_known_by$inferred_by_upper[i] <- NA
        ELF_known_by$inferred_by_lower[i] <- NA
      }else{
        ELF_known_by$SVL[i] <- ELF_SVL_data[which(ELF_SVL_data$ID == ELF_known_by[i,"ID"]), "SVL"]
        ELF_known_by$inferred_by[i] <- infer_birth_year(ELF_known_by[i,"ID"], location = "cass", ELF_SVL_data)[1]
        ELF_known_by$inferred_by_upper[i] <- infer_birth_year(ELF_known_by[i,"ID"], location = "cass", ELF_SVL_data)[2]
        ELF_known_by$inferred_by_lower[i] <- infer_birth_year(ELF_known_by[i,"ID"], location = "cass", ELF_SVL_data)[3]
      }
    },
    error = function(e){
      cat("Error occured at interation", i, ": ", conditionMessage(e), "\n")
    }
  )
}

plot(NULL, xlim = c(2012,2020), 
     ylim = c(2012,2020), 
     xlab = "known birth year", 
     ylab = "inferred birth year")
points(x = ELF_known_by$BirthYear, y = round(ELF_known_by$inferred_by), pch = 19, col = alpha("black", alpha = 0.25))
abline(a = 0, b = 1, lty = 2)

plot(NULL, xlim = c(2012,2020), 
     ylim = c(2012,2020), 
     xlab = "known birth year", 
     ylab = "inferred birth year")
points(x = ELF_known_by$BirthYear, y = ELF_known_by$inferred_by, pch = 19, col = alpha("green", alpha = 0.25))
points(x = PCC_known_by$BirthYear, y = PCC_known_by$inferred_by, pch = 19, col = alpha("blue", alpha = 0.25))

abline(a = 0, b = 1, lty = 2)
abline(a = 2, b = 1, lty = 2)
abline(a = -2, b = 1, lty = 2)
abline(a = 1, b = 1, lty = 2)
abline(a = -1, b = 1, lty = 2)


plot(NULL, xlim = c(-2, 2), ylim = c(10, 65),
     xlab = "known - inferred birth year", 
     ylab = "SVL")
points(x = as.numeric(ELF_known_by$BirthYear) - ELF_known_by$inferred_by, y = ELF_known_by$SVL, pch = 19, col = alpha("green", alpha = 0.25))
points(x = as.numeric(PCC_known_by$BirthYear) - PCC_known_by$inferred_by, y = PCC_known_by$SVL, pch = 19, col = alpha("blue", alpha = 0.25))


plot(NULL, xlim = c(-2, 2), ylim = c(-0.5, 0.5),
     xlab = "known - inferred birth year", 
     ylab = "inferred upper - lower")
points(x = as.numeric(ELF_known_by$BirthYear) - ELF_known_by$inferred_by, y = ELF_known_by$inferred_by_upper - ELF_known_by$inferred_by_lower, pch = 19, col = alpha("green", alpha = 0.25))
points(x = as.numeric(PCC_known_by$BirthYear) - PCC_known_by$inferred_by, y = PCC_known_by$inferred_by_upper - PCC_known_by$inferred_by_lower, pch = 19, col = alpha("blue", alpha = 0.25))


# AF by birth year ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
dim(ELF_results[[1]])
ELF_gt <- ELF_results[[1]]
ELF_gt[which(ELF_gt == -9)] <- NA
by_for_gt <- vector(length = nrow(ELF_gt))
for(i in 1:nrow(ELF_gt)){
  tryCatch(
    {
      by_for_gt[i] <- infer_birth_year(rownames(ELF_gt)[i], "cass", data = ELF_SVL_data)[1]
    },
    error = function(e){
      cat("Error occured at interation", i, ": ", conditionMessage(e), "\n")
      by_for_gt[i] <- NA
    }
  )
}

by_for_gt[which(by_for_gt ==0)] <- NA
by_for_gt <- round(by_for_gt)

by_for_gt_unique <- na.omit(unique(by_for_gt)[order(unique(by_for_gt))])

time_AF <- data.frame(matrix(nrow = ncol(ELF_gt), ncol = length(by_for_gt_unique)))
for(s in 1:ncol(ELF_gt)){
  for(y in 1:length(unique(by_for_gt_unique))){
    indices <- which(by_for_gt == by_for_gt_unique[y])
    time_AF[s,y] <- sum(ELF_gt[indices,s], na.rm = TRUE)/(2*length(na.omit(ELF_gt[indices,s]))) # need to handle NAs too 
  }
}

plot(NULL, xlim =c(1995, 2019), ylim = c(0,1), 
     xlab = "year", 
     ylab = "allele frequency")
for(s in 1:nrow(time_AF)){
  lines(x = by_for_gt_unique, y = time_AF[s,], col = alpha(1:10, alpha = 0.05))
  
}


# Apply to mystery pairs ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# for individuals in ambiguous PO pairs... 

# add by min and max relfecting their inferred ages +/- 1


# PCC

PCC_maybe_rel_ped <- GetMaybeRel(GenoM = PCC_gt_for_ped_drop2, LifeHistData=PCC_meta, Pedigree = ped.PCC.ped$Pedigree, Module = "ped", Err = error_rate)
PCC_maybe_po <- PCC_maybe_rel_ped$MaybeRel[PCC_maybe_rel_ped$MaybeRel$TopRel == "PO",]

PCC_PO_inds <- unique(c(PCC_maybe_po$ID1, PCC_maybe_po$ID2))

PCC_PO_by <- data.frame(matrix(data = NA, nrow = length(PCC_PO_inds), ncol = 2))
colnames(PCC_PO_by) <- c("ID", "by")
PCC_PO_by$ID <- PCC_PO_inds
for(i in 1:length(PCC_PO_inds)){
  tryCatch(
    {
      results <- infer_birth_year(PCC_PO_by$ID[i], location = "barry", PCC_SVL_data)
      PCC_PO_by$by[i] <- results[1]
      
    }, 
      error = function(e){
      cat("Error occured at interation", i, ": ", conditionMessage(e), "\n")
      PCC_PO_by$by[i] <- NA
      }
  )
}

PCC_PO_by$by.min <- round(PCC_PO_by$by) - 1 
PCC_PO_by$by.max <- round(PCC_PO_by$by) + 1 

date = "04132024"
save(PCC_PO_by, file = paste0("../pedigree_reconstruction/PCC_inferred_by_",date,".Robj"))


# ELF

ELF_maybe_rel_ped <- GetMaybeRel(GenoM = ELF_gt_for_ped_drop2, LifeHistData=ELF_meta, Pedigree = ped.ELF.ped$Pedigree, Module = "ped", Err = error_rate)
ELF_maybe_po <- ELF_maybe_rel_ped$MaybeRel[ELF_maybe_rel_ped$MaybeRel$TopRel == "PO",]

ELF_PO_inds <- unique(c(ELF_maybe_po$ID1, ELF_maybe_po$ID2))

ELF_PO_by <- data.frame(matrix(data = NA, nrow = length(ELF_PO_inds), ncol = 2))
colnames(ELF_PO_by) <- c("ID", "by")
ELF_PO_by$ID <- ELF_PO_inds
for(i in 1:length(ELF_PO_inds)){
  tryCatch(
    {
      results <- infer_birth_year(ELF_PO_by$ID[i], location = "barry", ELF_SVL_data)
      ELF_PO_by$by[i] <- results[1]
      
    }, 
    error = function(e){
      cat("Error occured at interation", i, ": ", conditionMessage(e), "\n")
      ELF_PO_by$by[i] <- NA
    }
  )
}

ELF_PO_by$by.min <- round(ELF_PO_by$by) - 1 
ELF_PO_by$by.max <- round(ELF_PO_by$by) + 1 

save(ELF_PO_by, file = paste0("../pedigree_reconstruction/ELF_inferred_by_",date,".Robj"))


# ------------------------------------------------------------------------------------------------------------









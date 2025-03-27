### Estimate individual longevity and extract birth year information from survival probabilities 
# ------------------------------------------------------------------------------------------------------------

library(stringr)
library(scales)
# load cumulative survival probabilities from Danielle's code output
# ------------------------------------------------------------------------------------------------------------

ELF_survival_data <- read.csv("../Massasauga survivorship/ELF_cumulative_survival.csv", header = T)[,-c(1,2)]
colnames(ELF_survival_data) <- str_remove(colnames(ELF_survival_data), "X")
PCC_survival_data <- read.csv("../Massasauga survivorship/PCCI_cumulative_survival.csv", header = T)[,-c(1,2)]
colnames(PCC_survival_data) <- str_remove(colnames(PCC_survival_data), "X")
# ------------------------------------------------------------------------------------------------------------

# estimate longevity using rounding rule
# ------------------------------------------------------------------------------------------------------------
ELF_death_years <- vector(length = nrow(ELF_survival_data))
# determine first year where survival probability is  0
  # deal with cases where length = 0
for(i in 1:nrow(ELF_survival_data)){
  rounded_probs <- round(ELF_survival_data[i,5:39], digits = 2)
  if(min(rounded_probs, na.rm = TRUE) > 0){
    ELF_death_years[i] <- NA
  }else{
    ELF_death_years[i] <- as.numeric(colnames(ELF_survival_data)[5:39][head(which(rounded_probs == 0), n = 1)])
  }
}

hist(as.numeric(ELF_death_years))

PCC_death_years <- vector(length = nrow(PCC_survival_data))
for(i in 1:nrow(PCC_survival_data)){
  rounded_probs <- round(PCC_survival_data[i,5:37], digits = 2)
  if(min(rounded_probs, na.rm = TRUE) > 0){
    PCC_death_years[i] <- NA
  }else{
    PCC_death_years[i] <- as.numeric(colnames(PCC_survival_data)[5:37][head(which(rounded_probs == 0), n = 1)])
  }
}


# ------------------------------------------------------------------------------------------------------------

# Find birth years from spreadsheet
# ------------------------------------------------------------------------------------------------------------
ELF_birth_years <- vector(length = nrow(ELF_survival_data))
for(i in 1:nrow(ELF_survival_data)){
  ELF_birth_years[i] <- colnames(ELF_survival_data)[5:39][head(which(!is.na(ELF_survival_data[i,5:39])), n = 1) -1]
}

PCC_birth_years <- vector(length = nrow(PCC_survival_data))
for(i in 1:nrow(PCC_survival_data)){
  PCC_birth_years[i] <- colnames(PCC_survival_data)[5:37][head(which(!is.na(PCC_survival_data[i,5:37])), n = 1) -1]
}

# ------------------------------------------------------------------------------------------------------------

# Make data frame with birth years and longevity 
# ------------------------------------------------------------------------------------------------------------

# use 2019 as last year DNA was taken from ELF, and 2021 from PCC 

ELF_dat <- cbind.data.frame(ELF_birth_years, ELF_death_years)
colnames(ELF_dat) <- c("estBirthYear", "estLastYear")
PCC_dat <- cbind.data.frame(PCC_birth_years, PCC_death_years)
colnames(PCC_dat) <- c("estBirthYear", "estLastYear")

ELF_dat$estLastYearPed <- ELF_dat$estLastYear
ELF_dat[which(is.na(ELF_dat$estLastYear)),"estLastYearPed"] <- 2019 # if our survival probability rule indicates an individual was alive at the end of the study, assign last year samples were collected for DNA extraction
ELF_dat[which(ELF_dat$estLastYear > 2019),"estLastYearPed"] <- 2019 #  # if our survival probability rule indicates an individual was likely died after the last year samples were collected for DNA extraction, last year samples were collected for DNA extraction
ELF_dat$longevity <- as.numeric(ELF_dat$estLastYear) - as.numeric(ELF_dat$estBirthYear)
ELF_dat$yearsContribute2Pedigree <- as.numeric(ELF_dat$estLastYearPed) - as.numeric(ELF_dat$estBirthYear)

PCC_dat$estLastYearPed <- PCC_dat$estLastYear
PCC_dat[which(is.na(PCC_dat$estLastYear)),"estLastYearPed"] <- 2021 # if our survival probability rule indicates an individual was alive at the end of the study, assign last year samples were collected for DNA extraction
PCC_dat[which(PCC_dat$estLastYear > 2021),"estLastYearPed"] <- 2021 #  # if our survival probability rule indicates an individual was likely died after the last year samples were collected for DNA extraction, last year samples were collected for DNA extraction
PCC_dat$longevity <- as.numeric(PCC_dat$estLastYear) - as.numeric(PCC_dat$estBirthYear)
PCC_dat$yearsContribute2Pedigree <- as.numeric(PCC_dat$estLastYearPed) - as.numeric(PCC_dat$estBirthYear)


# ------------------------------------------------------------------------------------------------------------

# visualize longevity 
# ------------------------------------------------------------------------------------------------------------

# technically, this should be referred to as the number of years an individual had the potential to contribute to the pedigree 
par(mfrow=c(1,2))
hist(ELF_dat$yearsContribute2Pedigree, main = "ELF", xlab = "years with potential to contribute to pedigree")
hist(PCC_dat$yearsContribute2Pedigree, main = "PCCI", xlab = "years with potential to contribute to pedigree")

# ------------------------------------------------------------------------------------------------------------

# match individual ID to DNA_IDs 
# ------------------------------------------------------------------------------------------------------------
ELF_dat$ID <- ELF_survival_data$ID
PCC_dat$ID <- PCC_survival_data$ID

PCC_survival_data[which(PCC_dat$yearsContribute2Pedigree > 16),]
PCC_dat[which(PCC_dat$yearsContribute2Pedigree > 16),]

# 11868056 --> large SVL, 63.4
# 27331044 --> large SVL, 63.3

load("../inbreeding_models/data_for_analyses_06252024.Robj") # data

head(data)
dim(data) # 1036   21

#data <- merge(x = data, y = ELF_dat, by.x = "id", by.y = "ID", all.x = TRUE, all.y = TRUE) # add ELF survival/longevity

PCC_dat <- cbind.data.frame(PCC_dat, PCC_survival_data$Fgrm, PCC_survival_data$Last.cap)
colnames(PCC_dat) <- c("estBirthYear", "estLastYear", "estLastYearPed", "longevity", "yearsContribute2Pedigree", "ID", "Fgrm", "Last.cap")

ELF_dat <- cbind.data.frame(ELF_dat, ELF_survival_data$Fgrm, ELF_survival_data$Last.cap)
colnames(ELF_dat) <- c("estBirthYear", "estLastYear", "estLastYearPed", "longevity", "yearsContribute2Pedigree", "ID", "Fgrm", "Last.cap")

temp_surv_dat <- rbind.data.frame(ELF_dat, PCC_dat)

data$Fgrm_rounded <- round(data$Fgrm, digits = 9)

# some error in merging dataframes due to rounding, fix here
temp_surv_dat[which(temp_surv_dat$ID == 27312577),"Fgrm"] <- -0.000092052
temp_surv_dat[which(temp_surv_dat$ID == 27076589),"Fgrm"] <- 0.000041382
temp_surv_dat[which(temp_surv_dat$ID == "ELF_854"),"Fgrm"] <- -0.000094367

# match based on Fgrm (easier than sorting out some PIT tags having 0s in front in some datasets and not others)
temp <- merge(x = data, y = temp_surv_dat, by.x = "Fgrm_rounded", by.y = "Fgrm", all.x = TRUE, all.y = TRUE)


head(temp)


plot(temp$offspring~temp$yearsContribute2Pedigree, pch = 19, col = alpha("black", alpha = 0.15), 
     xlab = "years contributing to pedigree", ylab = "number of offspring")

plot(temp$yearsContribute2Pedigree~temp$Fgrm, pch = 19, col = alpha("black", alpha = 0.15), 
     xlab = "Fgrm", ylab = "years contributing to pedigree")
# ------------------------------------------------------------------------------------------------------------

# How does birth year with Danielle's birth year rules change from my original estimates? 
# ------------------------------------------------------------------------------------------------------------

plot(round(temp$InferBirthYear), temp$estBirthYear, pch = 19, col = alpha("black", 0.1), 
     xlab = "original method", 
     ylab = "new method")
abline(a = 0, b = 1)

# ------------------------------------------------------------------------------------------------------------

# compare estimated last year with last year of capture
# ------------------------------------------------------------------------------------------------------------

temp[which(temp$Last.cap > temp$estLastYear),] # 0 



# Compile plots to share
# ------------------------------------------------------------------------------------------------------------
png("../figures/longevity_viz.png", width = 8, height = 5, units = "in", res = 400)
par(mfrow=c(2,3))
# age distributions 
hist(ELF_dat$yearsContribute2Pedigree, main = "ELF", xlab = "years with potential to contribute to pedigree")
hist(ELF_dat$longevity, main = "ELF", xlab = "longevity")
hist(PCC_dat$yearsContribute2Pedigree, main = "PCCI", xlab = "years with potential to contribute to pedigree")


# How does birth year with Danielle's birth year rules change from my original estimates? 
plot(round(temp$InferBirthYear), temp$estBirthYear, pch = 19, col = alpha("black", 0.1), 
     xlab = "birth year, original method", 
     ylab = "birth year, new method")
abline(a = 0, b = 1)

# how does longevity estimate correlate with Fgrm? 
plot(temp$yearsContribute2Pedigree~temp$Fgrm, pch = 19, col = alpha("black", alpha = 0.15), 
     xlab = "Fgrm", ylab = "years contributing to pedigree")

# how does longevity estimate corrlate with n. offspring? 
plot(temp$offspring~temp$yearsContribute2Pedigree, pch = 19, col = alpha("black", alpha = 0.15), 
     xlab = "years contributing to pedigree", ylab = "number of offspring")

dev.off()

# ------------------------------------------------------------------------------------------------------------

# Save data to use in pedigree reconstruction
# ------------------------------------------------------------------------------------------------------------

export_temp <- temp[,c("id", "estBirthYear", "estLastYear", "estLastYearPed", "BirthYear", "ID", "Last.cap")]
save(export_temp, file = "../pedigree_reconstruction/DRB_byEstimates.Robj")

plot(x = temp$BirthYear, y = temp$estBirthYear) # Danielle's method is more accurate! 
abline(a = 0, b = 1)

plot(x = temp$BirthYear, y = round(temp$InferBirthYear))
abline(a = 0, b = 1)


# what if we use different rounding rules?
# ------------------------------------------------------------------------------------------------------------
0.1
0.05
0.005



ELF_death_years_0.005 <- vector(length = nrow(ELF_survival_data))
ELF_death_years_0.05 <- vector(length = nrow(ELF_survival_data))
ELF_death_years_0.1 <- vector(length = nrow(ELF_survival_data))

# determine first year where survival probability is  0
# deal with cases where length = 0
for(i in 1:nrow(ELF_survival_data)){
  ELF_death_years_0.005[i] <- colnames(ELF_survival_data)[5:39][min(which(ELF_survival_data[i,5:39] <= 0.005))]
  ELF_death_years_0.05[i] <- colnames(ELF_survival_data)[5:39][min(which(ELF_survival_data[i,5:39] <= 0.05))]
  ELF_death_years_0.1[i] <- colnames(ELF_survival_data)[5:39][min(which(ELF_survival_data[i,5:39] <= 0.1))]
  
  # if(min(rounded_probs, na.rm = TRUE) > 0){
  #   ELF_death_years[i] <- NA
  # }else{
  #   ELF_death_years[i] <- as.numeric(colnames(ELF_survival_data)[5:39][head(which(rounded_probs == 0), n = 1)])
  # }
}

par(mfrow = c(1,3))
hist(as.numeric(ELF_death_years_0.005), main = "0.005")
hist(as.numeric(ELF_death_years_0.05), main = "0.05")
hist(as.numeric(ELF_death_years_0.1), main = "0.1")

ELF_longevity_0.005 <- as.numeric(ELF_death_years_0.005) - as.numeric(ELF_birth_years)
ELF_longevity_0.05 <- as.numeric(ELF_death_years_0.05) - as.numeric(ELF_birth_years)
ELF_longevity_0.1 <- as.numeric(ELF_death_years_0.1) - as.numeric(ELF_birth_years)

plot(ELF_death_years_0.005~ELF_birth_years)


hist(ELF_longevity_0.005,  main = "0.005")
hist(ELF_longevity_0.05,  main = "0.05")
hist(ELF_longevity_0.1,  main = "0.1")

sum(as.numeric(ELF_death_years_0.005) < ELF_survival_data$Last.cap, na.rm = T) / length(na.omit(ELF_death_years_0.005)) *100
sum(as.numeric(ELF_death_years_0.05) < ELF_survival_data$Last.cap, na.rm = T) / length(na.omit(ELF_death_years_0.05)) *100
sum(as.numeric(ELF_death_years_0.1) < ELF_survival_data$Last.cap, na.rm = T) / length(na.omit(ELF_death_years_0.1)) *100

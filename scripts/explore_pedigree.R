## File name: explore_pedigree.R
## Purpose: explore pedigrees for EMR project
## M. I. Clark, October 2022
## Last updated: 04/14/2024

## Loading ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

## load libraries
library(sequoia)
library(kinship2)
library(sp)
library(sf)
library(rgdal)
library(ggmap)
library(ggsn)
library(gganimate)
library("MetBrewer")
library(tidyverse) 
library(ggplot2)
library(tmap)
library(ape)
library(ncf)
library(labdsv)
library(FedData)
library(ggspatial)
library(gt)
library(broom)
library(raster)



## define color palette 
palette = "Archambault"
met.brewer(palette, n=6, type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = 6)


## date we want to analyze
# date <- 10182022
# date <- 11142022
# date <- "02102023"
# date = 040502023
# date = "05262023"
# date = "04132024"
date = "05072024"

## load data from make_prelim_pedigree.R
    # contain: 
            # 1. genotype information
            # 2. metadata
            # 3. duplicate individuals
            # 4. results from full pedigree run
            # 5. vcf filtering date

load(paste0("../pedigree_reconstruction/ELF_pedigree_results_", date, ".Robj"), verbose = T)
load(paste0("../pedigree_reconstruction/PCC_pedigree_results_", date, ".Robj"), verbose = T)

ELF_pedigree <- ELF_results[[4]]$Pedigree
PCC_pedigree <- PCC_results[[4]]$Pedigree

# metadata from pedigree reconstruction
ELF_meta <- ELF_results[[2]]
PCC_meta <- PCC_results[[2]]

# "extra" metadata for analyses
load(paste0("../pedigree_exploration/ELF_expanded_metaData.Robj"), verbose = T) # ELF_exp_metaData
load(paste0("../pedigree_exploration/PCC_expanded_metaData.Robj"), verbose = T) # PCC_exp_metaData

# coords
load(paste0("../pedigree_exploration/coords", "04132024", ".Robj"), verbose = TRUE)
PCC_coords <- coords[[1]]
ELF_coords <- coords[[2]]
# ------------------------------------------------------------------------------------------------------------

## Define custom functions ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
source("./pedigree_funcs.R") ## coordinate functions appear to be broken

# ------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------
## Correct duplicate names in pedigree  ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

ELF_pedigree$id[which(grepl("_P", ELF_pedigree$id))] <- paste0(unlist(strsplit(ELF_pedigree$id[which(grepl("_P", ELF_pedigree$id))], split = "_"))[1], "_", unlist(strsplit(ELF_pedigree$id[which(grepl("_P", ELF_pedigree$id))], split = "_"))[2])

PCC_pedigree$id[which(grepl("_P", PCC_pedigree$id))][1] <- paste0(unlist(strsplit(PCC_pedigree$id[which(grepl("_P", PCC_pedigree$id))], split = "_"))[1], "_", unlist(strsplit(PCC_pedigree$id[which(grepl("_P", PCC_pedigree$id))], split = "_"))[2])
PCC_pedigree$id[which(grepl("_P", PCC_pedigree$id))][1] <- paste0(unlist(strsplit(PCC_pedigree$id[which(grepl("_P", PCC_pedigree$id))], split = "_"))[1], "_", unlist(strsplit(PCC_pedigree$id[which(grepl("_P", PCC_pedigree$id))], split = "_"))[2])
PCC_pedigree$id[which(grepl("_P", PCC_pedigree$id))][1] <- paste0(unlist(strsplit(PCC_pedigree$id[which(grepl("_P", PCC_pedigree$id))], split = "_"))[1], "_", unlist(strsplit(PCC_pedigree$id[which(grepl("_P", PCC_pedigree$id))], split = "_"))[2])
# ------------------------------------------------------------------------------------------------------------

## Filter captive born individuals from ELF ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# make vector of individuals held for partuition 

# need to change: only want to filter offspring that were recorded because they were 
# born in captivity, not all offspring assigned to captive moms 

captive.inds <- list()
#captive.year <- vector()
count = 1
for(i in 1:length(ELF_exp_metaData)){
  if(any(ELF_exp_metaData[[i]]$held_captive == "Y")){
    print(i)
    captive.inds[[count]] <- subset(ELF_exp_metaData[[i]], held_captive == "Y", select = year)$year
    names(captive.inds)[count] <- names(ELF_exp_metaData)[[i]]
    print(subset(ELF_exp_metaData[[i]], held_captive == "Y", select = year)$year)
    count = count + 1
  }
}

# code above creates list where the list name is the name of the captive mom, and the contents of the list entry are the years she was held for partuition

captive.kids <- vector(mode = "list", length = length(captive.inds))
names(captive.kids) <- names(captive.inds)
# read in raw data to get mom info! 
ELF_rawData <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/ELF_LifeHistData.csv")[,1:45]

for(i in 1:length(captive.inds)){
  #print(ELF_rawData[grepl(strsplit(names(captive.inds)[i], "_")[[1]][2], ELF_rawData$Other.Notes),"Other.Notes"])
  #ELF_rawData[grepl(strsplit(names(captive.inds)[i], "_")[[1]][2], ELF_rawData$Other.Notes) & (grepl("mom", ELF_rawData$Other.Notes, ignore.case = TRUE) | grepl("mother", ELF_rawData$Other.Notes, ignore.case = TRUE)),"Other.Notes"]  
  captive.kids[[i]] <- ELF_rawData[grepl(strsplit(names(captive.inds)[i], "_")[[1]][2], ELF_rawData$Other.Notes) & (grepl("mom", ELF_rawData$Other.Notes, ignore.case = TRUE) | grepl("mother", ELF_rawData$Other.Notes, ignore.case = TRUE)),"ELF.ID"]
}

# filter didn't quite work for ELF_200
captive.kids[["ELF_200"]] <- captive.kids[["ELF_200"]][1:8]
# ELF_513 was held in lab for partuition in 2013, but passed 9-10 slugs, no neonates
captive.kids <- captive.kids[-which(names(captive.kids) == "ELF_513")]
captive.inds <- captive.inds[-which(names(captive.inds) == "ELF_513")]

# remove slug record
captive.kids[["ELF_161"]] <- captive.kids[["ELF_161"]][-6]

# find IDs of offspring assigned to captive mom and born in that year 

# make vector of their offspring determine which kids were only recorded as neonates in the lab and not recaptured
# captive.kids <- vector(mode = "list", length = length(captive.inds))
# names(captive.kids) <- names(captive.inds)

only_n_caps<- vector()
mult_caps<- vector() # individuals that have multiple captures (and thus were captured after birth, or were captured as adults)
need_help <- vector()
for(i in 1:length(captive.kids)){
  # captive.kids[[i]] <- find_children(rel.matrix, names(captive.inds)[i])
  #if(is.null(captive.kids[[i]])){
  #  print(paste(captive.inds[i], "had no offspring in pedigree"))
  #  break
  #}
  for (j in 1:length(captive.kids[[i]])){
    ind <- paste0("ELF_", captive.kids[[i]][j])
  #  if(!grepl("ELF", ind)){ # stop if this is a dummy individual
  #    break
  #  }else{
  #    # correct duplicate name
  #    if(grepl("_P", ind)){
  #      ind <- paste0(unlist(strsplit(ind, split = "_"))[1], "_", unlist(strsplit(ind, split = "_"))[2])
  #    }
      # index <- which(names(ELF_exp_metaData) == ind)
      if(length(ELF_exp_metaData[[ind]]$doy) == 1){ # only captured once
        # verify that age of capture is neonate and that the individual was a neonate in the year mom was held in lab
        if(any(ELF_exp_metaData[[ind]]$age.class == "N") & any(min(ELF_exp_metaData[[ind]]$year) == captive.inds[[i]])){
          print(paste0(ind, " was only captured as a neonate in year mom held in lab"))
          only_n_caps <- c(only_n_caps, ind)
        }else{
          if(any(ELF_exp_metaData[[ind]]$age.class == "A") | any(ELF_exp_metaData[[ind]]$age.class == "Y") | any(ELF_exp_metaData[[index]]$age.class == "J")){
            mult_caps <- c(mult_caps, ind)
          }else{
            need_help <- c(need_help, ind)
            print(paste0(ind, " was only captured once as a ", ELF_exp_metaData[[ind]]$age.class))
          }
        }
      }else{
        print(paste0(ind, " was captured more than once"))
        mult_caps <- c(mult_caps, ind)
      }
    }
  }

ELF_exp_metaData[["ELF_646"]] # born 2013
subset(ELF_pedigree, id == "ELF_646")
captive.inds[["ELF_394"]] # mom in lab 2012

# okay but some of only_n_caps might be parents or grandparents of individuals in the pedigree, indicating their survival past being a neonate

rel.matrix <- GetRelM(ELF_results[[4]]$Pedigree, GenBack = 1, Return = 'Matrix')
inferred.survivors <- vector()
for(i in 1:length(only_n_caps)){
  kids <- find_children(rel.matrix, only_n_caps[i])
  if(length(kids) > 0 ){
    print(paste0("inferred survival of ", only_n_caps[i]))
    inferred.survivors <- c(inferred.survivors, only_n_caps[i])
  }
}

only_n_caps <- subset(only_n_caps, !only_n_caps %in% inferred.survivors)

# length(only_n_caps)
# length(mult_caps)
# length(need_help)

# spot checking
# ELF_exp_metaData[[which(names(ELF_exp_metaData) == "ELF_760")]]

# only_n_caps is a vector of individuals who should not be included in reproductive output calculations because 
# they were only captured as neonates at birth

# filtered pedigree data object: 
ELF_pedigree_flt <- subset(ELF_pedigree, subset = !ELF_pedigree$id %in% only_n_caps)

save(ELF_pedigree_flt, file = paste0("../pedigree_exploration/EMR_ped_flt_", date, ".Robj"))
# ------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------
## PEDIGREE EXPLORATION
# ------------------------------------------------------------------------------------------------------------

### 1: How many individuals have parents?  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
percent_in_ped(PCC_pedigree, plot = TRUE)
# 10182022 : 
# "0.37 of individuals have no parents assigned"
# "0.52 of individuals have both parents assigned (genotyped or dummy)"
# "0.1 of individuals have one parent assigned (genotyped or dummy)"

# 11142022 : 
# "0.35 of individuals have no parents assigned"
# "0.53 of individuals have both parents assigned (genotyped or dummy)"
# "0.12 of individuals have one parent assigned (genotyped or dummy)"

# 02102023 : 
# "0.34 of individuals have no parents assigned"
# "0.54 of individuals have both parents assigned (genotyped or dummy)"
# "0.12 of individuals have one parent assigned (genotyped or dummy)"

# 05262023 :
# "0.32 of individuals have no parents assigned"
# "0.62 of individuals have both parents assigned (genotyped or dummy)"
# "0.06 of individuals have one parent assigned (genotyped or dummy)"


percent_in_ped(ELF_pedigree)
# 10182022 : 
# "0.08 of individuals have no parents assigned"
# "0.87 of individuals have both parents assigned (genotyped or dummy)"
# "0.05 of individuals have one parent assigned (genotyped or dummy)"

# 11142022 : 
# "0.08 of individuals have no parents assigned"
# "0.86 of individuals have both parents assigned (genotyped or dummy)"
# "0.05 of individuals have one parent assigned (genotyped or dummy)"

# 02102023 : 
# "0.07 of individuals have no parents assigned"
# "0.87 of individuals have both parents assigned (genotyped or dummy)"
# "0.05 of individuals have one parent assigned (genotyped or dummy)"

# 05262023 : 
# "0.05 of individuals have no parents assigned"
# "0.91 of individuals have both parents assigned (genotyped or dummy)"
# "0.03 of individuals have one parent assigned (genotyped or dummy)"


pdf(paste0("../pedigree_exploration/parents_assigned_", date, ".pdf"), width = 10, height = 6)
par(mfrow=c(1,2))
percent_in_ped(PCC_pedigree, plot = TRUE, legend = FALSE)
percent_in_ped(ELF_pedigree, plot = TRUE)
dev.off()


## talk figures

PCC_par <- percent_in_ped(PCC_pedigree, plot = FALSE, legend = FALSE) # length(geno.inds), parents.mia, both.parents, one.parent)
ELF_par <- percent_in_ped(ELF_pedigree, plot = FALSE, legend = FALSE)

met.brewer(palette, n=20, type="continuous")
MetBrewer::colorblind.friendly(palette)
more.colors <- met.brewer(palette, n = 20)

pdf(paste0("../pedigree_exploration/PCC_parents_assigned_legend_poster_", date, ".pdf"), width = 4, height = 6)
percents <- sapply(c(PCC_par[[2]], PCC_par[[3]], PCC_par[[4]]), function(x){x*100/PCC_par[[1]]})
barplot(as.matrix(percents), 
        xlim = c(0,4), 
        ylim = c(0,100), 
        col = more.colors[c(8,12,19)], 
        border = NA, 
        main = paste0("Barry County"), cex.axis = 1.5, cex.main = 1.5)
#legend("right", legend = rev(c("No parents assigned", "Both parent assigned", "One parents assigned")), 
#       pch = 15, col = rev(more.colors[c(8,12,19)]), 
#       bty="n", cex = 1.5)
dev.off()

pdf(paste0("../pedigree_exploration/ELF_parents_assigned_poster_", date, ".pdf"), width = 4, height = 6)
percents <- sapply(c(ELF_par[[2]], ELF_par[[3]], ELF_par[[4]]), function(x){x*100/ELF_par[[1]]})
barplot(as.matrix(percents), 
        xlim = c(0,4), 
        ylim = c(0,100), 
        col = more.colors[c(8,12,19)], 
        border = NA, 
        main = paste0("Cass County"), cex.axis = 1.5, cex.main = 1.5)
# legend("right", legend = c("No parents assigned", "One parent assigned", "Both parents assigned"), 
#        pch = 15, col = more.colors[c(8,12,19)], 
#        bty="n")
dev.off()


# ------------------------------------------------------------------------------------------------------------

### 2: Offspring per individual  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
PCC.off <- offspring_per_ind(PCC_pedigree)
range(PCC.off[[1]])
range(PCC.off[[2]])

ELF.off <- offspring_per_ind(ELF_pedigree_flt)
range(ELF.off[[1]])
range(ELF.off[[2]])

1-(length(unlist(PCC.off))/nrow(PCC_pedigree)) # 62% don't have kids

1-(length(unlist(ELF.off))/nrow(ELF_pedigree)) # 70% don't have kids

pdf(file = paste0("../pedigree_exploration/repro_success_report_", date, ".pdf"), height = 6, width = 6)
# offspring_per_ind(PCC_pedigree)
# offspring_per_ind(ELF_pedigree_flt)
par(mfrow=c(2,2))
hist(PCC.off[[1]], main = "PCCI dams", xlab = "number of offspring")
hist(PCC.off[[2]], main = "PCCI sires", xlab = "number of offspring")
hist(ELF.off[[1]], main = "ELF dams", xlab = "number of offspring")
hist(ELF.off[[2]], main = "ELF sires", xlab = "number of offspring")

dev.off()

# number of offpsring by birthyear
pdf(file = paste0("../pedigree_exploration/no_off_by_by_", date, ".pdf"), height = 6, width = 10)
par(mfrow=c(1,2))
#PCC
dat1 <- cbind.data.frame(PCC_results[[4]]$LifeHistSib$id, PCC_results[[4]]$LifeHistSib$BY.est)
rownames(dat1) <- PCC_results[[4]]$LifeHistSib$id

dat <- cbind.data.frame(names(unlist(PCC.off)), unlist(PCC.off))
colnames(dat) <- c("ID", "no.off")
dat <- dat[grepl("PCC", dat$ID),]
dat <- cbind.data.frame(dat, dat1[match(dat$ID, dat1$`PCC_results[[4]]$LifeHistSib$id`),2])
colnames(dat) <- c("ID", "no.off", "by.est")

plot(x = dat$by.est, y = dat$no.off, 
     pch = 19, 
     col = colors[1], 
     ylab = "number of offspring", 
     xlab = "estimated birth year", 
     main = "PCC")
# ELF
dat1 <- cbind.data.frame(ELF_results[[4]]$LifeHistSib$id, ELF_results[[4]]$LifeHistSib$BY.est)
rownames(dat1) <- ELF_results[[4]]$LifeHistSib$id

dat <- cbind.data.frame(names(unlist(ELF.off)), unlist(ELF.off))
colnames(dat) <- c("ID", "no.off")
dat <- dat[grepl("ELF", dat$ID),]
dat <- cbind.data.frame(dat, dat1[match(dat$ID, dat1$`ELF_results[[4]]$LifeHistSib$id`),2])
colnames(dat) <- c("ID", "no.off", "by.est")

plot(x = dat$by.est, y = dat$no.off, 
     pch = 19, 
     col = colors[2], 
     ylab = "number of offspring", 
     xlab = "estimated birth year", 
     main = "ELF")

dev.off()

## offspring combo
pdf(file = paste0("../pedigree_exploration/repro_success_site_", date, ".pdf"), height = 10, width = 10)
par(mfrow=c(2,2))
hist(PCC.off[[1]], 
     xlab = "number of offspring", 
     main = paste0("Barry County, dams"), col = colors[1], cex.lab = 2, cex.axis = 2, cex.main = 2)
# mtext("includes dummy individuals")
hist(PCC.off[[2]], 
     xlab = "number of offspring", 
     main = paste0("Barry County, sires"), col = colors[1], cex.lab = 2, cex.axis = 2, cex.main = 2)
# mtext("includes dummy individuals")

hist(ELF.off[[1]], 
     xlab = "number of offspring", 
     main = paste0("Cass County, dams"), col = colors[2], cex.lab = 2, cex.axis = 2, cex.main = 2)
# mtext("includes dummy individuals")
hist(ELF.off[[2]], 
     xlab = "number of offspring", 
     main = paste0("Cass County, sires"), col = colors[2], cex.lab = 2, cex.axis = 2, cex.main = 2)
# mtext("includes dummy individuals")
dev.off()




# ------------------------------------------------------------------------------------------------------------

### 3: Generations  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
PCC_gen <- getGenerations(Ped = PCC_pedigree)
range(PCC_gen)

ELF_gen <- getGenerations(Ped = ELF_pedigree)
range(ELF_gen)

pdf(file = paste0("../pedigree_exploration/gens_back_", date, ".pdf"), height = 6, width = 10)
par(mfrow=c(1,2))
hist(PCC_gen, 
     main = "PCC", 
     xlab = "generations", 
     col = colors[1])
hist(ELF_gen, 
     main = "ELF", 
     xlab = "generations", 
     col = colors[2])
dev.off()

# ------------------------------------------------------------------------------------------------------------

### 4: Who are these super successful reproducers?  --------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# think about sampling bias --> remove neonates who are never re-capped? 
    # how to get this information from the metadata? 

ELF_supers <- find_super_reproducers(ELF.off, "ELF", ELF_meta)
PCC_supers <- find_super_reproducers(PCC.off, "PCC", PCC_meta)

pdf(file = paste0("../pedigree_exploration/supers_v_by_max_", date, ".pdf"), height = 6, width = 10)
boxplot(c(ELF_supers$BY.max, PCC_supers$BY.max)~c(rep("ELF", nrow(ELF_supers)), rep("PCC", nrow(PCC_supers))), 
        ylab = "BY.max", 
        xlab = NULL, 
        main = "BY.max of most successful reproducers", 
        col = colors[c(1,2)])
dev.off()

pdf(file = paste0("../pedigree_exploration/supers_v_by_est_", date, ".pdf"), height = 6, width = 10)
boxplot(c(PCC_results[[4]]$LifeHistSib[match(x = PCC_supers$ID, table = PCC_results[[4]]$LifeHistSib$id),"BY.est"], 
          ELF_results[[4]]$LifeHistSib[match(x = ELF_supers$ID, table = ELF_results[[4]]$LifeHistSib$id),"BY.est"])~c(rep("PCC", nrow(PCC_supers)), rep("ELF", nrow(ELF_supers))), 
        ylab = "birth year estimate", 
        xlab = NULL, 
        main = "BY.max of most successful reproducers", 
        col = colors[c(1,2)])
dev.off()

pdf(file = paste0("../pedigree_exploration/supers_v_by_", date, ".pdf"), height = 6, width = 10)
par(mfrow=c(1,2))
boxplot(c(ELF_supers$BY.max, PCC_supers$BY.max)~c(rep("ELF", nrow(ELF_supers)), rep("PCC", nrow(PCC_supers))), 
        ylab = "BY.max", 
        xlab = NULL, 
        main = "BY.max of most successful reproducers", 
        col = colors[c(1,2)])

boxplot(c(PCC_results[[4]]$LifeHistSib[match(x = PCC_supers$ID, table = PCC_results[[4]]$LifeHistSib$id),"BY.est"], 
          ELF_results[[4]]$LifeHistSib[match(x = ELF_supers$ID, table = ELF_results[[4]]$LifeHistSib$id),"BY.est"])~c(rep("PCC", nrow(PCC_supers)), rep("ELF", nrow(ELF_supers))), 
        ylab = "birth year estimate", 
        xlab = NULL, 
        main = "BY.max of most successful reproducers", 
        col = colors[c(1,2)])
dev.off()

# Most super reproducers at ELF
## ELF_152
focal.ind = "ELF_152"
ELF_pedigree_flt[grepl(focal.ind, ELF_pedigree_flt$id),] # no parents assigned
unique(ELF_pedigree_flt[grepl(focal.ind, ELF_pedigree_flt$dam),]$sire) # 5 unique mates
ELF_meta[grepl(focal.ind, ELF_meta$ID),]
ELF_exp_metaData[[which(names(ELF_exp_metaData) == focal.ind)]]
  # captured 7 times, from 2011-2014


## ELF_135
focal.ind = "ELF_135"
ELF_pedigree_flt[grepl(focal.ind, ELF_pedigree_flt$id),] # no parents assigned
unique(ELF_pedigree_flt[grepl(focal.ind, ELF_pedigree_flt$dam),]$sire) # 6 unique mates
ELF_meta[grepl(focal.ind, ELF_meta$ID),] 
ELF_exp_metaData[[which(names(ELF_exp_metaData) == focal.ind)]]
  # captured 11 times from 2011-2015

## ELF_49
focal.ind = "ELF_49"
ELF_pedigree_flt[grepl(focal.ind, ELF_pedigree_flt$id),] # no parents assigned
unique(ELF_pedigree_flt[grepl(focal.ind, ELF_pedigree_flt$sire),]$dam) # at least 6 unique mates
ELF_meta[grepl(focal.ind, ELF_meta$ID),] 
ELF_exp_metaData[[which(names(ELF_exp_metaData) == focal.ind)]]
  # captured 2 times from 2010-2012
## ELF_150
focal.ind = "ELF_150"
ELF_pedigree_flt[grepl("ELF_150", ELF_pedigree_flt$id),] # assigned to inferred parents
unique(ELF_pedigree_flt[grepl("ELF_150", ELF_pedigree_flt$sire),]$dam) # 7 unique mates
ELF_meta[grepl(focal.ind, ELF_meta$ID),] 
ELF_exp_metaData[[which(names(ELF_exp_metaData) == focal.ind)]]
  # captured once in 2011

# ------------------------------------------------------------------------------------------------------------

### 5: What is the average number of mates for a given individual?  ------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

PCC_mates <- mate_no(PCC_pedigree)
ELF_mates <- mate_no(ELF_pedigree_flt)

pdf(file = paste0("../pedigree_exploration/number_of_mates_", date, ".pdf"), height = 6, width = 10)
par(mfrow=c(1,2))
hist(PCC_mates$no.mates, main = "PCC", col = colors[1], xlab = "No. of mates")
hist(ELF_mates$no.mates, main = "ELF", col = colors[2], xlab = "No. of mates")
dev.off()

mean(PCC_mates$no.mates)
mean(ELF_mates$no.mates)

mean(PCC_mates$no.offspring)
mean(ELF_mates$no.offspring)

pdf(file = paste0("../pedigree_exploration/number_of_mates_v_off_", date, ".pdf"), height = 6, width = 10)
par(mfrow=c(1,2))
plot(x = PCC_mates$no.mates, y = PCC_mates$no.offspring, 
     xlab = "Number of mates", 
     ylab = "Number of offspring", 
     main = "PCC", 
     pch = 19, 
     col = colors[1])
plot(x = ELF_mates$no.mates, y = ELF_mates$no.offspring, 
     xlab = "Number of mates", 
     ylab = "Number of offspring", 
     main = "ELF", 
     pch = 19, 
     col = colors[2])
dev.off()

# ------------------------------------------------------------------------------------------------------------

### 6: Where are super reproducers captured?  ------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# isolate super reproducers 
  
ELF_supers_coords_dams <- get_coords(ELF_exp_metaData, ELF_supers[which(ELF_supers$sex == 1),"ID"], 
                                     easting_col = "UTM.easting", 
                                     northing_col = "UTM.northing", 
                                     age_col = "age", 
                                     year_col = "year")
ELF_supers_coords_sires <- get_coords(ELF_exp_metaData, ELF_supers[which(ELF_supers$sex == 2),"ID"], 
                                      easting_col = "UTM.easting", 
                                      northing_col = "UTM.northing", 
                                      age_col = "age", 
                                      year_col = "year")

PCC_supers_coords_dams <- get_coords(PCC_exp_metaData, PCC_supers[which(PCC_supers$sex == 1),"ID"], 
                                     easting_col = "DD.easting", 
                                     northing_col = "DD.northing", 
                                     age_col = "age.class", 
                                     year_col = "year")
PCC_supers_coords_sires <- get_coords(PCC_exp_metaData, PCC_supers[which(PCC_supers$sex == 2),"ID"], 
                                     easting_col = "DD.easting", 
                                     northing_col = "DD.northing", 
                                     age_col = "age.class", 
                                     year_col = "year")

# plot
plot_snakes <- function(coords, easting = 2, northing = 3, 
                        convert.coords = TRUE, 
                        box = c(left = -110.0145, bottom = 41.942, right = -109.9960, top = 41.962), 
                        zoom = 12,
                        map.type = "terrain", 
                        plot.title = NULL, 
                        scale = TRUE){
  # testing vars
  # coords = ELF_coords
  # plot.title = "PCC super damns"
  # convert.coords = FALSE
  # box = c(left = -85.279, bottom = 42.532, right = -85.305, top = 42.5445)
  
  
  # clean coords
  coords_clean <- subset(coords, subset = !coords[,easting] == ".")
  coords_clean[,2:3] <- sapply(coords_clean[,c(easting, northing)], as.numeric)
  coords_clean <-  subset(coords_clean, subset = !is.na(coords_clean[,easting]))
  
  #convert to lat long
  if(convert.coords == TRUE){
    sputm <- SpatialPoints(coords_clean[c(easting, northing)], proj4string=CRS("+proj=utm +zone=12 +datum=WGS84")) 
    spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
  }else{
    spgeo<- SpatialPoints(coords_clean[c(easting, northing)], proj4string=CRS("+proj=longlat +datum=WGS84")) 
  }

  # plot(spgeo)
  dat <- SpatialPointsDataFrame(coords = spgeo, data = coords_clean[,c(1,4,5)])
  dat.df<-data.frame(dat)
  colnames(dat.df) <- c("ID", "age", "year", "long", "lat", "optional")
 
  map <- get_stamenmap(bbox = box,  
                          zoom = zoom, maptype = map.type)
  
  p <-ggmap(map) + geom_point(aes(x = long, y = lat, color = ID), 
                 data = dat.df, size = 2.5) + 
    geom_path(aes(x = long, y = lat, color = ID), 
              data = dat.df) +
    theme(legend.position = "none") + 
    labs(title=plot.title) + 
    scale_color_manual(values=met.brewer(palette, n = nrow(coords)))
    #geom_sf()
    #north(anndat.dat) +
  if(scale == TRUE){
    p +     scalebar(dat.df, location = "topleft", dist = 200,
                     dist_unit = "m", transform = TRUE, st.size = 3.25, 
                     nudge_y = -0.0002)
  }else{
    p
  }
}

pdf(file = "../pedigree_exploration/ELF_supers_map.pdf", width = 6, height = 6)
plot_snakes(ELF_supers_coords_sires, map.type = "terrain", zoom = 14, plot.title = "ELF super sires")
plot_snakes(ELF_supers_coords_dams, map.type = "terrain", zoom = 14, plot.title = "ELF super dams")
dev.off()

pdf(file = "../pedigree_exploration/PCC_supers_map.pdf", width = 6, height = 6)
plot_snakes(PCC_supers_coords_dams, plot.title = "PCC super damns", convert.coords = FALSE, 
            box = c(left = -85.300, bottom = 42.53, right = -85.285, top = 42.545), 
            zoom = 13, scale = FALSE)
plot_snakes(PCC_supers_coords_sires, plot.title = "PCC super damns", convert.coords = FALSE, 
            box = c(left = -85.300, bottom = 42.53, right = -85.285, top = 42.545), 
            zoom = 13, scale = FALSE)
dev.off()

#
# ------------------------------------------------------------------------------------------------------------

### 7: Pedigrees as family groups ------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

make_fam_peds <- function(pedigree){
  fams <- FindFamilies(pedigree)
  site <- unlist(strsplit(fams$id, split = "_"))[1]
  fam.IDs <- unique(fams$FID)
  # make trees for families
  for(fid in fam.IDs){
    ped <- fams[which(fams$FID == fid),]
    if(nrow(ped) > 1){
      PedP <- sequoia::PedPolish(ped, DropNonSNPd=FALSE,
                                 FillParents = TRUE)
      PedP$Sex <- with(PedP, ifelse(id %in% dam, "female",  "male"))
      # default to 'male' to avoid warning: "More than 25% of the gender values are
      #  'unknown'"
      
      Ped.fix <- with(PedP, kinship2::fixParents(id=id, dadid=sire, momid=dam,
                                                 sex=Sex))
      Ped.k <- with(Ped.fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))
      
      genotyped <- grep(site, Ped.k$id)
      not_geno <- grep(site, Ped.k$id, invert=TRUE)
      affected <- vector(length=length(Ped.k$id))
      affected[genotyped] <- 1
      affected[not_geno] <- 0
      plot(Ped.k, id = rep("", dim(PedP)[1]), symbolsize=1, 
           density = c(-1, 35, 65, 20), 
           affected = affected)
    }
  }
  # To get the output pedigree into kinship2 compatible format:
}  

pdf(file = paste0("../pedigree_exploration/ELF_fams.pdf"), height = 10, width = 10)
make_fam_peds(ELF_pedigree_flt)
dev.off()
  
pdf(file = paste0("../pedigree_exploration/PCC_fams.pdf"), height = 10, width = 10)
make_fam_peds(PCC_pedigree)
dev.off()

# ------------------------------------------------------------------------------------------------------------

### 8: Reproductive skew -----------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

pdf(file = paste0("../pedigree_exploration/estimated_by.pdf"), height = 6, width = 10)
par(mfrow=c(1,2))
hist(PCC_results[[4]]$LifeHistSib$BY.est, 
     main = "PCC", 
     xlab = "Estimated birth year", 
     col = colors[1])
hist(ELF_results[[4]]$LifeHistSib$BY.est, 
     main = "ELF", 
     xlab = "Estimated birth year", 
     col = colors[2])
dev.off()

pdf(file = paste0("../pedigree_exploration/estimated_by_cohorts.pdf"), height = 6, width = 10)
par(mfrow=c(1,2))
find_old_young_cohorts(PCC_results, percent.older = 0.5, percent.younger = 0.1)
find_old_young_cohorts(ELF_results, percent.older = 0.5, percent.younger = 0.1)
dev.off()

PCC_cohorts <- find_old_young_cohorts(PCC_results, percent.older = 0.25, percent.younger = 0.25)
PCC_lines <- find_descendants(PCC_results, PCC_cohorts[[1]])
total <- matrix(data = NA, nrow = length(PCC_lines), ncol = 2)
for(i in 1:length(PCC_lines)){
  count = 0
  if(sum(PCC_lines[[i]] %in% PCC_cohorts[[2]]) > 0) {
    count = count + 1
  }
  total[i,] <- cbind(i, count)
}
paste0((sum(total[,2])/ length(PCC_cohorts[[2]]) * 100), " percent of young cohort are direct descendants of old cohort")

ELF_cohorts <- find_old_young_cohorts(ELF_results, percent.older = 0.25, percent.younger = 0.1)
ELF_lines <- find_descendants(ELF_results, ELF_cohorts[[1]])
total <- matrix(data = NA, nrow = length(ELF_lines), ncol = 2)
for(i in 1:length(ELF_lines)){
  count = 0
  if(sum(ELF_lines[[i]] %in% ELF_cohorts[[2]]) > 0){
    count = count + 1
  }
  total[i,] <- cbind(i, count)
}
paste0((sum(total[,2])/ length(ELF_cohorts[[2]]) * 100), " percent of young cohort are direct descendants of old cohort")

# how many descendants does each "old" individual have?
pdf(file = paste0("../pedigree_exploration/total_desc_of_old_cohort.pdf"), height = 6, width = 10)
par(mfrow = c(1,2))
hist(unlist(lapply(PCC_lines, length)), main = "Oldest 25% of PCC", xlab = "Number of descendants", col = colors[1])
hist(unlist(lapply(ELF_lines, length)), main = "Oldest 25% of ELF", xlab = "Number of descendants", col = colors[2])
dev.off()
# ------------------------------------------------------------------------------------------------------------

### 9: Lines over time -----------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# get bounding box
coords_clean <- subset(ELF_coords, subset = !ELF_coords[,"easting"] == ".")
coords_clean[,2:3] <- sapply(coords_clean[,c("easting", "northing")], as.numeric)
coords_clean <-  subset(coords_clean, subset = !is.na(coords_clean[,"easting"]))

sputm <- SpatialPoints(coords_clean[c("easting", "northing")], proj4string=CRS("+proj=utm +zone=16 +datum=WGS84")) 
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
spgeo@bbox


register_google("") # add API key

pdf(file = "../pedigree_exploration/ELF_152.pdf", width = 12, height = 6)
plot_ELF_lines("ELF_152")
dev.off()

pdf(file = "../pedigree_exploration/ELF_135.pdf", width = 12, height = 6)
plot_ELF_lines("ELF_135")
dev.off()

pdf(file = "../pedigree_exploration/ELF_49.pdf", width = 12, height = 6)
plot_ELF_lines("ELF_49")
dev.off()


pdf(file = "../pedigree_exploration/ELF_150.pdf", width = 12, height = 6)
plot_ELF_lines("ELF_150")
dev.off()

ELF_pedigree[which(ELF_pedigree$sire == "ELF_49"),]

#### old code
grepl(ELF_lines[[ind]], ELF_coords$ID)
ELF_coords
coords <- get_coords(ELF_exp_metaData, c(ind, ELF_lines[[ind]]), easting_col = "UTM.easting", 
                       northing_col = "UTM.northing", age_col = "age", year_col = "year")

plot_snakes(coords, zoom = 14, scale = FALSE)

coords_clean <- coords[-which(coords$UTM.easting == "."),]

sapply(coords_clean, class)
coords_clean[,2:3] <- sapply(coords_clean[,2:3], as.numeric)

#convert to lat long
sputm <- SpatialPoints(coords_clean[2:3], proj4string=CRS("+proj=utm +zone=12 +datum=WGS84")) 
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))

plot(spgeo)

anndat <- SpatialPointsDataFrame(coords = spgeo, data = coords_clean[,c(1,4,5)])
anndat.dat<-data.frame(anndat)
colnames(anndat.dat) <- c("ID", "age", "year", "long", "lat", "optional")

bc_map <- get_stamenmap(bbox = c(left = -110.013, bottom = 41.949, right = -109.9993, top = 41.959),  
                        zoom = 14, maptype = "terrain")

years <- coords_clean$year
plot_cols <- met.brewer(palette, n = length(unique(coords_clean$ID)))

#pdf(file = "../pedigree_exploration/ELF_54_line/ELF_54_line_map.pdf", height = 6, width = 6)
#for(year in range(years)[1]:range(years)[2]){
#  data <- anndat.dat[which(anndat.dat$year == year),]
#  data <- data[!duplicated(data$ID),]
#  bc_map <- get_stamenmap(bbox = c(left = -110.013, bottom = 41.949, right = -109.9993, top = 41.959),  
#                          zoom = 14, maptype = "terrain")
#  p <-ggmap(bc_map)
#  map <- p + geom_point(aes(x = long, y = lat, color = ID), 
#                 data = data, size = 2.5, 
#                 show.legend = FALSE) + 
#    theme(legend.position="bottom") + 
#    labs(title=year) + 
#    scale_color_manual(values=plot_cols)
#  print(map)
#}
#dev.off()

data <- anndat.dat
bc_map <- get_stamenmap(bbox = c(left = -110.013, bottom = 41.948, right = -109.9993, top = 41.959),  
                        zoom = 14, maptype = "terrain")
p <-ggmap(bc_map)
map <- p + geom_point(aes(x = long, y = lat, color = ID), 
                      data = data, size = 2.5, 
                      show.legend = FALSE) + 
  theme(legend.position="bottom") + 
  scale_color_manual(values=plot_cols)
map.animation = map +
  transition_time(as.integer(years)) +
  labs(subtitle = "Year: {frame_time}") +
  shadow_wake(wake_length = 0.01)

animate(map.animation, height = 400, width = 400, fps = 20, duration = 30,
          end_pause = 10, res = 100)

anim_save(paste0("../pedigree_exploration/lines/",ind,"_line.gif"))


# ------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------

### 11: Isolate individuals with known ages -----------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# ELF_pedigree_flt and PCC_pedigree
#ELF_meta
#PCC_meta

ELF_by <- as.numeric(ELF_meta$BirthYear)
PCC_by <- as.numeric(PCC_meta$BirthYear)

ELF_inds_w_by <- ELF_meta$ID[which(!ELF_by == -9)]
PCC_inds_w_by <- PCC_meta$ID[which(!PCC_by == -9)]

ELF_inds_w_by <- cbind.data.frame(ELF_inds_w_by, ELF_meta$BirthYear[which(!ELF_by == -9)])
PCC_inds_w_by <- cbind.data.frame(PCC_inds_w_by, PCC_meta$BirthYear[which(!PCC_by == -9)])

colnames(ELF_inds_w_by) <- c("id","by")
colnames(PCC_inds_w_by) <- c("id", "by")

ELF_inds_w_by$by <- as.numeric(ELF_inds_w_by$by)
PCC_inds_w_by$by <- as.numeric(PCC_inds_w_by$by)
hist(ELF_inds_w_by$by)
hist(PCC_inds_w_by$by )
# ------------------------------------------------------------------------------------------------------------

## 12: Does parental relatedness predict number of offspring? ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
#PCC
# get pairwise matrix of relatedness for all individuals
PCC_rel_mat <- CalcRped(PCC_pedigree, OUT = "M")

# identify list of unique parent combos from pedigree
parents <- na.omit(PCC_pedigree[,c("dam", "sire")])
unique_rents <- unique(parents)

PCC_n_kids <- vector(length = length(unique_rents))
PCC_pair_rel <- vector(length = length(unique_rents))
for(p in 1:nrow(unique_rents)){
  # get number of offspring for each combo 
  PCC_n_kids[p] <- nrow(subset(parents, dam %in% unique_rents[p,] & sire %in% unique_rents[p,]))
  
  # get relatedness for each combo 
  PCC_pair_rel[p] <- PCC_rel_mat[which(rownames(PCC_rel_mat) == unique_rents[p,2]),which(colnames(PCC_rel_mat) == unique_rents[p,1])]
}
hist(PCC_n_kids)
hist(PCC_pair_rel)

plot(x = PCC_pair_rel, y = PCC_n_kids)


# get pairwise matrix of relatedness for all individuals
ELF_rel_mat <- CalcRped(ELF_pedigree, OUT = "M")

# identify list of unique parent combos from pedigree
ELF_parents <- na.omit(ELF_pedigree[,c("dam", "sire")])
ELF_unique_rents <- unique(ELF_parents)

ELF_n_kids <- vector(length = length(ELF_unique_rents))
ELF_pair_rel <- vector(length = length(ELF_unique_rents))
for(p in 1:nrow(ELF_unique_rents)){
  # get number of offspring for each combo 
  ELF_n_kids[p] <- nrow(subset(ELF_parents, dam %in% ELF_unique_rents[p,] & sire %in% ELF_unique_rents[p,]))
  
  # get relatedness for each combo 
  ELF_pair_rel[p] <- ELF_rel_mat[which(rownames(ELF_rel_mat) == ELF_unique_rents[p,2]),which(colnames(ELF_rel_mat) == ELF_unique_rents[p,1])]
}

hist(ELF_n_kids)
hist(ELF_pair_rel)
plot(x = ELF_pair_rel, y = ELF_n_kids)

table(PCC_n_kids, PCC_pair_rel)
# x y sample.size
# 0 1 30
# 0 2 27 
PCC_table <- table(PCC_n_kids, PCC_pair_rel)

pdf(file = paste0("../pedigree_exploration/offspring_vs_ped_rel_", date, ".pdf"), height = 4, width = 8)
par(mfrow = c(1,2))
plot(x = jitter(PCC_pair_rel, factor = 0.25), y = jitter(PCC_n_kids, factor = 0.25), xlab = "parent pedigree relatedness", ylab = "number of offspring", col = alpha(colors[1], alpha = 0.25), pch = 19)
plot(x = jitter(ELF_pair_rel), y = jitter(ELF_n_kids),  xlab = "parent pedigree relatedness", ylab = "number of offspring", col = alpha(colors[2], alpha = 0.25), pch = 19)
dev.off()

pdf(file = paste0("../pedigree_exploration/offspring_vs_ped_rel_poster_", date, ".pdf"), height = 5, width = 10)
par(mfrow = c(1,2))
cex_val = 2
plot(x = jitter(PCC_pair_rel, factor = 0.25), y = jitter(PCC_n_kids, factor = 0.25), xlab = "parent pedigree relatedness", ylab = "number of offspring", col = alpha(colors[1], alpha = 0.25), pch = 19, 
     cex = cex_val, cex.lab = cex_val, 
     cex.axis = cex_val, cex.main =cex_val)
plot(x = jitter(ELF_pair_rel), y = jitter(ELF_n_kids),  xlab = "parent pedigree relatedness", ylab = "number of offspring", col = alpha(colors[2], alpha = 0.25), pch = 19, cex = cex_val, cex.lab = cex_val, 
     cex.axis = cex_val, cex.main =cex_val)
dev.off()


pdf(file = paste0("../pedigree_exploration/pair_ped_rel_hist_", date, ".pdf"), height = 6, width = 12)
par(mfrow = c(1,2))
hist(PCC_pair_rel, col = colors[1], xlab = "parent pedigree relatedness", main = NULL)
hist(ELF_pair_rel, col = colors[2], xlab = "parent pedigree relatedness", main = NULL)
dev.off()

# who are the pairs with high pedigree relatedness from ELF? 

related_pairs <- ELF_unique_rents[which(ELF_pair_rel > 0.0),]

ELF_relationships <- GetRelM(Pedigree = ELF_pedigree, GenBack = 2, patmat = TRUE, Return = "Array")

ELF_relationship_cats <- matrix(data = NA, nrow = nrow(related_pairs), ncol=dim(ELF_relationships)[3])
for(a in 1:dim(ELF_relationships)[3]){
  for(i in 1:nrow(related_pairs)){
    ELF_relationship_cats[i,a] <- ELF_relationships[which(rownames(ELF_relationships) == related_pairs[i,1]),which(colnames(ELF_relationships) == related_pairs[i,2]),a]
  }
}
colnames(ELF_relationship_cats) <- dimnames(ELF_relationships)$Rel

# simplified
ELF_relationships_simple <- GetRelM(Pedigree = ELF_pedigree, GenBack = 2, patmat = FALSE, Return = "Matrix")

ELF_relationship_cats_simple <- vector(length = nrow(related_pairs))
  for(i in 1:nrow(related_pairs)){
    ELF_relationship_cats_simple[i] <- ELF_relationships_simple[which(rownames(ELF_relationships_simple) == related_pairs[i,1]),which(colnames(ELF_relationships_simple) == related_pairs[i,2])]
  }

related_pairs<- cbind.data.frame(related_pairs, ELF_relationship_cats_simple)

#related_pairs$rel.type <- c("GO/GP", "GO/GP", "GO/GP", "HS", "P/O", "P/O", "P/O", "HS", "U", "FC1", "HA/HN", "P/O", "FN/FA", "FN/FA", "HA/HN", "HA/HN", "HS", "FN/FA", "HA/HN", "HS", "GO/GP", "P/O", "U", "HS", "P/O", "P/O")
#related_pairs$rel.type <- c("GO/GP", "HS", "FC1", "HA/HN", "P/O", "FN/FA", "FN/FA", "P/O", "HA/HN", "HA/HN", "HS", "FN/FA", "HS", "FS", "P/O", "U", "HS")
related_pairs$rel.type <- c("P/O", "U", "HA/HN", "P/O", "GO/GP", "P/O", "U", "U", "GO/GP", "P/O", "U", "FN/FA", "U", "HA/HN", "U", "HS", "FS", "P/O", "U", "HS", "U")


pdf(file = paste0("../pedigree_exploration/related_pair_categories", date, ".pdf"), height = 6, width = 10)
barplot(table(ELF_relationship_cats_simple), col = colors[2])
dev.off()

pdf(file = paste0("../pedigree_exploration/related_pair_categories_types_", date, ".pdf"), height = 6, width = 10)
barplot(table(related_pairs$rel.type), col = colors[2])
dev.off()

plot(x = jitter(ELF_pair_rel), y = jitter(ELF_n_kids),  xlab = "parent pedigree relatedness", ylab = "number of offspring", col = alpha(colors[2], alpha = 0.25), pch = 19)

ELF_pair_rel_kids <- list(ELF_pair_rel, ELF_n_kids)

save(ELF_pair_rel_kids, file = "../pedigree_exploration/ELF_pair_rel_kids.Robj")

library(bbmle)
detFunc <- function(x,a,b){
  return(a*exp(-b*x))
}

# define an inverse link function
#	in this case, the identity
invLink <- function(z){
  return(exp(z))
  }


mydata <- data.frame("x"= ELF_pair_rel,"y"= ELF_n_kids) 
mod <- mle2(y ~ dpois(lambda=invLink(
  detFunc(x,a,b) )
), data=mydata,
start=list("a"=0,"b"=0))

summary(mod)
lines(mydata$x,detFunc(mydata$x,a=mod@coef[1],b=mod@coef[2]), col = "red")

plot(ELF_n_kids ~ ELF_pair_rel)
library(glmmTMB)
data <- cbind.data.frame(ELF_n_kids, ELF_pair_rel)
colnames(data) <- c("kids", "pair_rel")
mod1 <- glmmTMB(ELF_n_kids ~ ELF_pair_rel,
                ziformula = ~.,
                data = data,
                family = poisson)
testDispersion(mod1)
plot(simulateResiduals(mod1))

# ------------------------------------------------------------------------------------------------------------

## 13: Proportion of dummy individuals? ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
sum(grepl("PCC_", PCC_pedigree$id)) / length(PCC_pedigree$id)

sum(grepl("ELF_", ELF_pedigree$id)) / length(ELF_pedigree$id)

# ------------------------------------------------------------------------------------------------------------

## 14: Notable individuals -----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
#47/500 largest female
# not in pedigree

# 811 smallest female
# assigned ELF_58 as dam

# captured the most years: 65
# no parents assigned
# assigned 16 offspring

ind <- "ELF_65"
colors = met.brewer("Archambault", n = 6)
descendants <- unlist(GetDescendants(ind, ELF_pedigree_flt)[c(1,2,3)])
coords <- ELF_coords[ELF_coords$ID %in% descendants,]
coords[which(coords$ID %in% GetDescendants(ind, ELF_pedigree_flt)$offspring),"status"] <- "offspring"
coords[which(coords$ID %in% GetDescendants(ind, ELF_pedigree_flt)$grandoffspring),"status"] <- "grandoffspring"
coords[which(coords$ID %in% GetDescendants(ind, ELF_pedigree_flt)$id),"status"] <- "self"

# clean coords
coords_clean <- subset(coords, subset = !coords[,"easting"] == ".")
coords_clean[,2:3] <- sapply(coords_clean[,c("easting", "northing")], as.numeric)
coords_clean <-  subset(coords_clean, subset = !is.na(coords_clean[,"easting"]))


#convert to lat long
sputm <- SpatialPoints(coords_clean[c("easting", "northing")], proj4string=CRS("+proj=utm +zone=16 +datum=WGS84")) 
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))

dat <- SpatialPointsDataFrame(coords = spgeo, data = coords_clean[,c(1,4,5,6)])
dat.df<-data.frame(dat)
colnames(dat.df) <- c("ID", "age", "year", "status", "long", "lat")

box = c(left = -86.01221, bottom = 41.94555, right = -85.98966, top = 41.95972)
zoom = 13
map.type = "satellite"
#map <- get_stamenmap(bbox = box, zoom = zoom, maptype = map.type)
map <- get_map(location = box, maptype = "satellite", source = "google")

self <- dat.df[which(dat.df$status == "self"),]

self.map <- ggmap(map) + geom_point(aes(x = long, y = lat, color = status), 
                             data = self, size = 2.5) + 
  geom_path(aes(x = long, y = lat, color = status), 
            data = self) +
  theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title=NULL) + ylab("") + xlab("") +
  scale_color_manual(values=colors[c(6)])

off <- dat.df[which(dat.df$status == "offspring"),]

off.map <- ggmap(map) + geom_point(aes(x = long, y = lat, color = status), 
                                    data = off, size = 2.5) + 
  geom_path(aes(x = long, y = lat, color = status), 
            data = off) +
  theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title=NULL) + ylab("") + xlab("") +
  scale_color_manual(values=colors[c(5)])

grand <- dat.df[which(dat.df$status == "grandoffspring"),]

grand.map <- ggmap(map) + geom_point(aes(x = long, y = lat, color = status), 
                                   data = grand, size = 2.5) + 
  geom_path(aes(x = long, y = lat, color = status), 
            data = grand) +
  theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title=NULL) + ylab("") + xlab("") +
  scale_color_manual(values=colors[c(4)])

pdf(file = "../pedigree_exploration/ELF_65_self.pdf", width = 12, height = 6)
self.map
dev.off()

pdf(file = "../pedigree_exploration/ELF_65_offspring.pdf", width = 12, height = 6)
off.map
dev.off()

pdf(file = "../pedigree_exploration/ELF_65_grandoffspring.pdf", width = 12, height = 6)
grand.map
dev.off()


pdf(file = "../pedigree_exploration/ELF_65.pdf", width = 12, height = 6)
plot_ELF_lines("ELF_65")
dev.off()

# 12 yo: 271
# assigned parents! 
# no assigned offspring


# ------------------------------------------------------------------------------------------------------------

## 15: SVL vs reproductive output -----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# get reproductive output of only real individuals
ELF.dams <- ELF.off[[1]][grepl("ELF", names(ELF.off[[1]]))]
ELF.sires <- ELF.off[[2]][grepl("ELF", names(ELF.off[[2]]))]

PCC.dams <- PCC.off[[1]][grepl("PCC", names(PCC.off[[1]]))]
PCC.sires <- PCC.off[[2]][grepl("PCC", names(PCC.off[[2]]))]

# get SVL of each individual

extract_SVL <- function(inds, exp_metaData){
  # inds <- names(ELF.dams)
  # data_type <- "SVL"
  # exp_metaData <- ELF_exp_metaData
  avg_SVL <- vector(length = length(inds))
  for(i in 1:length(inds)){
    snake_data <- exp_metaData[[inds[i]]]
    avg_SVL[i] <- mean(suppressWarnings(as.numeric(snake_data[which(snake_data$age.class == "A"),]$SVL)), na.rm = TRUE)
  }
  return(avg_SVL)
}

ELF.dams.SVL <- extract_SVL(names(ELF.dams), ELF_exp_metaData)
ELF.sires.SVL <- extract_SVL(names(ELF.sires), ELF_exp_metaData)

PCC.dams.SVL <- extract_SVL(names(PCC.dams), PCC_exp_metaData)
PCC.sires.SVL <- extract_SVL(names(PCC.sires), PCC_exp_metaData)

# plot ELF
par(mfrow=c(2,2))
plot(as.vector(ELF.dams)~ELF.dams.SVL, main = "Cass Dams", xlab = "average adult SVL", ylab = "number of offspring", pch = 19, col = colors[4])
abline(lm(as.vector(ELF.dams)~ELF.dams.SVL), lty = 2, lwd = 2)
summary(lm(as.vector(ELF.dams)~ELF.dams.SVL))

plot(as.vector(ELF.sires)~ELF.sires.SVL, main = "Cass Sires", xlab = "average adult SVL", ylab = "number of offspring", pch = 19, col = colors[5])
abline(lm(as.vector(ELF.sires)~ELF.sires.SVL), lty = 2, lwd = 2 )
summary(lm(as.vector(ELF.sires)~ELF.sires.SVL))

# plot PCC
plot(as.vector(PCC.dams)~PCC.dams.SVL, main = "Barry Dams", xlab = "average adult SVL", ylab = "number of offspring", pch = 19, col = colors[4])
abline(lm(as.vector(PCC.dams)~PCC.dams.SVL), lty = 2, lwd = 1)
summary(lm(as.vector(PCC.dams)~PCC.dams.SVL))

plot(as.vector(PCC.sires)~PCC.sires.SVL, main = "Barry Sires", xlab = "average adult SVL", ylab = "number of offspring", pch = 19, col = colors[5])
abline(lm(as.vector(PCC.sires)~PCC.sires.SVL), lty = 2, lwd = 1)
summary(lm(as.vector(PCC.sires)~PCC.sires.SVL))


## 15.5: SVL vs inbreeding -----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# load inbreeding data 
load(file = "../vcf_filtering/pop_gen_stats.Robj", verbose = T)
pop_gen_stats$id <- str_remove(pop_gen_stats$id, "_P[0-9]")

dim(pop_gen_stats)
pop_gen_stats[match(names(ELF.dams), pop_gen_stats$id),]

ELF.dams.Fgrm <- pop_gen_stats[match(names(ELF.dams), pop_gen_stats$id),"Fgrm"]
ELF.sires.Fgrm <- pop_gen_stats[match(names(ELF.sires), pop_gen_stats$id),"Fgrm"]

PCC.dams.Fgrm <- pop_gen_stats[match(names(PCC.dams), pop_gen_stats$id),"Fgrm"]
PCC.sires.Fgrm <- pop_gen_stats[match(names(PCC.sires), pop_gen_stats$id),"Fgrm"]

# plot ELF
par(mfrow=c(2,2))
plot(ELF.dams.Fgrm~ELF.dams.SVL, main = "Cass Dams", xlab = "average adult SVL", ylab = "Fgrm", pch = 19, col = colors[4])
abline(lm(ELF.dams.Fgrm~ELF.dams.SVL), lty = 2, lwd = 2)
summary(lm(ELF.dams.Fgrm~ELF.dams.SVL))

# plot(ELF.sires.Fgrm~ELF.sires.SVL, main = "Cass Sires", xlab = "average adult SVL", ylab = "Fgrm", pch = 19, col = colors[5])
# abline(lm(ELF.sires.Fgrm~ELF.sires.SVL), lty = 2, lwd = 2 )
# summary(lm(ELF.sires.Fgrm~ELF.sires.SVL))
# get rid of outlier
plot(ELF.sires.Fgrm[-20]~ELF.sires.SVL[-20], main = "Cass Sires", xlab = "average adult SVL", ylab = "Fgrm", pch = 19, col = colors[5])
abline(lm(ELF.sires.Fgrm[-20]~ELF.sires.SVL[-20]), lty = 2, lwd = 2 )
summary(lm(ELF.sires.Fgrm[-20]~ELF.sires.SVL[-20]))

# plot PCC
plot(PCC.dams.Fgrm~PCC.dams.SVL, main = "Cass Dams", xlab = "average adult SVL", ylab = "Fgrm", pch = 19, col = colors[4])
abline(lm(PCC.dams.Fgrm~PCC.dams.SVL), lty = 2, lwd = 2)
summary(lm(PCC.dams.Fgrm~PCC.dams.SVL))

plot(PCC.sires.Fgrm~PCC.sires.SVL, main = "Cass Sires", xlab = "average adult SVL", ylab = "Fgrm", pch = 19, col = colors[5])
abline(lm(PCC.sires.Fgrm~PCC.sires.SVL), lty = 2, lwd = 2 )
summary(lm(PCC.sires.Fgrm~PCC.sires.SVL))

# all together
dams.Fgrm <- c(ELF.dams.Fgrm, PCC.dams.Fgrm)
sires.Fgrm <- c(ELF.sires.Fgrm, PCC.sires.Fgrm)

dams.SVL <- c(ELF.dams.SVL, PCC.dams.SVL)
sires.SVL <- c(ELF.sires.SVL, PCC.sires.SVL)

par(mfrow=c(1,2))
plot(dams.SVL~dams.Fgrm, main = "Dams", xlab = "Fgrm", ylab = "average adult SVL", pch = 19, col = colors[4])
abline(lm(dams.SVL~dams.Fgrm), lty = 2, lwd = 2)
summary(lm(dams.SVL~dams.Fgrm))

plot(sires.SVL[-20]~sires.Fgrm[-20], main = "Sires", xlab = "Fgrm", ylab = "average adult SVL", pch = 19, col = colors[4])
abline(lm(sires.SVL[-20]~sires.Fgrm[-20]), lty = 2, lwd = 2)
summary(lm(sires.SVL[-20]~sires.Fgrm[-20]))

## Max SVL
extract_max_SVL <- function(inds, exp_metaData){
  # inds <- names(ELF.dams)
  # data_type <- "SVL"
  # exp_metaData <- ELF_exp_metaData
  max_SVL <- vector(length = length(inds))
  for(i in 1:length(inds)){
    snake_data <- exp_metaData[[inds[i]]]
    if(any(snake_data$age.class == "A")){
      max_SVL[i] <- max(suppressWarnings(as.numeric(snake_data[which(snake_data$age.class == "A"),]$SVL)), na.rm = TRUE)
    }else{
      max_SVL[i] <- NA
    }
  }
  return(max_SVL)
}

ELF.dams.max.SVL <- extract_max_SVL(names(ELF.dams), ELF_exp_metaData)
ELF.sires.max.SVL <- extract_max_SVL(names(ELF.sires), ELF_exp_metaData)

PCC.dams.max.SVL <- extract_max_SVL(names(PCC.dams), PCC_exp_metaData)
PCC.sires.max.SVL <- extract_max_SVL(names(PCC.sires), PCC_exp_metaData)

dams.max.SVL <- c(ELF.dams.max.SVL, PCC.dams.max.SVL)
sires.max.SVL <- c(ELF.sires.max.SVL, PCC.sires.max.SVL)

par(mfrow=c(1,2))
plot(dams.max.SVL~dams.Fgrm, main = "Dams", xlab = "Fgrm", ylab = "average adult SVL", pch = 19, col = colors[4])
abline(lm(dams.max.SVL~dams.Fgrm), lty = 2, lwd = 2)
summary(lm(dams.max.SVL~dams.Fgrm))

plot(sires.max.SVL[-20]~sires.Fgrm[-20], main = "Sires", xlab = "Fgrm", ylab = "average adult SVL", pch = 19, col = colors[4])
abline(lm(sires.max.SVL[-20]~sires.Fgrm[-20]), lty = 2, lwd = 2)
summary(lm(sires.max.SVL[-20]~sires.Fgrm[-20]))

max.SVL <- c(sires.max.SVL, dams.max.SVL)
avg.SVL <- c(sires.SVL, dams.SVL)
Fgrm <- c(sires.Fgrm, dams.Fgrm)

par(mfrow=c(1,2))
plot(max.SVL[-20]~Fgrm[-20], main = "All", xlab = "Fgrm", ylab = "maximum adult SVL", pch = 19, col = colors[4])
abline(lm(max.SVL[-20]~Fgrm[-20]), lty = 2, lwd = 2)
summary(lm(max.SVL[-20]~Fgrm[-20]))

plot(avg.SVL[-20]~Fgrm[-20], main = "All", xlab = "Fgrm", ylab = "average adult SVL", pch = 19, col = colors[4])
abline(lm(avg.SVL[-20]~Fgrm[-20]), lty = 2, lwd = 2)
summary(lm(avg.SVL[-20]~Fgrm[-20]))

# keep outlier
par(mfrow=c(1,2))
plot(max.SVL~Fgrm, main = "All", xlab = "Fgrm", ylab = "maximum adult SVL", pch = 19, col = colors[4])
abline(lm(max.SVL~Fgrm), lty = 2, lwd = 2)
summary(lm(max.SVL~Fgrm))

plot(avg.SVL~Fgrm, main = "All", xlab = "Fgrm", ylab = "average adult SVL", pch = 19, col = colors[4])
abline(lm(avg.SVL~Fgrm), lty = 2, lwd = 2)
summary(lm(avg.SVL~Fgrm))

# ------------------------------------------------------------------------------------------------------------

## 15.75: Fgrm vs recaps -----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
extract_no_recaps <- function(inds, exp_metaData){
  # inds <- names(ELF.dams)
  # data_type <- "SVL"
  # exp_metaData <- ELF_exp_metaData
  no_recaps <- vector(length = length(inds))
  for(i in 1:length(inds)){
    snake_data <- exp_metaData[[inds[i]]]
    no_recaps[i] <- length(unique(snake_data$year))
  }
  return(no_recaps)
}

ELF_recaps <- extract_no_recaps(pop_gen_stats[which(pop_gen_stats$site == "ELF"), "id"], ELF_exp_metaData)
PCC_recaps <- extract_no_recaps(pop_gen_stats[which(pop_gen_stats$site == "PCC"), "id"], PCC_exp_metaData)

ELF_Fgrm <- pop_gen_stats[which(pop_gen_stats$site == "ELF"), "Fgrm"]
PCC_Fgrm <- pop_gen_stats[which(pop_gen_stats$site == "PCC"), "Fgrm"]

par(mfrow=c(1,2))
plot(ELF_recaps~ ELF_Fgrm, main = "All", xlab = "Fgrm", ylab = "number of years recaptured", pch = 19, col = colors[4])
abline(lm(ELF_recaps~ ELF_Fgrm), lty = 2, lwd = 2)
summary(lm(ELF_recaps~ ELF_Fgrm))

plot(PCC_recaps~PCC_Fgrm, main = "All", xlab = "Fgrm", ylab = "number of years recaptured", pch = 19, col = colors[4])
abline(lm(PCC_recaps~PCC_Fgrm), lty = 2, lwd = 2)
summary(lm(PCC_recaps~PCC_Fgrm))

# together
recaps <- c(ELF_recaps, PCC_recaps)
Fgrm_recaps <- c(ELF_Fgrm, PCC_Fg)
par(mfrow=c(1,1))
plot(recaps~Fgrm_recaps, main = "All", xlab = "Fgrm", ylab = "number of years recaptured", pch = 19, col = colors[4])
abline(lm(recaps~Fgrm_recaps), lty = 2, lwd = 2)
summary(lm(recaps~Fgrm_recaps))

# ------------------------------------------------------------------------------------------------------------

## 16: The age problem -----------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

CalcBYprobs(Pedigree = ELF_pedigree, LifeHistData=ELF_meta)

ELF.age.prior <- MakeAgePrior(Pedigree = ELF_pedigree, LifeHistData = ELF_meta, MinAgeParent = 1, Plot=T, Return = "all")

PlotAgePrior(ELF.age.prior$LR.RU.A)

ELF.by.probs <- CalcBYprobs(Pedigree = ELF_pedigree, LifeHistData = ELF_meta)

ELF.BirthYears <- cbind.data.frame(ELF_meta$ID, as.numeric(ELF_meta$BirthYear))
colnames(ELF.BirthYears) <- c("id", "BirthYear")

for(i in 1:nrow(ELF.by.probs)){
  ind <- rownames(ELF.by.probs)[i]
  if(sum(ELF.by.probs[i,]) > 0){
    if(length(which(ELF.by.probs[i,]==max(ELF.by.probs[i,]))) == 1){
      by <- as.numeric(colnames(ELF.by.probs)[which(ELF.by.probs[i,]==max(ELF.by.probs[i,]))])
      ELF.BirthYears[which(ELF.BirthYears$id == ind),"BirthYear"] <- by
    }else if((length(which(ELF.by.probs[i,]==max(ELF.by.probs[i,]))) %% 2) == 1){ # odd numbers
      by <- median(as.numeric(colnames(ELF.by.probs)[which(ELF.by.probs[i,]==max(ELF.by.probs[i,]))]))
      ELF.BirthYears[which(ELF.BirthYears$id == ind),"BirthYear"] <- by
    }else if((length(which(ELF.by.probs[i,]==max(ELF.by.probs[i,]))) %% 2) == 0){ # even numbers
      by <- mean(as.numeric(colnames(ELF.by.probs)[which(ELF.by.probs[i,]==max(ELF.by.probs[i,]))]))
      ELF.BirthYears[which(ELF.BirthYears$id == ind),"BirthYear"] <- by
    }
  }
}

# Is the relationship between number of offspring and SVL of ELF dams still significant if BirthYear is a effect in the model 
ELF.dams.by.svl <- cbind.data.frame(names(ELF.dams), as.vector(ELF.dams), ELF.dams.SVL)
colnames(ELF.dams.by.svl) <- c("id", "no_offspring", "SVL")

ELF.dams.by.svl$age_est <- 2018 - ELF.dams.by.svl$BirthYear
ELF.dams.by.svl[which(ELF.dams.by.svl$age_est > 9),"age_est"] <- 9 # capping ages at 9 

ELF.dams.by.svl$BirthYear <- ELF.BirthYears[match(ELF.dams.by.svl$id, ELF.BirthYears$id),"BirthYear"]
ELF.dams.by.svl <- ELF.dams.by.svl[which(ELF.dams.by.svl$BirthYear != -9),]

plot(ELF.dams.by.svl$no_offspring~ELF.dams.by.svl$age_est)

mod <- lm(no_offspring ~ SVL + age_est, data = ELF.dams.by.svl)
summary(mod)

# okay, what about heterozygosit? 
# using "het" from filter_vcf_popgen.R
ELF.dams.by.svl.het <- ELF.dams.by.svl
ELF.dams.by.svl.het$het <- het[match(ELF.dams.by.svl$id, names(het))]

plot(no_offspring~het, data = ELF.dams.by.svl.het)

mod1 <- lm(no_offspring ~ het + age_est + SVL, data = ELF.dams.by.svl.het)
mod2 <- lm(no_offspring ~ age_est + SVL, data = ELF.dams.by.svl.het)
AICtab(mod1, mod2)

plot(no_offspring~age_est, data = ELF.dams.by.svl.het) # I should include all individuals we have heterozygosity for, including thoes with no offpsring

hist(het[which(grepl("ELF", names(het)))], breaks = 50)
abline(v= mean(het[which(grepl("ELF", names(het)))]), col = "black", lty = 2)
hist(ELF.dams.by.svl.het$het, add = T, col = "red")
abline(v= mean(ELF.dams.by.svl.het$het), col = "red", lty = 2)

ELF.BirthYears[match(names(ELF.dams), ELF.BirthYears$id),]

ELF.dams.by <- cbind.data.frame(ELF.BirthYears[match(names(ELF.dams), ELF.BirthYears$id),], as.vector(ELF.dams))

ELF.dams.by <- ELF.dams.by[which(ELF.dams.by$BirthYear != -9),]
colnames(ELF.dams.by) <- c("id", "BirthYear", "offspring")




# ------------------------------------------------------------------------------------------------------------

# other questions: 
# [x] average number of offspring
# [x] average number of mates
# [] average number of reproductive years
# [] is the number of offspring from a parent correlated with offspring being captured as an adult?
# [] do sibling groups have similar recapture probabilities? 
# [] are there changes in reproductive success over time? 

# to do: 
# revise metadata and record number of recaptures and dates of recaptures 
# also capture GPS coords 

# make lists of sibilings 
# 








## load libraries
library(sequoia)
library(kinship2)
library(sp)
library(sf)
#library(rgdal)
library(ggmap)
#library(ggsn)
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

## Define custom functions ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
percent_in_ped <- function(pedigree, plot = TRUE, legend = TRUE){
  inds <- pedigree$id
  site <- unlist(strsplit(inds, split = "_"))[1]
  geno.inds <- inds[grepl(site, inds)]
  
  parents.mia <- 0
  both.parents <- 0
  one.parent <- 0
  
  for(i in 1:length(geno.inds)){
    parents <- pedigree[i,2:3]
    if(sum(is.na(parents)) == 0){ # both parents called
      both.parents <- both.parents + 1
    }else if(sum(is.na(parents)) == 1){ # one parent called
      one.parent <- one.parent + 1
    }else if(sum(is.na(parents)) == 2){ # both parents NA
      parents.mia <- parents.mia + 1 
    }
  }
  print(paste0(round(parents.mia / length(geno.inds), 2), " of individuals have no parents assigned")) 
  print(paste0(round(both.parents / length(geno.inds), 2), " of individuals have both parents assigned (genotyped or dummy)"))
  print(paste0(round(one.parent / length(geno.inds), 2), " of individuals have one parent assigned (genotyped or dummy)"))
  
  if(plot == TRUE){
    
    percents <- sapply(c(parents.mia, one.parent, both.parents), function(x){x*100/sum(c(parents.mia, one.parent, both.parents))})
    barplot(as.matrix(percents), 
            xlim = c(0,4), 
            ylim = c(0,100), 
            col = colors[c(1,3,2)], 
            border = NA, 
            main = paste0("Parents assigned at ", site))
    if(legend == TRUE){
      legend("right", legend = c("No parents assigned", "One parent assigned", "Both parents assigned"), 
             pch = 15, col = colors[c(1,3,2)], 
             bty="n")
    }
  }
  
  return(list(length(geno.inds), parents.mia, both.parents, one.parent))
}

offspring_per_ind <- function(pedigree, plot = TRUE){
  # testing vars
  #  pedigree <-  PCC_pedigree
  
  # start function
  inds <- pedigree$id
  site <- unlist(strsplit(inds, split = "_"))[1]
  
  dams <- pedigree$dam
  sires <- pedigree$sire
  
  u.dams <- na.omit(unique(dams)) # vector of unique dams 
  u.sires <- na.omit(unique(sires))  # vector of unique sires 
  # u.dams <- u.dams[grepl(site, u.dams)]   # remove dummy individuals
  
  offspring.no.dams <- table(dams)
  offspring.no.sires <- table(sires)
  
  if(plot == TRUE){
    par(mfrow=c(1,2))
    hist(offspring.no.dams, 
         xlab = "number of offspring", 
         main = paste0(site, " dams"), col = colors[4])
    mtext("includes dummy individuals")
    hist(offspring.no.sires, 
         xlab = "number of offspring", 
         main = paste0(site, " sires"), col = colors[5])
    mtext("includes dummy individuals")
  } # end plots
  return(list(offspring.no.dams, offspring.no.sires)) 
} # end function

find_super_reproducers <- function(offspring.list, site, meta.data){
  top.repro <- vector(mode = "list", length = 2)
  for(i in 1:2){
    inds <- offspring.list[[i]]
    inds <- inds[grepl(site, names(inds))] # remove dummy individuals
    top.repro[[i]] <- inds[order(inds)][round(length(inds)*.90):length(inds)] # get top 10% of reproducers? 
  }
  # find metadata
  data <- unlist(top.repro)
  top.repro.df <- data.frame(matrix(data = NA, nrow = length(data), ncol = 4))
  colnames(top.repro.df) <- c("ID", "sex", "offspring.no", "BY.max")
  for(ind in 1:length(names(data))){
    top.repro.df$ID[ind] <- names(data)[ind]
    top.repro.df$sex[ind] <- meta.data[which(meta.data$ID == names(data)[ind]),"Sex"]
    top.repro.df$offspring.no[ind] <- data[ind]
    top.repro.df$BY.max[ind] <- as.integer(meta.data[which(meta.data$ID == names(data)[ind]),"BY.max"])
  }
  return(top.repro.df)
}

mate_no <- function(pedigree){
  inds <- pedigree$id
  site <- unlist(strsplit(inds, split = "_"))[1]
  inds <- inds[grepl(site, inds)]
  
  mates <- pedigree[,2:3]
  mates <- mates[!is.na(mates$dam) | !is.na(mates$sire),] # get rid of when both parents are NAs
  percent_unique <- dim(unique(mates))[1] / dim(mates)[1]
  print(paste0(percent_unique, " percent of pairings are unique"))
  
  no.mates <- data.frame(matrix(data=NA, nrow = length(inds), ncol = 3))
  colnames(no.mates) <- c("ID", "no.mates", "no.offspring")
  for (i in 1:length(inds)){
    no.mates[i,]$ID <- inds[i] 
    if(!is.na(any(mates == inds[i]))){
      no.mates[i,]$no.mates <- nrow(unique(mates[which(mates == inds[i], arr.ind=TRUE)[,1],]))
      no.mates[i,]$no.offspring <- nrow(mates[which(mates == inds[i], arr.ind=TRUE)[,1],])
    }else{
      no.mates[i,]$no.mates <- 0
      no.mates[i,]$no.offspring <- 0
    }
  }
  return(no.mates)
} # end function

get_coords <- function(exp_metaData, target_inds, 
                       easting_col, northing_col, 
                       age_col = "age.class", 
                       year_col = "year", 
                       doy = FALSE){
  # testing vars
  #   exp_metaData <- PCC_exp_metaData
  #   target_inds <-  PCC_pedigree$id
  
  # remove dummy inds
  site <- unlist(strsplit(target_inds, split = "_"))[1]
  target_inds <- target_inds[grepl(site, target_inds)] # removes dummy individuals if using a list from a pedigree
  
  # remove tech replicate names
  target_inds[which(grepl("_P", target_inds))] <- str_remove(target_inds[which(grepl("_P", target_inds))], "_P[0-9]")
  
  coords <- data.frame()
  count = 1
  for(i in 1:length(target_inds)){
    snake_data <- exp_metaData[[which(names(exp_metaData) == target_inds[i])]]
    if(dim(snake_data)[1] > 0){
      for(j in 1:nrow(snake_data)){
        coords[count,1] <- target_inds[i]
        coords[count,2] <- snake_data[,easting_col][j]
        coords[count,3] <- snake_data[,northing_col,][j]
        coords[count,4] <- snake_data[,age_col][j]
        coords[count,5] <- snake_data[,year_col][j]
        
        if(doy == TRUE){
          coords[count,6] <- snake_data[,"doy"][j]
          }
        
        count <- count + 1
      }
    }
  }
  colnames(coords) <- c("ID", "easting", "northing", "age", "year")
  coords$easting <- as.numeric(coords$easting)
  coords$northing <- as.numeric(coords$northing)
  
  return(coords)
}

find_children <- function(relatedness_matrix, individual){
  relationships <- relatedness_matrix[,which(colnames(relatedness_matrix) == individual)] # relationships relative to ind (ind is the ___ of rowname)
  children <- names(which(relationships == "MP"))
  return(children)
}

find_descendants <- function(seq_results, inds_list){
  # testing vars
  # seq_results <- PCC_results
  # inds_list <- PCC_cohorts[[1]]
  
  site <- unlist(strsplit(seq_results[[4]]$Pedigree$id, split = "_"))[1]
  rel.matrix <- GetRelM(seq_results[[4]]$Pedigree, GenBack = 1, Return = 'Matrix')
  
  # how many descendants does each ind have?
  family_line <- vector(mode = "list", length = length(inds_list))
  for (ind in inds_list){
    descendants <- c()
    children <- find_children(rel.matrix, ind)
    if(length(children > 0)){
      descendants <- c(descendants, children)
      for (child in children){
        grandkids <- find_children(rel.matrix, child)
        if(length(grandkids) > 0){
          descendants <- c(descendants, grandkids)
          for(grandkid in grandkids){
            ggkids <- find_children(rel.matrix, grandkid)
            if(length(ggkids) > 0){
              descendants <- c(descendants, ggkids)
            }else(print("the broken line"))
          }
        }else(print("the broken line"))
      }
    }
    if(!is.null(descendants)){
      family_line[[which(inds_list == ind)]] <- descendants
    } 
    names(family_line)[[which(inds_list == ind)]] <- ind
  } # end loop through inds_list
  return(family_line)
}

find_old_young_cohorts <- function(seq_results, percent.older, percent.younger, plot = TRUE){
  site <- unlist(strsplit(seq_results[[4]]$Pedigree$id, split = "_"))[1]
  
  by <- seq_results[[4]]$LifeHistSib$BY.est
  y.young <- max(by, na.rm=T) - round(length(seq(to = range(by, na.rm=T)[1], 
                                                 from = range(by, na.rm=T)[2]))*percent.younger)
  y.old <- min(by, na.rm=T) + round(length(seq(to = range(by, na.rm=T)[1], 
                                               from = range(by, na.rm=T)[2]))*percent.older)
  
  old.inds <- seq_results[[4]]$LifeHistSib[which(seq_results[[4]]$LifeHistSib$BY.est <= y.old),"id"]
  young.inds <- seq_results[[4]]$LifeHistSib[which(seq_results[[4]]$LifeHistSib$BY.est >= y.young),"id"]
  
  if(plot == TRUE){
    hist(by, 
         main = site, 
         xlab = "Estimated birth year", 
         col = colors[1])
    abline(v=y.young, col = "black", lty = 2, lwd = 2)
    abline(v=y.old, col = "black", lty = 2, lwd = 2)
    
  }
  return(list(old.inds, young.inds))
}

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
    theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
    labs(title=plot.title) + ylab("") + xlab("") + 
    scale_color_manual(values=met.brewer(palette, n = length(unique(dat.df$ID))))
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

plot_ELF_lines <- function(ind){
  #ind <- "ELF_135"
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
  
  #hist(coords_clean$easting)
  #coords_clean[which(coords_clean$easting == 583745),]
  if(ind == "ELF_135"){
    coords_clean <- coords_clean[-35,]
  } # quick fix, need to change in metadata
  
  #convert to lat long
  sputm <- SpatialPoints(coords_clean[c("easting", "northing")], proj4string=CRS("+proj=utm +zone=16 +datum=WGS84")) 
  spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
  
  dat <- SpatialPointsDataFrame(coords = spgeo, data = coords_clean[,c(1,4,5,7)])
  dat.df<-data.frame(dat)
  colnames(dat.df) <- c("ID", "age", "year", "status", "long", "lat")
  
  box = c(left = -86.01221, bottom = 41.94555, right = -85.98966, top = 41.95972)
  zoom = 13
  map.type = "satellite"
  map <- get_map(location = box, maptype = "satellite", source = "google")
  
  p <- ggmap(map) + geom_point(aes(x = long, y = lat, color = status), 
                               data = dat.df, size = 2.5) + 
    theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
    labs(title=paste0(ind, " and descendants")) + ylab("") + xlab("") +
    scale_color_manual(values=colors[c(4,5,6)])
   for(i in 1:length(dat.df$ID)){
     if(nrow(dat.df[which(dat.df$ID == dat.df$ID[i]),]) > 1){
       p <- p + geom_path(aes(x = long, y = lat, color = status), 
                          data = dat.df[which(dat.df$ID == dat.df$ID[i]),])
     }
   }
   p
  # p + scalebar(location = "bottomleft", dist = 200,
  #              dist_unit = "m", transform = TRUE, st.size = 3.25, 
  #              nudge_y = -0.0002, 
  #              x.min = box[1]*.99999, 
  #              x.max = box[3], 
  #              y.min = box[2]*1.00003, 
  #              y.max = box[4])
}

plot_PCC_lines <- function(ind){
  #ind <- "PCC_39"
  colors = met.brewer("Archambault", n = 6)
  descendants <- unlist(GetDescendants(ind, PCC_pedigree)[c(1,2,3)])
  coords <- PCC_coords[PCC_coords$ID %in% descendants,]
  coords[which(coords$ID %in% GetDescendants(ind, PCC_pedigree)$offspring),"status"] <- "offspring"
  coords[which(coords$ID %in% GetDescendants(ind, PCC_pedigree)$grandoffspring),"status"] <- "grandoffspring"
  coords[which(coords$ID %in% GetDescendants(ind, PCC_pedigree)$id),"status"] <- "self"
  
  # clean coords
  coords_clean <- subset(coords, subset = !coords[,"easting"] == ".")
  coords_clean[,2:3] <- sapply(coords_clean[,c("easting", "northing")], as.numeric)
  coords_clean <-  subset(coords_clean, subset = !is.na(coords_clean[,"easting"]))
  
  #convert to spatial object
  spgeo<- SpatialPoints(coords_clean[c("easting", "northing")], proj4string=CRS("+proj=longlat +datum=WGS84")) 
  
  dat <- SpatialPointsDataFrame(coords = spgeo, data = coords_clean[,c(1,4,5,6)])
  dat.df<-data.frame(dat)
  colnames(dat.df) <- c("ID", "age", "year", "status", "long", "lat")
  
  box = c(left = -85.301, bottom = 42.529, right = -85.284, top = 42.5476)
  zoom = 13
  map.type = "satellite"
  #map <- get_stamenmap(bbox = box, zoom = zoom, maptype = map.type)
  map <- get_map(location = box, maptype = "satellite", source = "google")
  
  p <- ggmap(map) + geom_point(aes(x = long, y = lat, color = status), 
                               data = dat.df, size = 2.5) + 
    theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
    labs(title=paste0(ind, " and descendants")) + ylab("") + xlab("") +
    scale_color_manual(values=colors[c(4,5,6)])

  for(i in 1:length(dat.df$ID)){
    if(nrow(dat.df[which(dat.df$ID == dat.df$ID[i]),]) > 1){
    p <- p + geom_path(aes(x = long, y = lat, color = status), 
        data = dat.df[which(dat.df$ID == dat.df$ID[i]),])
    }
  }
  p
}

infer_birth_year <- function(ind, location, data){
  # ind <- ELF_known_by[i,"ID"]
  # location = "cass"
  # data = ELF_SVL_data
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


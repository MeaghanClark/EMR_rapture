### Make figures for EMR inbreeding paper

## To do: make relM with dummy individuals to look at offspring vs ped inbreeding 

### Load packages 
# ------------------------------------------------------------------------------------------------------------
library(mapproj)
library(kinship2)
library(sequoia)
library(stringr)
library(RColorBrewer)
library(MetBrewer)
library(scales)
library(ggplot2)
library(gridExtra)

# ------------------------------------------------------------------------------------------------------------

plot.pedigree.mod <- function (x, id = x$id, status = x$status, affected = x$affected, 
                               cex = 1, col = 1, symbolsize = 1, branch = 0.6, packed = TRUE, 
                               align = c(1.5, 2), width = 8, density = c(-1, 35, 65, 20), 
                               angle = c(90, 65, 40, 0), keep.par = FALSE, subregion, pconnect = 0.5, 
                               ...) 
{
  Call <- match.call()
  n <- length(x$id)
  if (n < 3) {
    stop("Cannot plot pedigree with fewer than 3 subjects")
  }
  if (is.null(status)) 
    status <- rep(0, n)
  else {
    if (!all(status == 0 | status == 1)) 
      stop("Invalid status code")
    if (length(status) != n) 
      stop("Wrong length for status")
  }
  if (!missing(id)) {
    if (length(id) != n) 
      stop("Wrong length for id")
  }
  if (is.null(affected)) {
    affected <- matrix(0, nrow = n)
  }
  else {
    if (is.matrix(affected)) {
      if (nrow(affected) != n) 
        stop("Wrong number of rows in affected")
      if (is.logical(affected)) 
        affected <- 1 * affected
      if (ncol(affected) > length(angle) || ncol(affected) > 
          length(density)) 
        stop("More columns in the affected matrix than angle/density values")
    }
    else {
      if (length(affected) != n) 
        stop("Wrong length for affected")
      if (is.logical(affected)) 
        affected <- as.numeric(affected)
      if (is.factor(affected)) 
        affected <- as.numeric(affected) - 1
    }
    if (max(affected, na.rm = TRUE) > min(affected, na.rm = TRUE)) {
      affected <- matrix(affected - min(affected, na.rm = TRUE), 
                         nrow = n)
    }
    else {
      affected <- matrix(affected, nrow = n)
    }
    affected[is.na(affected)] <- -1
    if (!all(affected == 0 | affected == 1 | affected == 
             -1)) 
      stop("Invalid code for affected status")
  }
  if (length(col) == 1) 
    col <- rep(col, n)
  else if (length(col) != n) 
    stop("Col argument must be of length 1 or n")
  subregion2 <- function(plist, subreg) {
    if (subreg[3] < 1 || subreg[4] > length(plist$n)) 
      stop("Invalid depth indices in subreg")
    lkeep <- subreg[3]:subreg[4]
    for (i in lkeep) {
      if (!any(plist$pos[i, ] >= subreg[1] & plist$pos[i, 
      ] <= subreg[2])) 
        stop(paste("No subjects retained on level", 
                   i))
    }
    nid2 <- plist$nid[lkeep, ]
    n2 <- plist$n[lkeep]
    pos2 <- plist$pos[lkeep, ]
    spouse2 <- plist$spouse[lkeep, ]
    fam2 <- plist$fam[lkeep, ]
    if (!is.null(plist$twins)) 
      twin2 <- plist$twins[lkeep, ]
    for (i in 1:nrow(nid2)) {
      keep <- which(pos2[i, ] >= subreg[1] & pos2[i, ] <= 
                      subreg[2])
      nkeep <- length(keep)
      n2[i] <- nkeep
      nid2[i, 1:nkeep] <- nid2[i, keep]
      pos2[i, 1:nkeep] <- pos2[i, keep]
      spouse2[i, 1:nkeep] <- spouse2[i, keep]
      fam2[i, 1:nkeep] <- fam2[i, keep]
      if (!is.null(plist$twins)) 
        twin2[i, 1:nkeep] <- twin2[i, keep]
      if (i < nrow(nid2)) {
        tfam <- match(fam2[i + 1, ], keep, nomatch = 0)
        fam2[i + 1, ] <- tfam
        if (any(spouse2[i, tfam] == 0)) 
          stop("A subregion cannot separate parents")
      }
    }
    n <- max(n2)
    out <- list(n = n2[1:n], nid = nid2[, 1:n, drop = F], 
                pos = pos2[, 1:n, drop = F], spouse = spouse2[, 
                                                              1:n, drop = F], fam = fam2[, 1:n, drop = F])
    if (!is.null(plist$twins)) 
      out$twins <- twin2[, 1:n, drop = F]
    out
  }
  plist <- align.pedigree(x, packed = packed, width = width, 
                          align = align)
  if (!missing(subregion)) 
    plist <- subregion2(plist, subregion)
  xrange <- range(plist$pos[plist$nid > 0])
  maxlev <- nrow(plist$pos)
  frame()
  oldpar <- par(xpd = TRUE, ...)
  psize <- par("pin")
  stemp1 <- strwidth("ABC", units = "inches", cex = cex) * 
    2.5/3
  stemp2 <- strheight("1g", units = "inches", cex = cex)
  stemp3 <- max(strheight(id, units = "inches", cex = cex))
  ht1 <- psize[2]/maxlev - (stemp3 + 1.5 * stemp2)
  if (ht1 <= 0) 
    stop("Labels leave no room for the graph, reduce cex")
  ht2 <- psize[2]/(maxlev + (maxlev - 1)/2)
  wd2 <- 0.8 * psize[1]/(0.8 + diff(xrange))
  boxsize <- symbolsize * min(ht1, ht2, stemp1, wd2)
  hscale <- (psize[1] - boxsize)/diff(xrange)
  vscale <- (psize[2] - (stemp3 + stemp2/2 + boxsize))/max(1, 
                                                           maxlev - 1)
  boxw <- boxsize/hscale
  boxh <- boxsize/vscale
  labh <- stemp2/vscale
  legh <- min(1/4, boxh * 1.5)
  par(usr = c(xrange[1] - boxw/2, xrange[2] + boxw/2, maxlev + 
                boxh + stemp3 + stemp2/2, 1))
  circfun <- function(nslice, n = 50) {
    nseg <- ceiling(n/nslice)
    theta <- -pi/2 - seq(0, 2 * pi, length = nslice + 1)
    out <- vector("list", nslice)
    for (i in 1:nslice) {
      theta2 <- seq(theta[i], theta[i + 1], length = nseg)
      out[[i]] <- list(x = c(0, cos(theta2)/2), y = c(0, 
                                                      sin(theta2)/2) + 0.5)
    }
    out
  }
  polyfun <- function(nslice, object) {
    zmat <- matrix(0, ncol = 4, nrow = length(object$x))
    zmat[, 1] <- object$x
    zmat[, 2] <- c(object$x[-1], object$x[1]) - object$x
    zmat[, 3] <- object$y
    zmat[, 4] <- c(object$y[-1], object$y[1]) - object$y
    ns1 <- nslice + 1
    theta <- -pi/2 - seq(0, 2 * pi, length = ns1)
    x <- y <- double(ns1)
    for (i in 1:ns1) {
      z <- (tan(theta[i]) * zmat[, 1] - zmat[, 3])/(zmat[, 
                                                         4] - tan(theta[i]) * zmat[, 2])
      tx <- zmat[, 1] + z * zmat[, 2]
      ty <- zmat[, 3] + z * zmat[, 4]
      inner <- tx * cos(theta[i]) + ty * sin(theta[i])
      indx <- which(is.finite(z) & z >= 0 & z <= 1 & inner > 
                      0)
      x[i] <- tx[indx]
      y[i] <- ty[indx]
    }
    nvertex <- length(object$x)
    temp <- data.frame(indx = c(1:ns1, rep(0, nvertex)), 
                       theta = c(theta, object$theta), x = c(x, object$x), 
                       y = c(y, object$y))
    temp <- temp[order(-temp$theta), ]
    out <- vector("list", nslice)
    for (i in 1:nslice) {
      rows <- which(temp$indx == i):which(temp$indx == 
                                            (i + 1))
      out[[i]] <- list(x = c(0, temp$x[rows]), y = c(0, 
                                                     temp$y[rows]) + 0.5)
    }
    out
  }
  if (ncol(affected) == 1) {
    polylist <- list(square = list(list(x = c(-1, -1, 1, 
                                              1)/2, y = c(0, 1, 1, 0))), circle = list(list(x = 0.5 * 
                                                                                              cos(seq(0, 2 * pi, length = 50)), y = 0.5 * sin(seq(0, 
                                                                                                                                                  2 * pi, length = 50)) + 0.5)), diamond = list(list(x = c(0, 
                                                                                                                                                                                                           -0.5, 0, 0.5), y = c(0, 0.5, 1, 0.5))), triangle = list(list(x = c(0, 
                                                                                                                                                                                                                                                                              -0.56, 0.56), y = c(0, 1, 1))))
  }
  else {
    nc <- ncol(affected)
    square <- polyfun(nc, list(x = c(-0.5, -0.5, 0.5, 0.5), 
                               y = c(-0.5, 0.5, 0.5, -0.5), theta = -c(3, 5, 7, 
                                                                       9) * pi/4))
    circle <- circfun(nc)
    diamond <- polyfun(nc, list(x = c(0, -0.5, 0, 0.5), 
                                y = c(-0.5, 0, 0.5, 0), theta = -(1:4) * pi/2))
    triangle <- polyfun(nc, list(x = c(-0.56, 0, 0.56), 
                                 y = c(-0.5, 0.5, -0.5), theta = c(-2, -4, -6) * 
                                   pi/3))
    polylist <- list(square = square, circle = circle, diamond = diamond, 
                     triangle = triangle)
  }
  drawbox <- function(x, y, sex, affected, status, col, polylist, 
                      density, angle, boxw, boxh) {
    for (i in 1:length(affected)) {
      if (affected[i] == 0) {
        polygon(x + (polylist[[sex]])[[i]]$x * boxw, 
                y + (polylist[[sex]])[[i]]$y * boxh, col = NA, 
                border = col)
      }
      if (affected[i] == 1) {
        polygon(x + (polylist[[sex]])[[i]]$x * boxw, 
                y + (polylist[[sex]])[[i]]$y * boxh, col = col, 
                border = col, density = density[i], angle = angle[i])
      }
      if (affected[i] == -1) {
        polygon(x + (polylist[[sex]])[[i]]$x * boxw, 
                y + (polylist[[sex]])[[i]]$y * boxh, col = NA, 
                border = col)
        midx <- x + mean(range(polylist[[sex]][[i]]$x * 
                                 boxw))
        midy <- y + mean(range(polylist[[sex]][[i]]$y * 
                                 boxh))
        points(midx, midy, pch = "?", cex = min(1, cex * 
                                                  2/length(affected)))
      }
    }
    if (status == 1) 
      segments(x - 0.6 * boxw, y + 1.1 * boxh, x + 0.6 * 
                 boxw, y - 0.1 * boxh, )
  }
  sex <- as.numeric(x$sex)
  for (i in 1:maxlev) {
    for (j in seq_len(plist$n[i])) {
      k <- plist$nid[i, j]
      drawbox(plist$pos[i, j], i, sex[k], affected[k, 
      ], status[k], col[k], polylist, density, angle, 
      boxw, boxh)
      text(plist$pos[i, j], i + boxh + labh * 0.7, id[k], 
           cex = cex, adj = c(0.5, 1), ...)
    }
  }
  maxcol <- ncol(plist$nid)
  for (i in 1:maxlev) {
    tempy <- i + boxh/2
    if (any(plist$spouse[i, ] > 0)) {
      temp <- (1:maxcol)[plist$spouse[i, ] > 0]
      segments(plist$pos[i, temp] + boxw/2, rep(tempy, 
                                                length(temp)), plist$pos[i, temp + 1] - boxw/2, 
               rep(tempy, length(temp)))
      temp <- (1:maxcol)[plist$spouse[i, ] == 2]
      if (length(temp)) {
        tempy <- tempy + boxh/10
        segments(plist$pos[i, temp] + boxw/2, rep(tempy, 
                                                  length(temp)), plist$pos[i, temp + 1] - boxw/2, 
                 rep(tempy, length(temp)))
      }
    }
  }
  for (i in 2:maxlev) {
    zed <- unique(plist$fam[i, ])
    zed <- zed[zed > 0]
    for (fam in zed) {
      xx <- plist$pos[i - 1, fam + 0:1]
      parentx <- mean(xx)
      who <- (plist$fam[i, ] == fam)
      if (is.null(plist$twins)) 
        target <- plist$pos[i, who]
      else {
        twin.to.left <- (c(0, plist$twins[i, who])[1:sum(who)])
        temp <- cumsum(twin.to.left == 0)
        tcount <- table(temp)
        target <- rep(tapply(plist$pos[i, who], temp, 
                             mean), tcount)
      }
      yy <- rep(i, sum(who))
      segments(plist$pos[i, who], yy, target, yy - legh)
      if (any(plist$twins[i, who] == 1)) {
        who2 <- which(plist$twins[i, who] == 1)
        temp1 <- (plist$pos[i, who][who2] + target[who2])/2
        temp2 <- (plist$pos[i, who][who2 + 1] + target[who2])/2
        yy <- rep(i, length(who2)) - legh/2
        segments(temp1, yy, temp2, yy)
      }
      if (any(plist$twins[i, who] == 3)) {
        who2 <- which(plist$twins[i, who] == 3)
        temp1 <- (plist$pos[i, who][who2] + target[who2])/2
        temp2 <- (plist$pos[i, who][who2 + 1] + target[who2])/2
        yy <- rep(i, length(who2)) - legh/2
        text((temp1 + temp2)/2, yy, "?")
      }
      segments(min(target), i - legh, max(target), i - 
                 legh)
      if (diff(range(target)) < 2 * pconnect) 
        x1 <- mean(range(target))
      else x1 <- pmax(min(target) + pconnect, pmin(max(target) - 
                                                     pconnect, parentx))
      y1 <- i - legh
      if (branch == 0) 
        segments(x1, y1, parentx, (i - 1) + boxh/2)
      else {
        y2 <- (i - 1) + boxh/2
        x2 <- parentx
        ydelta <- ((y2 - y1) * branch)/2
        segments(c(x1, x1, x2), c(y1, y1 + ydelta, y2 - 
                                    ydelta), c(x1, x2, x2), c(y1 + ydelta, y2 - 
                                                                ydelta, y2))
      }
    }
  }
  arcconnect <- function(x, y) {
    xx <- seq(x[1], x[2], length = 15)
    yy <- seq(y[1], y[2], length = 15) + (seq(-7, 7))^2/98 - 
      0.5
    lines(xx, yy, lty = 2)
  }
  uid <- unique(plist$nid)
  for (id in unique(uid[uid > 0])) {
    indx <- which(plist$nid == id)
    if (length(indx) > 1) {
      tx <- plist$pos[indx]
      ty <- ((row(plist$pos))[indx])[order(tx)]
      tx <- sort(tx)
      # for (j in 1:(length(indx) - 1)) arcconnect(tx[j + 
      #                                                 0:1], ty[j + 0:1])
    }
  }
  ckall <- x$id[is.na(match(x$id, x$id[plist$nid]))]
  if (length(ckall > 0)) 
    cat("Did not plot the following people:", ckall, "\n")
  if (!keep.par) 
    par(oldpar)
  tmp <- match(1:length(x$id), plist$nid)
  invisible(list(plist = plist, x = plist$pos[tmp], y = row(plist$pos)[tmp], 
                 boxw = boxw, boxh = boxh, call = Call))
}

source("./pedigree_funcs.R") 

get_dist_matrix <- function(relM, distM, relationship){
  if(sum(colnames(relM) != colnames(distM)) > 0){
    print("the names of the relatedness matrix do not match the distance matrix!")
  }else{
    pairs <- which(relM == relationship, arr.ind = TRUE)
    dist <- vector(length = nrow(pairs))
    for(i in 1:nrow(pairs)){
      dist[i] <- distM[pairs[i,][1],pairs[i,][2]]
    }
    return(dist)
  }
}

est_age_from_curve <- function(L, Loo, k = 0.378, Lo = 19.24){
  # function returns an age estimate based on a given SVL and asymptotic size
  if(L < Loo){
    return(-(log((Loo-L)/(Loo-Lo))) / (k))
  }else{
    return(NULL) 
  }
}

infer_birth_year <- function(ind, location, ind_SVL, ind_year, ind_sex){
  # function for inferring birth year from SVL based on growth curve
  ind_year <- as.numeric(ind_year)
  if(ind_sex == 1){
    ind_sex <- "female"
  }else if(ind_sex == 2){
    ind_sex <- "male"
  }else{
    #print("missing sex information")
    ind_sex <- NA
  }
  if(!is.na(ind_sex)){
    ind_age <- est_age_from_curve(L = ind_SVL, Loo = subset(asymptotic_L, site == location & sex == ind_sex, select = L_est))
    
    if(is.null(ind_age)){
      ind_age = NA  # return NA if individual was above asymptotic size
    }
    return(as.numeric(ind_year - ind_age))
  }else{
    return(NA) #  return NA if no sex information
  }
}

colFunc <- function (x, cols, nCols, valRange) {
  if (is.null(valRange)) {
    valRange <- c(min(x), max(x))
  }
  cols <- (grDevices::colorRampPalette(cols))(nCols)[findInterval(x, 
                                                                  seq(valRange[1], valRange[2], length.out = nCols))]
  return(cols)
}

### Load data objects 
# ------------------------------------------------------------------------------------------------------------

# [1] load popgen data from vcf_filtering_popgen.R
load(file = "../vcf_filtering/pop_gen_stats_02082024.Robj", verbose = T)
pop_gen_stats$id <- str_remove(pop_gen_stats$id, "_P[0-9]+")

load(file = "../vcf_filtering/pwp_02082024.Robj", verbose = T)
PCC_PWP <- PWP[[1]]
ELF_PWP <- PWP[[2]]

load(file = "../vcf_filtering/PCA_loadings_05212024.Robj", verbose = T)
ELF_pca <- pca[[1]]
PCC_pca <- pca[[2]]

load(file = "../vcf_filtering/joint_pca.RObj", verbose = T)

# [2] load pedigrees from make_pedigrees.R
load("../pedigree_reconstruction/ELF_pedigree_results_06262024.Robj", verbose = TRUE)
load("../pedigree_reconstruction/PCC_pedigree_results_06262024.Robj", verbose = TRUE)

PCC_pedigree <- PCC_results[[4]]$Pedigree
PCC_pedigree$id <- str_remove(PCC_pedigree$id, "_P[0-9]+")
PCC_results[[2]]$ID <- str_remove(PCC_results[[2]]$ID, "_P[0-9]+")
PCC_results[[2]] <- PCC_results[[2]][-286,]
PCC_meta <- PCC_results[[2]]

ELF_pedigree <- ELF_results[[4]]$Pedigree
ELF_pedigree$id <- str_remove(ELF_pedigree$id, "_P[0-9]+")
ELF_results[[2]] <- ELF_results[[2]][-779,]
ELF_results[[2]]$ID <- str_remove(ELF_results[[2]]$ID, "_P[0-9]+")
ELF_meta <- ELF_results[[2]]

load("../pedigree_reconstruction/elf_conf_06262024.robj", verbose = T) 
load("../pedigree_reconstruction/pcc_conf_06262024.robj", verbose = T)

# [3] load capture metadata 
load(paste0("../pedigree_exploration/ELF_expanded_metaData.Robj"), verbose = T) # ELF_exp_metaData
load(paste0("../pedigree_exploration/PCC_expanded_metaData.Robj"), verbose = T) # PCC_exp_metaData

# [4] read in raw data to get mom info for ELF pedigree filtering
ELF_rawData <- read.csv("~/Desktop/EMR_rapture/pedigree_reconstruction/ELF_LifeHistData.csv")[,1:45]

# [5] Load birth year estimates and years contributing from DRB

load("../pedigree_reconstruction/DRB_byEstimates.Robj", verbose = T) # export_temp

# make data frame of asymptotic size values
# asymptotic_L <- data.frame(matrix(NA, nrow = 4, ncol = 5))
# colnames(asymptotic_L) <- c("site", "sex", "lower_CI", "upper_CI", "L_est")
# asymptotic_L[1,] <- c("cass", "female", 57.5, 59.8, 58.7) # 58.7 cm (95% CI = 57.5 – 59.8) 
# asymptotic_L[2,] <- c("cass", "male", 60.3, 63.3, 61.8) # 61.8 cm (95% CI = 60.3 – 63.3)
# asymptotic_L[3,] <- c("barry", "female", 58.9, 61.6, 60.3) # 60.3 cm (95% CI = 58.9 – 61.6)
# asymptotic_L[4,] <- c("barry", "male", 61.8, 65.1, 63.4) # 63.4 cm (95% CI = 61.8 – 65.1) 
# 
# asymptotic_L$lower_CI <- as.numeric(asymptotic_L$lower_CI)
# asymptotic_L$upper_CI <- as.numeric(asymptotic_L$upper_CI)
# asymptotic_L$L_est <- as.numeric(asymptotic_L$L_est)

# ------------------------------------------------------------------------------------------------------------

### define color palette for figures
# ------------------------------------------------------------------------------------------------------------

palette = "Archambault"
met.brewer(palette, n=6, type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = 6)

# ------------------------------------------------------------------------------------------------------------

### Calculate spatial things: coordinates, centroid coordinates, pairwise distances  
# ------------------------------------------------------------------------------------------------------------
## coordinates

PCC_coords <- get_coords(exp_metaData =PCC_exp_metaData, target_inds =subset(PCC_pedigree, grepl("PCC", id), select = id)$id, easting_col ="DD.easting", northing_col ="DD.northing", doy = FALSE)
ELF_coords <- get_coords(exp_metaData =ELF_exp_metaData, target_inds =subset(ELF_pedigree, grepl("ELF", id), select = id)$id, easting_col ="UTM.easting", northing_col ="UTM.northing", doy = TRUE)

# quality filter coordinates: 
PCC_coords_flt <- na.omit(PCC_coords)
PCC_coords_flt$easting[PCC_coords_flt$easting > 0] <- PCC_coords_flt$easting[PCC_coords_flt$easting > 0]*-1 # easting should be negative 

PCC_coords_flt <- PCC_coords_flt[-which(PCC_coords_flt$easting < -85.3583),] # PCC_268, outside study bounds
PCC_coords_flt <- PCC_coords_flt[-which(PCC_coords_flt$northing < 42.518),] # PCC_303, outside study bounds

colnames(ELF_coords) <- c("ID", "easting", "northing", "age", "year", "doy")
ELF_coords_flt <- na.omit(ELF_coords)
ELF_coords_flt <- ELF_coords_flt[-which(ELF_coords_flt$easting > 583500),] # ELF_361, coord is ecolab

## centroid coordinates

# make spatial objects from coords
PCC_sf <- st_as_sf(PCC_coords_flt, coords = c("easting", "northing"))
PCC_sf <- st_set_crs(PCC_sf, "+proj=longlat +datum=WGS84")
plot(PCC_sf)
max(st_distance(PCC_sf)) # max distance between individual captures at BARRY 2042.665


ELF_sf <- st_as_sf(ELF_coords_flt, coords = c("easting", "northing"))
ELF_sf <- st_set_crs(ELF_sf, "+proj=utm +zone=16 +datum=WGS84")
plot(ELF_sf)
max(st_distance(ELF_sf)) # max distance between individual captures at CASS: 1725.077 m 

# covert to multipoints when appropriate (multiple capture events)
PCC_multi_sf <- aggregate(PCC_sf["geometry"], by = list(PCC_sf$ID), FUN = function(x) {
  st_multipoint(as.matrix(x))
})

ELF_multi_sf <- aggregate(ELF_sf["geometry"], by = list(ELF_sf$ID), FUN = function(x) {
  st_multipoint(as.matrix(x))
})

# Calculate centroids for each group
PCC_centroids <- st_centroid(PCC_multi_sf)
ELF_centroids <- st_centroid(ELF_multi_sf)

## Calculate pairwise distances between individuals 
PCC_centroid_dist <- st_distance(PCC_centroids) # dim: 256, 256
colnames(PCC_centroid_dist) <- PCC_centroids$Group.1
rownames(PCC_centroid_dist) <- PCC_centroids$Group.1

ELF_centroid_dist <- st_distance(ELF_centroids) # dim: 621, 621 
colnames(ELF_centroid_dist) <- ELF_centroids$Group.1
rownames(ELF_centroid_dist) <- ELF_centroids$Group.1

# ------------------------------------------------------------------------------------------------------------

### Calculate pedigree things: number of offspring, relatedness matrix, estimate birth year from SVL 
# ------------------------------------------------------------------------------------------------------------
## filter pedigree relationships at both sites to remove types of relationships that have low confidence probabilities 

# at ELF: genotyped dam, no sire, dummy focal (0.842 confidence probability) 
subset(ELF_pedigree, ((grepl("M", id) | grepl("F0", id)) & grepl("ELF", dam) & is.na(sire))) # 3 relationships to remove
ELF_pedigree_flt <- subset(ELF_pedigree, !((grepl("M", id) | grepl("F0", id)) & grepl("ELF", dam) & is.na(sire))) # 887 rows 

# genotyped sire, dummy dam, dummy focal (0 confidence probability )
subset(ELF_pedigree_flt, ((grepl("M", id) | grepl("F0", id)) & grepl("ELF", sire) & grepl("F0", dam))) # 3 relationships to remove
ELF_pedigree_flt <- subset(ELF_pedigree_flt, !((grepl("M", id) | grepl("F0", id)) & grepl("ELF", sire) & grepl("F0", dam))) # 884 rows 

# at PCC: genotyped dam, no sire, dummy focal (0.25 confidence probability) 
subset(PCC_pedigree, ((grepl("M", id) | grepl("F0", id)) & grepl("PCC", dam) & is.na(sire))) # 2 relationships to remove
PCC_pedigree_flt <- subset(PCC_pedigree, !((grepl("M", id) | grepl("F0", id)) & grepl("PCC", dam) & is.na(sire))) # 307 rows 

## Additional filtering: filter ELF pedigree to remove individuals born in captivity and not captured again 
# code  creates list where the list name is the name of the captive mom, and the contents of the list entry are the years she was held for partuition
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
captive.kids <- vector(mode = "list", length = length(captive.inds))
names(captive.kids) <- names(captive.inds)

for(i in 1:length(captive.inds)){
  captive.kids[[i]] <- ELF_rawData[grepl(strsplit(names(captive.inds)[i], "_")[[1]][2], ELF_rawData$Other.Notes) & (grepl("mom", ELF_rawData$Other.Notes, ignore.case = TRUE) | grepl("mother", ELF_rawData$Other.Notes, ignore.case = TRUE)),"ELF.ID"]
}

# filter didn't quite work for ELF_200
captive.kids[["ELF_200"]] <- captive.kids[["ELF_200"]][1:8]
# ELF_513 was held in lab for partuition in 2013, but passed 9-10 slugs, no neonates
captive.kids <- captive.kids[-which(names(captive.kids) == "ELF_513")]
captive.inds <- captive.inds[-which(names(captive.inds) == "ELF_513")]

# remove slug record
captive.kids[["ELF_161"]] <- captive.kids[["ELF_161"]][-6]

# make vector of their offspring determine which kids were only recorded as neonates in the lab and not recaptured

only_n_caps<- vector()
mult_caps<- vector() # individuals that have multiple captures (and thus were captured after birth, or were captured as adults)
need_help <- vector()
for(i in 1:length(captive.kids)){
  for (j in 1:length(captive.kids[[i]])){
    ind <- paste0("ELF_", captive.kids[[i]][j])
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

ELF_pedigree_flt2 <- subset(ELF_pedigree_flt, subset = !ELF_pedigree$id %in% only_n_caps) # 739 rows

## Calculate number of offspring using filtered pedigrees 
PCC_off <- offspring_per_ind(PCC_pedigree_flt, plot = T)
ELF_off <- offspring_per_ind(ELF_pedigree_flt2, plot = T)

## make relatedness matrix 
# using unfiltered pedigree here because we still have spatial information for some captive born individuals that should be informative
ELF_relM <- GetRelM(Pedigree = ELF_pedigree_flt, GenBack = 2, Return = "Matrix", patmat = TRUE)
ELF_relM <- ELF_relM[grepl("ELF", colnames(ELF_relM)),grepl("ELF", colnames(ELF_relM))] # remove dummy individuals

PCC_relM <- GetRelM(Pedigree = PCC_pedigree_flt, GenBack = 2, Return = "Matrix", patmat = TRUE)
PCC_relM <- PCC_relM[grepl("PCC", colnames(PCC_relM)),grepl("PCC", colnames(PCC_relM))] # remove dummy individuals


## make kinship matrix 
ELF_kinM <- CalcRped(ELF_pedigree_flt, OUT = 'M')
PCC_kinM <- CalcRped(PCC_pedigree_flt, OUT = 'M')

## calculate distances between parent/offspring and parent/parents

# filter individuals without spatial data from relM
PCC_relM_flt <- PCC_relM[match(colnames(PCC_centroid_dist), colnames(PCC_relM)), match(colnames(PCC_centroid_dist), colnames(PCC_relM))]
ELF_relM_flt <- ELF_relM[match(colnames(ELF_centroid_dist), colnames(ELF_relM)), match(colnames(ELF_centroid_dist), colnames(ELF_relM))]

# dam-offspring
PCC_dam_off_dist <- get_dist_matrix(relM = PCC_relM_flt, distM = PCC_centroid_dist, relationship = "M")
ELF_dam_off_dist <- get_dist_matrix(relM = ELF_relM_flt, distM = ELF_centroid_dist, relationship = "M")

# sire-offspring
PCC_sire_off_dist <- get_dist_matrix(relM = PCC_relM_flt, distM = PCC_centroid_dist, relationship = "P")
ELF_sire_off_dist <- get_dist_matrix(relM = ELF_relM_flt, distM = ELF_centroid_dist, relationship = "P")

# parent-parent
PCC_parents <- na.omit(PCC_pedigree_flt[,c("dam", "sire")])
PCC_unique_rents <- unique(PCC_parents)
PCC_unique_rents <- PCC_unique_rents[grepl("PCC_", PCC_unique_rents[,1]) & grepl("PCC_", PCC_unique_rents[,2]),] # no dummy individuals

PCC_parent_dist <- vector(length = nrow(PCC_unique_rents))
for(i in 1:nrow(PCC_unique_rents)){
  tryCatch(
    {
      PCC_parent_dist[i] <-  PCC_centroid_dist[which(rownames(PCC_centroid_dist) == PCC_unique_rents[i,1]),
                                           which(colnames(PCC_centroid_dist) == PCC_unique_rents[i,2])]
    }, 
    error = function(e){
      cat("Error occured at interation", i, ": ", conditionMessage(e), "\n")
    }
  )
}

# filter dummy individuals 
ELF_pedigree_flt_nD <- ELF_pedigree_flt[grepl("ELF", ELF_pedigree_flt$id),] # remove dummy individuals

ELF_parents <- na.omit(ELF_pedigree_flt_nD[,c("dam", "sire")]) # ELF_pedigree_flt_nD not defined yet! 
ELF_unique_rents <- unique(ELF_parents)
ELF_unique_rents <- ELF_unique_rents[grepl("ELF_", ELF_unique_rents[,1]) & grepl("ELF_", ELF_unique_rents[,2]),] # no dummy individuals

ELF_parent_dist <- vector(length = nrow(ELF_unique_rents))
for(i in 1:nrow(ELF_unique_rents)){
  tryCatch(
    {
      ELF_parent_dist[i] <-  ELF_centroid_dist[which(rownames(ELF_centroid_dist) == ELF_unique_rents[i,1]),
                                               which(colnames(ELF_centroid_dist) == ELF_unique_rents[i,2])]
    }, 
    error = function(e){
      ELF_parent_dist[i] <- NA
      cat("Error occured at interation", i, ": ", conditionMessage(e), "\n")
      
    }
  )
}
ELF_parent_dist[c(79, 90)] <- NA
# errors thrown because of missing spatial data for ELF_409 and ELF_584
# both only captured as a neonate, no spatial coords in dataset

# distance between unrelated individuals with no offspring in pedigree

# add parents to relationships in relM so that these pairs are not considered "unrelated"
PCC_relM_flt_par <- PCC_relM_flt
for(i in 1:nrow(PCC_unique_rents)){
  first <- which(colnames(PCC_relM_flt) == PCC_unique_rents[i,1])
  second <- which(colnames(PCC_relM_flt) == PCC_unique_rents[i,2])
  PCC_relM_flt_par[first,second] <- "PAR"
  PCC_relM_flt_par[second,first] <- "PAR"
}

ELF_relM_flt_par <- ELF_relM_flt
for(i in 1:nrow(ELF_unique_rents)){
  first <- which(colnames(ELF_relM_flt) == ELF_unique_rents[i,1])
  second <- which(colnames(ELF_relM_flt) == ELF_unique_rents[i,2])
  ELF_relM_flt_par[first,second] <- "PAR"
  ELF_relM_flt_par[second,first] <- "PAR"
}

# make unsymmetrical relatedness matrices 
PCC_relM_flt_par_asym <- PCC_relM_flt_par
PCC_relM_flt_par_asym[lower.tri(PCC_relM_flt_par_asym)] <- NA
PCC_unrel_pairs <- which(PCC_relM_flt_par_asym == "U", arr.ind = TRUE)

ELF_relM_flt_par_asym <- ELF_relM_flt_par
ELF_relM_flt_par_asym[lower.tri(ELF_relM_flt_par_asym)] <- NA
ELF_unrel_pairs <- which(ELF_relM_flt_par_asym == "U", arr.ind = TRUE)

# get unrelated distances
PCC_unrelated_dist <- vector(length = nrow(PCC_unrel_pairs))
for(i in 1:nrow(PCC_unrel_pairs)){
  PCC_unrelated_dist[i] <-  PCC_centroid_dist[PCC_unrel_pairs[i,][1],PCC_unrel_pairs[i,][2]]
}
ELF_unrelated_dist <- vector(length = nrow(ELF_unrel_pairs))
for(i in 1:nrow(ELF_unrel_pairs)){
  ELF_unrelated_dist[i] <-  ELF_centroid_dist[ELF_unrel_pairs[i,][1],ELF_unrel_pairs[i,][2]]
}

## estimate birth year from SVL: PCC_by_est and ELF_by_est
#PCC 
PCC_pedigree_flt_nDums <- PCC_pedigree_flt[grepl("PCC", PCC_pedigree_flt$id),] # remove dummy individuals
PCC_by_est <- as.data.frame(matrix(nrow = nrow(PCC_pedigree_flt_nDums), ncol = 3))
colnames(PCC_by_est) <- c("ID", "SVL", "InferBirthYear")
PCC_by_est$ID <- PCC_pedigree_flt_nDums$id
PCC_by_est <- merge(x = PCC_by_est, y = PCC_results[[2]], all.x = TRUE, by = "ID")

# for(i in 1:length(PCC_by_est$ID)){
#   #tryCatch(
#     #{
#       ind = PCC_by_est$ID[i]
#       no_caps <- dim(PCC_exp_metaData[[ind]])[1]
#       if(no_caps > 1){
#         # chose first capture year... 
#         focal_year <- PCC_exp_metaData[[ind]][which(PCC_exp_metaData[[ind]]$year == min(PCC_exp_metaData[[ind]]$year)),"year"]
#         SVL <- PCC_exp_metaData[[ind]][which(PCC_exp_metaData[[ind]]$year == min(PCC_exp_metaData[[ind]]$year)),"SVL"]
#         
#         # if the snake was captured more than once in the focal_year... 
#         if(length(SVL) > 1){
#           SVL <- mean(as.numeric(SVL), na.rm = T) # take the mean SVL value
#         }
#         if(length(focal_year) > 1){
#           focal_year <- focal_year[1] # and make sure length(focal_year) == 1
#         }
#         # if there isn't an SVL measurement for that year... 
#         if(is.na(as.numeric(subset(PCC_exp_metaData[[ind]], year == focal_year, select = SVL)$SVL))){
#           # chose the smallest SVL that isn't NA, and redefine SVL and focal_year
#           SVL <- PCC_exp_metaData[[ind]][which(PCC_exp_metaData[[ind]]$SVL == min(as.numeric(PCC_exp_metaData[[ind]]$SVL), na.rm = T)),"SVL"]
#           focal_year <-  PCC_exp_metaData[[ind]][which(PCC_exp_metaData[[ind]]$SVL == min(as.numeric(PCC_exp_metaData[[ind]]$SVL), na.rm = T)),"year"]
#         }
#       }else if (no_caps ==1){
#         SVL <- as.numeric(PCC_exp_metaData[[ind]]$SVL)
#         focal_year <- as.numeric(PCC_exp_metaData[[ind]]$year)
#       }
#       if(!is.na(SVL)){
#         PCC_by_est[i, "InferBirthYear"] <- infer_birth_year(ind, "barry", as.numeric(SVL), focal_year, subset(PCC_meta, ID == ind, select = Sex)$Sex)
#         PCC_by_est[i, "SVL"] <- SVL
#       }else{
#         PCC_by_est[i, "InferBirthYear"] <- NA
#         PCC_by_est[i, "SVL"] <- SVL
#       }
#     #}, 
#     #warning = function(e){
#     #  cat("warning occured at interation", i, ": ", conditionMessage(e), "\n")
#     #}
#   #)
# }  # warnings checked and can be ignored in this case

#ELF
# ELF_pedigree_flt_nD <- ELF_pedigree_flt[grepl("ELF", ELF_pedigree_flt$id),] # remove dummy individuals
ELF_by_est <- as.data.frame(matrix(nrow = nrow(ELF_pedigree_flt_nD), ncol = 3))
colnames(ELF_by_est) <- c("ID", "SVL", "InferBirthYear")
ELF_by_est$ID <- ELF_pedigree_flt_nD$id
ELF_by_est <- merge(x = ELF_by_est, y = ELF_results[[2]], all.x = TRUE, by = "ID")

# for(i in 1:length(ELF_by_est$ID)){
# #for(i in 1:141){
# 
#   #tryCatch(
#     #{
#       ind = ELF_by_est$ID[i]
#       no_caps <- dim(ELF_exp_metaData[[ind]])[1]
#       if(no_caps > 1){
#         # chose first capture year... 
#         focal_year <- ELF_exp_metaData[[ind]][which(ELF_exp_metaData[[ind]]$year == min(ELF_exp_metaData[[ind]]$year)),"year"]
#         SVL <- ELF_exp_metaData[[ind]][which(ELF_exp_metaData[[ind]]$year == min(ELF_exp_metaData[[ind]]$year)),"SVL"]
#         
#         # if the snake was captured more than once in the focal_year... 
#         if(length(SVL) > 1){
#           SVL <- mean(as.numeric(SVL), na.rm = T) # take the mean SVL value
#         }
#         if(length(focal_year) > 1){
#           focal_year <- focal_year[1] # and make sure length(focal_year) == 1
#         }
#         # if there isn't an SVL measurement for that year... 
#         if(is.na(as.numeric(SVL))){
#           # chose the smallest SVL that isn't NA, and redefine SVL and focal_year
#           SVL <- ELF_exp_metaData[[ind]][which(as.numeric(ELF_exp_metaData[[ind]]$SVL) == min(as.numeric(ELF_exp_metaData[[ind]]$SVL), na.rm = T)),"SVL"]
#           focal_year <-  ELF_exp_metaData[[ind]][which(as.numeric(ELF_exp_metaData[[ind]]$SVL) == min(as.numeric(ELF_exp_metaData[[ind]]$SVL), na.rm = T)),"year"]
#           if(length(SVL) == 0){
#             SVL <- NA
#           }
#         }
#       }else if (no_caps ==1){
#         SVL <- as.numeric(ELF_exp_metaData[[ind]]$SVL)
#         focal_year <- as.numeric(ELF_exp_metaData[[ind]]$year)
#       }
#       if(!is.na(SVL) & !is.nan(SVL)){
#         ELF_by_est[i, "InferBirthYear"] <- infer_birth_year(ind, "cass", as.numeric(SVL), focal_year, subset(ELF_meta, ID == ind, select = Sex)$Sex)
#         ELF_by_est[i, "SVL"] <- SVL
#       }else{
#         ELF_by_est[i, "InferBirthYear"] <- NA
#         ELF_by_est[i, "SVL"] <- SVL
#       }
#     #}, 
#     #warning = function(e){
#     #  cat("warning occured at interation", i, ": ", conditionMessage(e), "\n")
#     #}
#   #)
# }  # warnings checked and can be ignored in this case

# ------------------------------------------------------------------------------------------------------------

### Combine data objects
# ------------------------------------------------------------------------------------------------------------
# combine into one data object! 
# id, dam, sire, site, sex, centroid lat, centroid lon, Fgrm, heterozygosity, no. offspring, 
data <- pop_gen_stats
# data <- data[,c("id", "dam", "sire", "site", "sex", "centroid", "Fgrm", "het", "n_offspring")]

# add data from pedigrees
# are there any individuals in the genetic data not in the pedigree? 
subset(data, site == "PCC", select = id)$id[!subset(data, site == "PCC", select = id)$id %in% PCC_pedigree_flt$id] # is the genetic individual in the pedigree? 
subset(data, site == "ELF", select = id)$id[!subset(data, site == "ELF", select = id)$id %in% ELF_pedigree_flt$id] # is the genetic individual in the pedigree? 

# add additional important information to data 
mega_ped <- rbind.data.frame(PCC_pedigree_flt, ELF_pedigree_flt) # note, using FULL ELF pedigree here, not filtered
data <- merge(x = data, y = mega_ped, by = "id", all.x = TRUE)

mega_meta <- rbind.data.frame(PCC_by_est, ELF_by_est)

data <- merge(x = data, y = mega_meta, by.x = "id", by.y = "ID", all.x = TRUE)

mega_coords <- rbind.data.frame(PCC_centroids, ELF_centroids)

data <- merge(x = data, y = mega_coords, by.x = "id", by.y = "Group.1", all.x = TRUE)

#data <- merge(x = data, y = data_flt, by = "id", all.x = TRUE)
#data <- data[-816,]
mega_off <- cbind.data.frame(c(names(unlist(PCC_off)), names(unlist(ELF_off))), c(unlist(PCC_off), unlist(ELF_off)))
colnames(mega_off) <- c("id", "offspring")
data <- merge(x = data, y = mega_off, by = "id", all.x = TRUE)

data[which(is.na(data$offspring)), "offspring"] <- 0

data <- data[,c("id", "dam", "sire", "site", "Sex", "geometry", "Fgrm", "het", "offspring", "BirthYear", "BY.min", "BY.max", "InferBirthYear", "SVL")]
colnames(data) <- c("id", "dam", "sire", "site", "sex", "geometry", "Fgrm", "het", "offspring", "BirthYear", "BY.min", "BY.max", "InferBirthYear", "SVL")

# add pedigree inbreeding 
PCC_ord_ped <- pedigree::orderPed(PCC_pedigree[,1:3])
PCC_ped_inbreeding <- pedigree::calcInbreeding(PCC_pedigree[order(PCC_ord_ped),1:3])
PCC_ped_inbreeding <- cbind.data.frame(PCC_pedigree[order(PCC_ord_ped),"id"], PCC_ped_inbreeding)
colnames(PCC_ped_inbreeding) <- c("id", "ped_F")
ELF_ord_ped <- pedigree::orderPed(ELF_pedigree[,1:3])
ELF_ped_inbreeding <-pedigree::calcInbreeding(ELF_pedigree[order(ELF_ord_ped),1:3])
ELF_ped_inbreeding <- cbind.data.frame(ELF_pedigree[order(ELF_ord_ped),"id"], ELF_ped_inbreeding)
colnames(ELF_ped_inbreeding) <- c("id", "ped_F")

mega_ped_inbreeding <- rbind.data.frame(PCC_ped_inbreeding, ELF_ped_inbreeding)

data <- merge(x = data, y = mega_ped_inbreeding, by = "id", all.x = TRUE)

mega_pca<- cbind.data.frame(rownames(joint_PCA[[1]]), joint_PCA[[1]][,1:6]) # first 6 PCs

colnames(mega_pca)[1] <- "id"
data <- merge(x = data, y = mega_pca, by="id", all.x = TRUE)
colnames(data)[16:21] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

data <- data[-which(data$Fgrm > 0.8),]

data <- merge(x = data, y = export_temp, by = "id", all.x = TRUE)
sum(!data$BirthYear.x == data$BirthYear.y, na.rm = T) # can remove BirthYear.y

data <- subset(data, select = -BirthYear.y)

data$yearsContribute2Ped <- as.numeric(data$estLastYearPed) - as.numeric(data$estBirthYear)

data[which(data$sex == -9), "sex"] <- NA

### Missing data
# ------------------------------------------------------------------------------------------------------------
# sex 
dim(subset(data, is.na(sex))) # 4 individuals with no sex information

# This is unfortunate if yearsContribute2Ped > 0 
subset(data, is.na(sex) & yearsContribute2Ped > 0) # 0 

# data_update <- data
# for(i in 1:nrow(data)){
#   if(is.na(data[i,"sex"])){
#     if(sum(grepl(data[i,"id"], ELF_meta$ID)) >0){
#       if(!ELF_meta[which(ELF_meta$ID==data[i,"id"]),"Sex"] == -9){
#         data_update[i,"sex"] <- ELF_meta[which(ELF_meta$ID==data[i,"id"]),"Sex"][1]
#       }
#     }else if(sum(grepl(data[i,"id"], PCC_meta$ID)) >0){
#       if(!PCC_meta[which(PCC_meta$ID==data[i,"id"]),"Sex"] == -9){
#         data_update[i,"sex"] <- PCC_meta[which(PCC_meta$ID==data[i,"id"]),"Sex"] 
#       }
#     }
#   }
# }
# 
# dim(subset(data_update, is.na(sex) & yearsContribute2Ped > 0)) # 0! 
# dim(subset(data_update, is.na(sex))) # 4 individuals with no sex information
# 
# data <- data_update

save(data, file = "../inbreeding_models/data_for_analyses_07132024.Robj")

load("../inbreeding_models/data_for_analyses_07132024.Robj", verbose = T)
# ------------------------------------------------------------------------------------------------------------


#### random stuff 
# for each parent-offspring pair, plot Fgrm vs LLRdam/sire score
mega_ped
data
# add LLR to data
data_temp <- merge(x = data, y = mega_ped, by = "id", all.x = TRUE)

plot(data_temp$LLRdam~data_temp$Fgrm)
abline(lm(data_temp$LLRdam~data_temp$Fgrm))
plot(data_temp$LLRsire~data_temp$Fgrm)
abline(lm(data_temp$LLRsire~data_temp$Fgrm))

#### check negative LLRsire
data_temp[which(data_temp$LLRsire < 0), c("dam.x", "LLRdam", "LLRsire", "LLRpair")]

plot(log(data_temp$offspring+1)~data_temp$Fgrm)
abline(lm(log(data_temp$offspring+1)~data_temp$Fgrm))
lm(log(data_temp$offspring+1)~data_temp$Fgrm)

### Stats for the paper
# ------------------------------------------------------------------------------------------------------------
## percent of individuals with Fgrm > 0.05
(sum(subset(data, site == "PCC", select = Fgrm)$Fgrm > 0.05) / length(subset(data, site == "PCC", select = Fgrm)$Fgrm)) *100
(sum(subset(data, site == "ELF", select = Fgrm)$Fgrm > 0.05) / length(subset(data, site == "ELF", select = Fgrm)$Fgrm)) *100

# range of Fgrm at each site 
range(subset(data, site == "PCC", select = Fgrm)$Fgrm)
range(subset(data, site == "ELF", select = Fgrm)$Fgrm)


# range of offspring at each site 
range(subset(data, site == "PCC", select = offspring)$offspring)
range(subset(data, site == "ELF", select = offspring)$offspring)

# how many individuals have parents? 
percent_in_ped(subset(data, site == "PCC", select = c(id, dam, sire)), plot = FALSE)
percent_in_ped(subset(data, site == "ELF", select = c(id, dam, sire)), plot = FALSE)

# percent of individuals with no offspring? 
sum(subset(data, site == "PCC", select = offspring)$offspring == 0) / length(subset(data, site == "PCC", select = offspring)$offspring)
sum(subset(data, site == "ELF", select = offspring)$offspring == 0) / length(subset(data, site == "ELF", select = offspring)$offspring)

# by sex 
sum(subset(data, site == "PCC" & sex == 1, select = offspring)$offspring == 0) / length(subset(data, site == "PCC" & sex == 1, select = offspring)$offspring)
sum(subset(data, site == "PCC" & sex == 1, select = offspring)$offspring == 0) 
sum(subset(data, site == "PCC" & sex == 2, select = offspring)$offspring == 0) / length(subset(data, site == "PCC" & sex == 2, select = offspring)$offspring)
sum(subset(data, site == "PCC" & sex == 2, select = offspring)$offspring == 0) 

sum(subset(data, site == "ELF" & sex == 1, select = offspring)$offspring == 0) / length(subset(data, site == "ELF" & sex == 1, select = offspring)$offspring)
sum(subset(data, site == "ELF" & sex == 1, select = offspring)$offspring == 0)
sum(subset(data, site == "ELF" & sex == 2, select = offspring)$offspring == 0) / length(subset(data, site == "ELF" & sex == 2, select = offspring)$offspring)
sum(subset(data, site == "ELF" & sex == 2, select = offspring)$offspring == 0)

# how many generations? 
PCC_gen <- getGenerations(PCC_pedigree)
range(PCC_gen)
ELF_gen <- getGenerations(ELF_pedigree)
range(ELF_gen)

# pedigree confidence 
PCC.ped.conf
ELF.ped.conf

# pedigree relatedness  and kinship
range(subset(data, site == "PCC", select = ped_F))
range(subset(data, site == "ELF", select = ped_F))

ELF_pair_rel <- vector(length = length(ELF_unique_rents))
ELF_kin_rel <- vector(length = length(ELF_unique_rents))
for(i in 1:nrow(ELF_unique_rents)){
  ELF_pair_rel[i] <- ELF_relM[which(rownames(ELF_relM) == ELF_unique_rents[i,2]),which(colnames(ELF_relM) == ELF_unique_rents[i,1])]
  ELF_kin_rel[i] <- ELF_kinM[which(rownames(ELF_kinM) == ELF_unique_rents[i,2]),which(colnames(ELF_kinM) == ELF_unique_rents[i,1])]
}

cbind.data.frame(ELF_pair_rel, ELF_kin_rel)

PCC_pair_rel <- vector(length = length(PCC_unique_rents))
PCC_kin_rel <- vector(length = length(PCC_unique_rents))
for(i in 1:nrow(PCC_unique_rents)){
  PCC_pair_rel[i] <- PCC_relM[which(rownames(PCC_relM) == PCC_unique_rents[i,2]),which(colnames(PCC_relM) == PCC_unique_rents[i,1])]
  PCC_kin_rel[i] <- PCC_kinM[which(rownames(PCC_kinM) == PCC_unique_rents[i,2]),which(colnames(PCC_kinM) == PCC_unique_rents[i,1])]
}

# ELF_269
subset(ELF_pedigree, sire == "ELF_269")
ELF_269_kids <- subset(ELF_pedigree, sire == "ELF_269", select = "id")$id

subset(data, id == "ELF_473")
subset(data, id == "ELF_534")
subset(data, id == "ELF_666")
subset(data, id == "ELF_774")
subset(data, id == "ELF_828")


# ------------------------------------------------------------------------------------------------------------

### Offspring figure for talks 
# ------------------------------------------------------------------------------------------------------------

pdf(file = "../figures/offspring_wZeros.pdf", width = 8, height = 4.5)
par(mfrow=c(1,2))
par(mar = c(4, 4.5, 4, 3))
barplot(table(subset(data, site == "ELF", "offspring")$offspring), 
        xlab = "Offspring", ylab = "Number of individuals",
        main = "Cass", 
        col = colors[2], 
        axis.lty = 1, xlim = c(0,22), ylim = c(0,610), 
        cex.lab = 1.25)
barplot(table(subset(data, site == "PCC", "offspring")$offspring), 
        xlab = "Offspring", ylab = "Number of individuals",
        main = "Barry", 
        col = colors[1], 
        axis.lty = 1, xlim = c(0,10), ylim = c(0,200), 
        cex.lab = 1.25)
dev.off()

# hist(data$offspring/data$yearsContribute2Ped)
# hist(subset(data, site == "ELF", offspring)$offspring/ subset(data, site == "ELF", yearsContribute2Ped)$yearsContribute2Ped)
# hist(subset(data, site == "PCC", offspring)$offspring/ subset(data, site == "PCC", yearsContribute2Ped)$yearsContribute2Ped)
# ------------------------------------------------------------------------------------------------------------

### Figure 1: Inbreeding
# ------------------------------------------------------------------------------------------------------------

ELF_fgrm <- subset(data, site == "ELF", select = Fgrm)
rownames(ELF_fgrm) <- subset(data, site == "ELF", select = id)$id
PCC_fgrm <- subset(data, site == "PCC", select = Fgrm)
rownames(PCC_fgrm) <- subset(data, site == "PCC", select = id)$id

# pedigree set up 
# PCC
PCC_polish_ped <- sequoia::PedPolish(PCC_pedigree, DropNonSNPd=FALSE,
                                     FillParents = TRUE)
# fill in sex information
# add sex information from ELF_meta
# if "1" --> female, else, male
PCC_polish_ped$Sex <- ifelse(PCC_meta[match(PCC_polish_ped$id, PCC_meta$ID),"Sex"] == 1, yes = "female", no = "male")

# check that all genotyped individuals have sex information
ifelse(is.na(PCC_meta[match(PCC_polish_ped$id, PCC_meta$ID),"Sex"]), yes = NA, no = "FINE")

# fill in sex info for dummy individuals
PCC_polish_ped[grepl("F0", PCC_polish_ped$id), "Sex"] <- "female" # females
PCC_polish_ped[grepl("M0", PCC_polish_ped$id), "Sex"] <- "male" # males

# Fix the sex of parents, add parents that are missing from the pedigree
PCC_fix_ped <- with(PCC_polish_ped, kinship2::fixParents(id=id, dadid=sire, momid=dam,
                                                         sex=Sex))

PCC_fgrm_ped_order <- PCC_fgrm[match(PCC_fix_ped$id, rownames(PCC_fgrm)),1]

PCC_Ped_k <- with(PCC_fix_ped, kinship2::pedigree(id, dadid, momid, sex, missid=0))

genotyped <- !is.na(PCC_fgrm_ped_order)
not_geno <- is.na(PCC_fgrm_ped_order)
PCC_genotyped <- vector(length=length(PCC_Ped_k$id))
PCC_genotyped[genotyped] <- 1
PCC_genotyped[not_geno] <- 0

# ELF
# "polish" pedigree-- assign "dummy" individuals to focal inds with only 1 parent assigned
ELF_polish_ped <- sequoia::PedPolish(ELF_pedigree, DropNonSNPd=FALSE,
                                     FillParents = TRUE)
# fill in sex information
# add sex information from ELF_meta
  # if "1" --> female, else, male
ELF_polish_ped$Sex <- ifelse(ELF_meta[match(ELF_polish_ped$id, ELF_meta$ID),"Sex"] == 1, yes = "female", no = "male")

# check that all genotyped individuals have sex information
ifelse(is.na(ELF_meta[match(ELF_polish_ped$id, ELF_meta$ID),"Sex"]), yes = NA, no = "FINE")

# fill in sex info for dummy individuals
ELF_polish_ped[grepl("F0", ELF_polish_ped$id), "Sex"] <- "female" # females
ELF_polish_ped[grepl("M0", ELF_polish_ped$id), "Sex"] <- "male" # males

# Fix the sex of parents, add parents that are missing from the pedigree
ELF_fix_ped <- with(ELF_polish_ped, kinship2::fixParents(id=id, dadid=sire, momid=dam,
                                                         sex=Sex))
# match Fgrm to pedigree object
ELF_fgrm_ped_order <- ELF_fgrm[match(ELF_fix_ped$id, rownames(ELF_fgrm)),1]

# create pedigree object
ELF_Ped_k <- with(ELF_fix_ped, kinship2::pedigree(id, dadid, momid, sex, missid=0))

# for untrimmed pedigree: 
genotyped <- grep("ELF", ELF_Ped_k$id)
not_geno <- grep("ELF", ELF_Ped_k$id, invert=TRUE)
ELF_genotyped <- vector(length=length(ELF_Ped_k$id))
ELF_genotyped[genotyped] <- 1
ELF_genotyped[not_geno] <- 0

# for trimmed pedigree: 
# ELF_Ped_k_trimmed <- pedigree.shrink(ELF_Ped_k, avail = ELF_genotyped)
# ELF_fgrm_ped_order_trimmed <- ELF_fgrm[match(ELF_Ped_k_trimmed$pedObj$id, rownames(ELF_fgrm)),1]
# 
# genotyped <- grep("ELF", ELF_Ped_k_trimmed$pedObj$id)
# not_geno <- grep("ELF", ELF_Ped_k_trimmed$pedObj$id, invert=TRUE)
# ELF_genotyped <- vector(length=length(ELF_Ped_k_trimmed$pedObj$id))
# ELF_genotyped[genotyped] <- 1
# ELF_genotyped[not_geno] <- 0

# define colors scales 
# one color scale for both plots! 

ped_colors <- met.brewer("Hiroshige", n = 20) # define palette
#ped_colors_alt <- ped_colors[c(2,4,8,9,10:20)]
ped_colors_alt <- ped_colors[c(20 - c(2,4,8,9,10:20))]

all_fgrm <- subset(data, select = Fgrm)$Fgrm

# define evenly spaced colors 
x2 <- seq(length.out = length(all_fgrm), from = range(all_fgrm, na.rm = T)[1], to = range(all_fgrm, na.rm = T)[2])
c2 <- colFunc(x2, ped_colors_alt, nCols = length(x2), valRange = range(all_fgrm, na.rm =T))
#plot(x2, col = c2, xlim = c(0,1900))

# define PCC colors
PCC_cols <- colFunc(PCC_fgrm_ped_order, ped_colors_alt, nCols = length(all_fgrm), valRange = range(all_fgrm, na.rm =T))
PCC_cols[is.na(PCC_fgrm_ped_order)] <- "gray55"
# points(x = PCC_fgrm_ped_order, col = PCC_cols)

# define ELF colors
# for trimmed pedigree
ELF_cols <- colFunc(ELF_fgrm_ped_order_trimmed, ped_colors_alt, nCols = length(all_fgrm), valRange = range(all_fgrm, na.rm =T))
ELF_cols[is.na(ELF_fgrm_ped_order_trimmed)] <- "gray55"
# for untrimmed pedigree:
ELF_cols <- colFunc(ELF_fgrm_ped_order, ped_colors_alt, nCols = length(all_fgrm), valRange = range(all_fgrm, na.rm =T))
ELF_cols[is.na(ELF_fgrm_ped_order)] <- "gray55"

# points(x = 500:(499+length(ELF_fgrm_ped_order)), ELF_fgrm_ped_order, col = ELF_cols)
# abline(h=c(0,0.1, 0.2, 0.3), col = "gray", lty = 2)


#pdf(file = "../figures/fig_1_nolines_trim.pdf", width = 8, height = 10)
# png(file = "../figures/fig_1_nolines_trim.png", width = 8, height = 10, units = "in", res = 500)

# mat <- matrix(rep(c(1,1,1,1,2,3,3,3,3), 3), ncol = 3)
# mat <- matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,
#                 1,1,2,2,3,3,4,4,5,5,6,6, 
#                 rep(c(rep(7, 6), rep(8, 6)), 4), 
#                 rep(c(rep(9, 5), 10, 10, rep(11, 5)), 6)), ncol = 12)
# 
# layout(mat)

# output files separately to arrange in affinity designer
# 1: Barry inbreeding
pdf(file = "../figures/fig_1/Barry_Fgrm.pdf", width = 4, height = 4)
hist(subset(data, site == "PCC", select = Fgrm)$Fgrm, 
     xlab = "Fgrm", 
     main = "Barry")
dev.off()
# 2: Cass inbreeding
pdf(file = "../figures/fig_1/Cass_Fgrm.pdf", width = 4, height = 4)
hist(subset(data, site == "ELF", select = Fgrm)$Fgrm, 
     xlab = "Fgrm", 
     main = "Cass")
dev.off()

# Full Barry pedigree 
pdf(file = "../figures/fig_1/Barry_fullPed.pdf", width = 12, height = 6)
plot.pedigree(PCC_Ped_k, id = rep("", length(PCC_Ped_k$id)), symbolsize=1, 
              density = c(-1, 35, 65, 20), 
              affected = PCC_genotyped, col = PCC_cols, 
              width = 10)
dev.off()

# Full Cass pedigree 
pdf(file = "../figures/fig_1/Cass_fullPed.pdf", width = 12, height = 6)
plot.pedigree(ELF_Ped_k, id = rep("", length(ELF_Ped_k$id)), symbolsize=1, 
              density = c(-1, 35, 65, 20), 
              affected = ELF_genotyped, col = ELF_cols, 
              width = 10)
dev.off()

# select family groups for plotting examples 

# ELF example 1 # ------------------------------------------------------------------------------------------------------------
ELF_pair_rel[96] # pair of full siblings, and a pair of PHS or MHS? 
# focal pair: ELF_65, ELF_350
# parents: F0033, M0041
# focal offspring: ELF_921
# selection of other mates of ELF_65: ELF_392
# selection of other offspring of ELF_65: ELF_686, ELF_685, ELF_684, ELF_683, ELF_682, ELF_682, ELF_681
# selection of other mates of ELF_350: ELF_329
# selection of other offspring of ELF_350: ELF_515
subset(ELF_pedigree, dam == "ELF_65")

ped1_inds <- c("ELF_65", "ELF_350", "F0033", "M0041", "ELF_921", "ELF_392", "ELF_686", "ELF_685", "ELF_684", "ELF_683", "ELF_682", "ELF_682", "ELF_681", 
               "ELF_329", "ELF_515")
# offspring: 
# make pedigree data frame with just desired relatives 
ELF_ex_ped1 <- subset(ELF_pedigree_flt_nD, id %in% ped1_inds)

# polish
ELF_ex_ped1_polish <- sequoia::PedPolish(ELF_ex_ped1, DropNonSNPd=FALSE,
                                     FillParents = TRUE)
# add sex information 
ELF_ex_ped1_polish$Sex <- ifelse(ELF_meta[match(ELF_ex_ped1_polish$id, ELF_meta$ID),"Sex"] == 1, yes = "female", no = "male")

# fill in sex info for dummy individuals
ELF_ex_ped1_polish[grepl("F0", ELF_ex_ped1_polish$id), "Sex"] <- "female" # females
ELF_ex_ped1_polish[grepl("M0", ELF_ex_ped1_polish$id), "Sex"] <- "male" # males

# Fix the sex of parents, add parents that are missing from the pedigree
ELF_ex_ped1_fix <- with(ELF_ex_ped1_polish, kinship2::fixParents(id=id, dadid=sire, momid=dam,
                                                         sex=Sex))

# get pedigree object
ELF_ex_ped1_k <- with(ELF_ex_ped1_fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))

# trim pedigree, F0037, M0042, M0050, F0028
ELF_ex_ped1_k <- ELF_ex_ped1_k[-c(15,17,18,13)]

# get Fgrm in correct order
ELF_ex_ped1_Fgrm <- ELF_fgrm[match(ELF_ex_ped1_k$id, rownames(ELF_fgrm)),1]

# get colors 
ELF_ped1_cols <- colFunc(ELF_ex_ped1_Fgrm, ped_colors_alt, nCols = length(all_fgrm), valRange = range(all_fgrm, na.rm =T))
ELF_ped1_cols[is.na(ELF_ex_ped1_Fgrm)] <- "gray55"

# generate genotype matrix 
genotyped <- !is.na(ELF_ex_ped1_Fgrm)
not_geno <- is.na(ELF_ex_ped1_Fgrm)
ELF_ped1_genotyped <- vector(length=length(ELF_ex_ped1_k$id))
ELF_ped1_genotyped[genotyped] <- 1
ELF_ped1_genotyped[not_geno] <- 0

# plot 
# plot.pedigree(ELF_ex_ped1_k, symbolsize=1, 
#               density = c(-1, 35, 65, 20), 
#               affected = ELF_ped1_genotyped, col = ELF_ped1_cols, 
#               width = 100, 
#               align = c(1,2))
# 
pdf(file = "../figures/fig_1/Cass_ex1_ped.pdf", width = 4, height = 4)
plot.pedigree(ELF_ex_ped1_k, id = rep("", length(ELF_ex_ped1_k$id)), symbolsize=1, 
             density = c(-1, 35, 65, 20), 
             affected = ELF_ped1_genotyped, col = ELF_ped1_cols, 
             width = 100, 
             align = c(6,2))
dev.off()
# ELF example 2 #------------------------------------------------------------------------------------------------------------
# ELF_909 ELF_717 ELF_409

ELF_pair_rel[90]
ELF_unique_rents[90,] # ELF_717 ELF_409
# focal pair: ELF_717, ELF_409
# focal offspring: ELF_909
# parents of ELF_717: ELF_383, M0010
# parents of ELF_409: ELF_383, ELF_62

# other mate of ELF_717: NA
# other offspring of ELF_717: NA
# selected other mate of ELF_409: NA
# selected other offspring of ELF_409: NA

ped2_inds <- c("ELF_717", "ELF_409", "ELF_909", "ELF_383", "M0010", "ELF_62")

# offspring: 
# make pedigree data frame with just desired relatives 
ELF_ex_ped2 <- subset(ELF_pedigree_flt2, id %in% ped2_inds)

# polish
ELF_ex_ped2_polish <- sequoia::PedPolish(ELF_ex_ped2, DropNonSNPd=FALSE,
                                         FillParents = TRUE)
# add sex information 
ELF_ex_ped2_polish$Sex <- ifelse(ELF_meta[match(ELF_ex_ped2_polish$id, ELF_meta$ID),"Sex"] == 1, yes = "female", no = "male")

# fill in sex info for dummy individuals
ELF_ex_ped2_polish[grepl("F0", ELF_ex_ped2_polish$id), "Sex"] <- "female" # females
ELF_ex_ped2_polish[grepl("M0", ELF_ex_ped2_polish$id), "Sex"] <- "male" # males

# Fix the sex of parents, add parents that are missing from the pedigree
ELF_ex_ped2_fix <- with(ELF_ex_ped2_polish, kinship2::fixParents(id=id, dadid=sire, momid=dam,
                                                                 sex=Sex))

# get pedigree object
ELF_ex_ped2_k <- with(ELF_ex_ped2_fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))

# trim pedigree, M0007, F0022, M0046, F0034
ELF_ex_ped2_k <- ELF_ex_ped2_k[-c(8, 6, 10, 7)]

# get Fgrm in correct order
ELF_ex_ped2_Fgrm <- ELF_fgrm[match(ELF_ex_ped2_k$id, rownames(ELF_fgrm)),1]

# get colors 
ELF_ped2_cols <- colFunc(ELF_ex_ped2_Fgrm, ped_colors_alt, nCols = length(all_fgrm), valRange = range(all_fgrm, na.rm =T))
ELF_ped2_cols[is.na(ELF_ex_ped2_Fgrm)] <- "gray55"

# generate genotype matrix 
genotyped <- !is.na(ELF_ex_ped2_Fgrm)
not_geno <- is.na(ELF_ex_ped2_Fgrm)
ELF_ped2_genotyped <- vector(length=length(ELF_ex_ped2_k$id))
ELF_ped2_genotyped[genotyped] <- 1
ELF_ped2_genotyped[not_geno] <- 0

# plot 
plot.pedigree(ELF_ex_ped2_k, symbolsize=1, 
              density = c(-1, 35, 65, 20), 
              affected = ELF_ped2_genotyped, col = ELF_ped2_cols, 
              width = 100, 
              align = c(6,2))

pdf(file = "../figures/fig_1/Cass_ex2_ped.pdf", width = 4, height = 4)
plot.pedigree(ELF_ex_ped2_k, id = rep("", length(ELF_ex_ped2_k$id)), symbolsize=1, 
              density = c(-1, 35, 65, 20), 
              affected = ELF_ped2_genotyped, col = ELF_ped2_cols, 
              width = 100, 
              align = c(6,2))
dev.off()


# PCC example 1 #------------------------------------------------------------------------------------------------------------
PCC_pair_rel[13]
PCC_unique_rents[13,] # PCC_109 PCC_116
# focal pair: PCC_109 PCC_116
# focal offspring: PCC_150
# parents of PCC_109: F0004 M0009
# parents of PCC_116: F0004 M0009

# other mate of PCC_109: M0028
# other offspring of PCC_109: PCC_281

# other mate of PCC_109: NA
# other offspring of PCC_109: PCC_56

# selected other mate of PCC_116: NA
# selected other offspring of PCC_116: NA

ped1_inds <- c("PCC_109", "PCC_116", "PCC_150", "F0004", "M0009", "PCC_56")

# offspring: 
# make pedigree data frame with just desired relatives 
PCC_ex_ped1 <- subset(PCC_pedigree_flt, id %in% ped1_inds)

# polish
PCC_ex_ped1_polish <- sequoia::PedPolish(PCC_ex_ped1, DropNonSNPd=FALSE,
                                         FillParents = TRUE)
# add sex information 
PCC_ex_ped1_polish$Sex <- ifelse(PCC_meta[match(PCC_ex_ped1_polish$id, PCC_meta$ID),"Sex"] == 1, yes = "female", no = "male")

# fill in sex info for dummy individuals
PCC_ex_ped1_polish[grepl("F0", PCC_ex_ped1_polish$id), "Sex"] <- "female" # females
PCC_ex_ped1_polish[grepl("M0", PCC_ex_ped1_polish$id), "Sex"] <- "male" # males

# Fix the sex of parents, add parents that are missing from the pedigree
PCC_ex_ped1_fix <- with(PCC_ex_ped1_polish, kinship2::fixParents(id=id, dadid=sire, momid=dam,
                                                                 sex=Sex))

# get pedigree object
PCC_ex_ped1_k <- with(PCC_ex_ped1_fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))

# trim pedigree, M0007, F0022, M0046, F0034
PCC_ex_ped1_k <- PCC_ex_ped1_k[-c()]

# get Fgrm in correct order
PCC_ex_ped1_Fgrm <- PCC_fgrm[match(PCC_ex_ped1_k$id, rownames(PCC_fgrm)),1]

# get colors 
PCC_ped1_cols <- colFunc(PCC_ex_ped1_Fgrm, ped_colors_alt, nCols = length(all_fgrm), valRange = range(all_fgrm, na.rm =T))
PCC_ped1_cols[is.na(PCC_ex_ped1_Fgrm)] <- "gray55"

# generate genotype matrix 
genotyped <- !is.na(PCC_ex_ped1_Fgrm)
not_geno <- is.na(PCC_ex_ped1_Fgrm)
PCC_ped1_genotyped <- vector(length=length(PCC_ex_ped1_k$id))
PCC_ped1_genotyped[genotyped] <- 1
PCC_ped1_genotyped[not_geno] <- 0

# plot 
plot.pedigree(PCC_ex_ped1_k, symbolsize=1, 
              density = c(-1, 35, 65, 20), 
              affected = PCC_ped1_genotyped, col = PCC_ped1_cols, 
              width = 100)

plot(PCC_ex_ped1_k)
plot.pedigree(PCC_ex_ped1_k)

pdf(file = "../figures/fig_1/Barry_ex1_ped.pdf", width = 4, height = 4)
plot.pedigree(PCC_ex_ped1_k, id = rep("", length(PCC_ex_ped1_k$id)), symbolsize=1, 
              density = c(-1, 35, 65, 20), 
              affected = PCC_ped1_genotyped, col = PCC_ped1_cols, 
              width = 100, 
              align = c(6,2))
dev.off()


# ELF example 3 #------------------------------------------------------------------------------------------------------------
# choose highly inbred individual 
subset(data, site == "ELF" & Fgrm > 0.2) # ELF_609, 

# ELF_609, ELF_610, ELF_612
subset(ELF_pedigree, sire == "M0003")

# ELF_326 M0003
# focal pair: ELF_326 M0003
# focal offspring: ELF_607, ELF_608, ELF_609, ELF_610, ELF_611, ELF_612, ELF_646, ELF_664
# parents of ELF_326: ELF_157 M0020
# parents of M0003: ELF_326 M0038


ped3_inds <- c("ELF_326", "M0003", "ELF_607", "ELF_608", "ELF_609", "ELF_610", "ELF_611", "ELF_612", "ELF_157", "M0020", "M0038")

# offspring: 
# make pedigree data frame with just desired relatives 
ELF_ex_ped3 <- subset(ELF_pedigree_flt, id %in% ped3_inds)

# polish
ELF_ex_ped3_polish <- sequoia::PedPolish(ELF_ex_ped3, DropNonSNPd=FALSE,
                                         FillParents = TRUE)
# add sex information 
ELF_ex_ped3_polish$Sex <- ifelse(ELF_meta[match(ELF_ex_ped3_polish$id, ELF_meta$ID),"Sex"] == 1, yes = "female", no = "male")

# fill in sex info for dummy individuals
ELF_ex_ped3_polish[grepl("F0", ELF_ex_ped3_polish$id), "Sex"] <- "female" # females
ELF_ex_ped3_polish[grepl("M0", ELF_ex_ped3_polish$id), "Sex"] <- "male" # males

# Fix the sex of parents, add parents that are missing from the pedigree
ELF_ex_ped3_fix <- with(ELF_ex_ped3_polish, kinship2::fixParents(id=id, dadid=sire, momid=dam,
                                                                 sex=Sex))

# get pedigree object
ELF_ex_ped3_k <- with(ELF_ex_ped3_fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))

# trim pedigree, M0007, F0022, M0046, F0034
ELF_ex_ped3_k <- ELF_ex_ped3_k[-c()]

# get Fgrm in correct order
ELF_ex_ped3_Fgrm <- ELF_fgrm[match(ELF_ex_ped3_k$id, rownames(ELF_fgrm)),1]

# get colors 
ELF_ped3_cols <- colFunc(ELF_ex_ped3_Fgrm, ped_colors_alt, nCols = length(all_fgrm), valRange = range(all_fgrm, na.rm =T))
ELF_ped3_cols[is.na(ELF_ex_ped3_Fgrm)] <- "gray55"

# generate genotype matrix 
genotyped <- !is.na(ELF_ex_ped3_Fgrm)
not_geno <- is.na(ELF_ex_ped3_Fgrm)
ELF_ped3_genotyped <- vector(length=length(ELF_ex_ped3_k$id))
ELF_ped3_genotyped[genotyped] <- 1
ELF_ped3_genotyped[not_geno] <- 0

# plot 
plot.pedigree(ELF_ex_ped3_k, symbolsize=1, 
              density = c(-1, 35, 65, 20), 
              affected = ELF_ped3_genotyped, col = ELF_ped3_cols, 
              width = 100)

plot(ELF_ex_ped3_k)
plot.pedigree(ELF_ex_ped3_k)

pdf(file = "../figures/fig_1/Cass_ex3_ped.pdf", width = 4, height = 4)
plot.pedigree(ELF_ex_ped3_k, id = rep("", length(ELF_ex_ped3_k$id)), symbolsize=1, 
              density = c(-1, 35, 65, 20), 
              affected = ELF_ped3_genotyped, col = ELF_ped3_cols, 
              width = 100, 
              align = c(6,2))
dev.off()

###---- 

# pedigree assignment
# ------------------------------------------------------------------------------------------------------------
get_ped_percents <- function(pedigree){
  # pedigree <- subset(data, site == "PCC", select = c(id, dam, sire))
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
  print(paste0(round(parents.mia / length(geno.inds), 3), " of individuals have no parents assigned")) 
  print(paste0(round(both.parents / length(geno.inds), 3), " of individuals have both parents assigned (genotyped or dummy)"))
  print(paste0(round(one.parent / length(geno.inds), 3), " of individuals have one parent assigned (genotyped or dummy)"))
  
  percents <- sapply(c(parents.mia, one.parent, both.parents), function(x){x*100/sum(c(parents.mia, one.parent, both.parents))})
  
  return(percents)
}

ELF_ped_percents <- get_ped_percents(subset(data, site == "ELF", select = c(id, dam, sire)))
PCC_ped_percents <- get_ped_percents(subset(data, site == "PCC", select = c(id, dam, sire)))

barplot(as.matrix(cbind.data.frame(PCC_ped_percents, ELF_ped_percents)), 
        xlim = c(0,100), 
        ylim = c(0,2), 
        col = colors[c(1,3,2)], 
        border = NA, 
        horiz = T, 
        names.arg = c("Barry", "Cass"), 
        xlab = "percent assigned")


pdf(file = "../figures/fig_1/ped_assignment.pdf", width = 9.76, height = 6)
barplot(as.matrix(cbind.data.frame(PCC_ped_percents, ELF_ped_percents)), 
        xlim = c(0,100), 
        ylim = c(0,3), 
        col = c("gray", "gray50", "gray20"),
        border = NA, 
        horiz = T, 
        names.arg = c("Barry", "Cass"), 
        xlab = "percent assigned", 
        cex.lab = 2, 
        cex.names = 2, 
        cex.axis = 2)

legend("topright", legend = c("No parents assigned", "One parent assigned", "Both parents assigned"), 
       pch = 15, col = c("gray", "gray50", "gray20"), 
       bty="n",
       cex =2)

dev.off()

# number of offspring, include thoes with no offspring 
pdf(file = "../figures/fig_1/off_counts.pdf", width = 6.8, height = 3)
par(mfrow = c(1,2))
hist(subset(data, site == "PCC", select = offspring)$offspring, main = "Barry", 
     xlab = "number of assigned offspring")
hist(subset(data, site == "ELF", select = offspring)$offspring, main = "Cass", 
     xlab = "number of assigned offspring", 
     ylab = "")
dev.off()
#6
# par(oma=c(2,2,2,2)) 
# 
# par(mar=c(3,0,3,2)) 
# 
# # plot barry pedigree
# 
# plot.pedigree(PCC_Ped_k, id = rep("", length(PCC_Ped_k$id)), symbolsize=1, 
#               density = c(-1, 35, 65, 20), 
#               affected = PCC_genotyped, col = PCC_cols, 
#               width = 10)
# mtext("A", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 1.5)
# 
# par(mar=c(4,1,2,1))

#colRast <- as.raster(matrix(c2,nrow=1,ncol=length(c2)))
#plot(NULL, xlim = range(all_fgrm), ylim = c(0,10), axes = F, xlab = "")
#axis(side = 1, at = c(-0.1, 0, 0.1, 0.2, 0.3, 0.315), labels = TRUE, main = "Fgrm", xlim = range(all_fgrm))
#xleft = range(all_fgrm)[1]
#xright = range(all_fgrm)[2]
#ybottom = 0
#ytop= 10
#rasterImage(colRast,xleft,ybottom,xright,ytop,interpolate=TRUE)
#mtext("B", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 1.5)
#
#par(mar=c(3,0,3,2)) 
# plot cass pedigree
#plot.pedigree(ELF_Ped_k_trimmed$pedObj, id = rep("", length(ELF_Ped_k_trimmed$pedObj$id)), cex=2, 
#              density = c(-1, 35, 65, 20), 
#              affected = ELF_genotyped, col = ELF_cols, 
#              width = 10)
#
## plot.pedigree.mod(ELF_Ped_k, id = rep("", length(ELF_Ped_k$id)), symbolsize=2, 
#               density = c(-1, 35, 65, 20), 
#               affected = ELF_genotyped, col = ELF_cols, 
#               width = 1000)
#mtext("C", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 1.5)

#dev.off()

# figures for talks

pdf(file = "../figures/barry_only_nolines_trim.png", width = 8, height = 6)
plot.pedigree.mod(PCC_Ped_k, id = rep("", length(PCC_Ped_k$id)), symbolsize=1, 
                  density = c(-1, 35, 65, 20), 
                  affected = PCC_genotyped, col = PCC_cols, 
                  width = 10)
dev.off()

pdf(file = "../figures/scale.pdf", width = 12, height = 2)
colRast <- as.raster(matrix(c2,nrow=1,ncol=length(c2)))
plot(NULL, xlim = range(all_fgrm), ylim = c(0,10), axes = F, xlab = "", ylab = "")
axis(side = 1, at = c(-0.1, 0, 0.1, 0.2, 0.3, 0.315), labels = TRUE, main = "Fgrm", xlim = range(all_fgrm))
xleft = range(all_fgrm)[1]
xright = range(all_fgrm)[2]
ybottom = 0
ytop= 10
rasterImage(colRast,xleft,ybottom,xright,ytop,interpolate=TRUE)
dev.off()

pdf(file = "../figures/barry_only_nolines_trim_gray.png", width = 8, height = 6)

plot.pedigree.mod(PCC_Ped_k, id = rep("", length(PCC_Ped_k$id)), symbolsize=1, 
                  density = c(-1, 35, 65, 20), 
                  affected = PCC_genotyped, col = "gray40", 
                  width = 10)
dev.off()
# # plot map
# maps::map(database= "state", regions = "michigan", col = "black", resolution = 0, xlim = c(box[1], box[3]), ylim = c(box[2], box[4]), lforce = "s")
# maps::map(database= "state", regions = c("wisconsin", "illinois", "ohio", "indiana"), add = TRUE, col = "gray10", resolution = 0, xlim = c(box[1], box[3]), ylim = c(box[2], box[4]), lforce = "s")
# maps::map(database= "county", regions = "michigan,barry", add = TRUE, col = colors[1], fill = T)
# maps::map(database= "county", regions = "michigan,cass", add = TRUE, col = colors[2], fill = T)
# mtext("B", side = 3, adj = 0, outer = FALSE, cex = 2, line = 1.5)
# text(x = -85.17943, y = 43.76155, labels = "Michigan", cex = 1.5, font = 3)

# # plot barry fgrm
# par(mar=c(5,0,4,2))
# PCC_fgrm_colors <- met.brewer("OKeeffe2", n = length(PCC_fgrm$Fgrm))
# pal <- colorRampPalette(PCC_fgrm_colors)
# 
# hist(PCC_fgrm$Fgrm, 
#      breaks = 1e2,
#      border = NA,
#      xlab = "Fgrm", 
#      main = NULL, 
#      col = pal(1e2)) # fix colors! 
# mtext("C", side = 3, adj = 0, outer = FALSE, cex = 2, line = 1.5)
# 
# # plot cass fgrm 
# par(mar=c(5,0,4,2))
# ELF_fgrm_colors <- met.brewer("OKeeffe2", n = length(ELF_fgrm$Fgrm))
# pal <- colorRampPalette(ELF_fgrm_colors)
# 
# hist(ELF_fgrm$Fgrm[ELF_fgrm$Fgrm < 0.8], 
#      breaks = 1e2,
#      border = NA,
#      xlab = "Fgrm", 
#      main = NULL, 
#      col = pal(1e2)) # fix colors! 
# mtext("D", side = 3, adj = 0, outer = FALSE, cex = 2, line = 1.5)

# ------------------------------------------------------------------------------------------------------------

# non-random mating tests # ------------------------------------------------------------------------------------------------------------
# trying with genetic relatedness using KING method (calculated in get_genetic_relatedness.R)
# should also try with dummy individuals included (pedigree relatedness) -- would need to get vector of all
  # pairs of parents, not just sequenced individuals 

# load genetic relatedness data 
load(file = "/Users/meaghan/Desktop/EMR_rapture/vcf_stats/ELF_ibd.Robj")

# are parents more related than would be expected by chance? 
ELF_unique_rents

# genetic relatedness of parents vs geographic distance between parents (parent IBD figure) 

# assemble data object

ELF_unique_rents_rel <- ELF_unique_rents
ELF_unique_rents_rel$geo_dist <- ELF_parent_dist

# genetic relatedness
for(i in 1:nrow(ELF_unique_rents)){
  ids <- ELF_unique_rents[i,]
  dat1 <- ELF_ibd[which(ELF_ibd$ID1 == as.character(ids[1]) | ELF_ibd$ID2 == as.character(ids[1])),]
  # dat1[which(ELF_ibd$ID1 == as.character(ids[2]) | ELF_ibd$ID2 == as.character(ids[2])), "kinship"]
  if(length(na.omit(dat1[which(ELF_ibd$ID1 == as.character(ids[2]) | ELF_ibd$ID2 == as.character(ids[2])), "kinship"])) > 0){
    ELF_unique_rents_rel$KING[i] <- na.omit(dat1[which(ELF_ibd$ID1 == as.character(ids[2]) | ELF_ibd$ID2 == as.character(ids[2])), "kinship"])
  }else{
    ELF_unique_rents_rel$KING[i]  <- NA
  }
}

# PWP
for(i in 1:nrow(ELF_unique_rents)){
  ELF_unique_rents_rel$PWP[i] <- ELF_PWP[which(rownames(ELF_PWP) == ELF_unique_rents[i,1]), which(colnames(ELF_PWP) == ELF_unique_rents[i,2])]
}


# plot 
plot(KING~geo_dist, data = ELF_unique_rents_rel)
abline(lm(KING~geo_dist, data = ELF_unique_rents_rel))

plot(PWP~geo_dist, data = ELF_unique_rents_rel, xlab = "geo dist", ylab = "pi")
abline(lm(PWP~geo_dist, data = ELF_unique_rents_rel))
abline(v = mean(ELF_unique_rents_rel$geo_dist, na.rm = T), lty = 2)
abline(v = mean(ELF_unrelated_dist, na.rm = T), lty = 2, col = "orange")
legend("bottomright", legend = c("avg dist, parents", "avg dist, non-parents"), lty = 2, col = c("black", "orange"))


plot(PWP~KING, data = ELF_unique_rents_rel) # not sure that I had enohg loci for genetic relatedness from KING


# sample per grid square instead of randomly 


# old code
# get distribution of relatedness of observed parents (with genotypes)

# get null distribution of relatedness between pairs in a population, checking that the pairing would have been possible (individuals alive at the same time)
# use pairwise pi for now

# make function
random_mating_permut_test <- function(obs_vec, PWP_mat, data, reps, KS=TRUE, test_alt = "greater"){
  # testing 
  # obs_vec <- ELF_unique_rents_rel$KING
  # PWP_mat <- ELF_ibd_mat
  # data <- subset(data, site=="ELF" & !is.na(estBirthYear))
  # reps <- 100
  ### 
  test_outcomes <- vector(length = reps)
  
  # start permutations
  for(j in 1:reps){
    # build sample vector
    count <- 1
    sample_PWP <- vector(length = length(obs_vec))
    while(count <= length(obs_vec)){
      # choose focal individual
      focal_ind <- sample(1:nrow(data), size = 1)
      
      # choose random "mate"
      mate <- sample(1:nrow(data), size = 1)
      
      # make sure mate is opposite sex 
      if(data[mate, "sex"] != data[focal_ind, "sex"]){
        
        if(focal_ind != mate & !is.na(match(data[focal_ind,"id"], rownames(PWP_mat))) &  !is.na(match(data[mate,"id"], colnames(PWP_mat)))){
          # check if reproduction possible
          focal_repy <- as.numeric(data[focal_ind,"estBirthYear"]) + 3
          focal_dy <- as.numeric(data[focal_ind,"estLastYearPed"])
          
          mate_repy <- as.numeric(data[mate,"estBirthYear"]) + 3
          mate_dy <- as.numeric(data[mate,"estLastYearPed"])
          
          if(mate_repy < focal_dy & focal_repy < mate_dy){
            sample_PWP[count] <- PWP_mat[match(data[focal_ind,"id"], rownames(PWP_mat)), match(data[mate,"id"], colnames(PWP_mat))]
            count = count + 1
          }
        }
      }
    }
    # record KS test outcome for rep j 
    if(KS == TRUE){
      test_outcomes[j] <- ks.test(x = obs_vec, y = sample_PWP, alternative = test_alt, exact = NULL)$p.value
    }else if(KS == FALSE){
      test_outcomes[j] <- t.test(x = obs_vec, y = sample_PWP, alternative = test_alt)$p.value
    }
  }
  return(test_outcomes)
}

# get PWP of parents 

ELF_unique_rents
ELF_unique_rents_rel <- ELF_unique_rents

for(i in 1:nrow(ELF_unique_rents)){
  ELF_unique_rents_rel$PWP[i] <- ELF_PWP[which(rownames(ELF_PWP) == ELF_unique_rents[i,1]), which(colnames(ELF_PWP) == ELF_unique_rents[i,2])]
}

# no missing data 
for(i in 1:nrow(ELF_unique_rents)){
  ELF_unique_rents_rel$PWP_nMD[i] <- ELF_PWP_nMD[which(rownames(ELF_PWP) == ELF_unique_rents[i,1]), which(colnames(ELF_PWP) == ELF_unique_rents[i,2])]
}


mean(ELF_unique_rents_rel$PWP)
mean(ELF_PWP)

# PCC_KS_tests <- random_mating_permut_test(obs_vec = PCC_parent_PWP, PWP_mat = PCC_PWP, data = subset(data, site=="PCC" & !is.na(estBirthYear)), reps = 1000)
ELF_KS_tests <- random_mating_permut_test(obs_vec = ELF_unique_rents_rel$PWP, PWP_mat = ELF_PWP, data = subset(data, site=="ELF" & !is.na(estBirthYear)), reps = 100)
ELF_T_tests <- random_mating_permut_test(obs_vec = ELF_unique_rents_rel$PWP, PWP_mat = ELF_PWP, data = subset(data, site=="ELF" & !is.na(estBirthYear)), reps = 100, KS = FALSE)

sum(ELF_KS_tests <= 0.05) / 100 # not significant

# try genetic relatedness from KING 

for(i in 1:nrow(ELF_unique_rents)){
  ids <- ELF_unique_rents[i,]
  dat1 <- ELF_ibd[which(ELF_ibd$ID1 == as.character(ids[1]) | ELF_ibd$ID2 == as.character(ids[1])),]
  # dat1[which(ELF_ibd$ID1 == as.character(ids[2]) | ELF_ibd$ID2 == as.character(ids[2])), "kinship"]
  if(length(na.omit(dat1[which(ELF_ibd$ID1 == as.character(ids[2]) | ELF_ibd$ID2 == as.character(ids[2])), "kinship"])) > 0){
    ELF_unique_rents_rel$KING[i] <- na.omit(dat1[which(ELF_ibd$ID1 == as.character(ids[2]) | ELF_ibd$ID2 == as.character(ids[2])), "kinship"])
  }else{
    ELF_unique_rents_rel$KING[i]  <- NA
  }
}
# ELF_65 ELF_350 have low genetic relatedness but were marked as full sibs

which(ELF_pedigree$id == "ELF_65")
  
ELF_pedigree[which(ELF_pedigree$id == "ELF_65"),]
ELF_pedigree[which(ELF_pedigree$id == "ELF_350"),] # both parents are dummy individuals 

# covert ELF_ibd into matrix
ELF_ibd_mat <- xtabs(kinship ~ ID1 + ID2, data = ELF_ibd)

ELF_T_kinship <- random_mating_permut_test(ELF_unique_rents_rel$KING, ELF_ibd_mat, 
                                            data = subset(data, site=="ELF" & !is.na(estBirthYear)), 
                                            reps = 100, 
                                            KS = FALSE)

sum(ELF_T_kinship <= 0.05) / 100 # not significant

# try pedigree relationships 
for(i in 1:nrow(ELF_unique_rents)){
  ELF_unique_rents_rel$relationship[i] <- ELF_relM[which(rownames(ELF_relM) == ELF_unique_rents[i,1]), which(colnames(ELF_relM) == ELF_unique_rents[i,2])]
}

# or kinship matrix 
for(i in 1:nrow(ELF_unique_rents)){
  ELF_unique_rents_rel$kinship[i] <- ELF_kinM[which(rownames(ELF_kinM) == ELF_unique_rents[i,1]), which(colnames(ELF_kinM) == ELF_unique_rents[i,2])]
}
hist(ELF_unique_rents_rel$kinship)

ELF_KS_kinship <- random_mating_permut_test(ELF_unique_rents_rel$kinship, ELF_kinM, 
                                            data = subset(data, site=="ELF" & !is.na(estBirthYear)), 
                                            reps = 100, 
                                            KS = FALSE)

sum(ELF_KS_kinship <= 0.05) / 100 # not significant, if anything, parents are less related than expected in our dataset, 
  # perhaps because offspring of closely related parents are not detected 

# Parents in our dataset are not more related than expected in a panmictic population 

# Spatial distance nonparametric test 
ELF_Ttest_dist <- random_mating_permut_test(ELF_parent_dist, ELF_centroid_dist, 
                                         data = subset(data, site=="ELF" & !is.na(estBirthYear)), 
                                         reps = 100, 
                                         KS = FALSE, 
                                         test_alt = "less") # alternative hyp. here is "less" because we want to 
                                                            # ask if obs. dist less than random

sum(ELF_Ttest_dist <= 0.05) / 100 # all test have significant p-values, indicating that parents are 
# found closer together than would be expected in a panmictic population 

# is there a way to incorporate dummy individuals? 




# ------------------------------------------------------------------------------------------------------------

### Figure 2: Spatial stuff 
# ------------------------------------------------------------------------------------------------------------

# blank map colored by PC1 for both sites
PCC_pca <- data.frame(ID = rownames(PCC_pca), as.data.frame(PCC_pca))
PCC_combo <- merge(PCC_centroids, PCC_pca, by.x = "Group.1", by.y = "ID", all = FALSE)

PCC_combo$lon <- unlist(PCC_combo$geometry)[c(TRUE,FALSE)]
PCC_combo$lat <- unlist(PCC_combo$geometry)[c(FALSE,TRUE)]

# define colors
PCC_palette = "Hokusai3" # "Hiroshige" #"Hokusai3" #"OKeeffe2"
PCC_pca_colors <- rev(met.brewer(PCC_palette, n = 20))

PCC_cols_lat <- colFunc(PCC_combo$lat, PCC_pca_colors, nCols = length(PCC_combo$lat), valRange = range(PCC_combo$lat, na.rm =T))
PCC_cols_lon <- colFunc(PCC_combo$lon, PCC_pca_colors, nCols = length(PCC_combo$lon), valRange = range(PCC_combo$lon, na.rm =T))

plot(x = PCC_combo$V1, y = PCC_combo$V2, pch = 19, col = PCC_cols_lat)
plot(x = PCC_combo$V1, y = PCC_combo$V2, pch = 19, col = PCC_cols_lon)

PCC_cols_PC1 <- colFunc(PCC_combo$V1, PCC_pca_colors, nCols = length(PCC_combo$lat), valRange = range(PCC_combo$V1, na.rm =T))
plot(PCC_combo$geometry, col = PCC_cols_PC1, pch = 19)

PCC_cols_PC2 <- colFunc(PCC_combo$V2, PCC_pca_colors, nCols = length(PCC_combo$lat), valRange = range(PCC_combo$V2, na.rm =T))
plot(PCC_combo$geometry, col = PCC_cols_PC2, pch = 19)

ELF_pca <- data.frame(ID = rownames(ELF_pca), as.data.frame(ELF_pca))
ELF_combo <- merge(ELF_centroids, ELF_pca, by.x = "Group.1", by.y = "ID", all = FALSE)

ELF_combo$lon <- unlist(ELF_combo$geometry)[c(TRUE,FALSE)]
ELF_combo$lat <- unlist(ELF_combo$geometry)[c(FALSE,TRUE)]

ELF_palette =  "Hokusai3" # "Hiroshige" #"" # Hiroshige "OKeeffe2" 
ELF_pca_colors <- rev(met.brewer(ELF_palette, n = 20))

ELF_cols_lat <- colFunc(ELF_combo$lat, ELF_pca_colors, nCols = length(ELF_combo$lat), valRange = range(ELF_combo$lat, na.rm =T))
ELF_cols_lon <- colFunc(ELF_combo$lon, ELF_pca_colors, nCols = length(ELF_combo$lon), valRange = range(ELF_combo$lon, na.rm =T))

plot(x = ELF_combo$V1, y = ELF_combo$V2, pch = 19, col = ELF_cols_lat)
plot(x = ELF_combo$V1, y = ELF_combo$V2, pch = 19, col = ELF_cols_lon)

ELF_cols_PC1 <- colFunc(ELF_combo$V1, ELF_pca_colors, nCols = length(ELF_combo$lat), valRange = range(ELF_combo$V1, na.rm =T))
plot(ELF_combo$geometry, col = ELF_cols_PC1, pch = 19)

ELF_cols_PC2 <- colFunc(ELF_combo$V2, ELF_pca_colors, nCols = length(ELF_combo$lat), valRange = range(ELF_combo$V2, na.rm =T))
plot(ELF_combo$geometry, col = ELF_cols_PC2, pch = 19)

ELF_cols_PC3 <- colFunc(ELF_combo$V3, ELF_pca_colors, nCols = length(ELF_combo$lat), valRange = range(ELF_combo$V3, na.rm =T))
plot(ELF_combo$geometry, col = ELF_cols_PC3, pch = 19)

ELF_cols_PC4 <- colFunc(ELF_combo$V4, ELF_pca_colors, nCols = length(ELF_combo$lat), valRange = range(ELF_combo$V4, na.rm =T))
plot(ELF_combo$geometry, col = ELF_cols_PC4, pch = 19)

ELF_pca

dist_palette = "Hiroshige"

colors_alt <- met.brewer(dist_palette, n = 10)

pdf(file = "../figures/fig_2_altcolors2.pdf", width = 8, height = 7)
mat <- matrix(c(1, 1, 1, 3, 3, 1, 1, 1, 3, 3, 1, 1, 1, 3, 3, 2, 2, 2, 4, 4, 2, 2, 2, 4, 4, 2, 2, 2, 4, 4), ncol = 6, nrow = 5)
layout(mat)
par(mar = c(4, 4.2, 4, 4))
par(oma=c(0, 0, 0, 0) + 0.1)

# Barry, plot of PCA colored by lat
plot(x = PCC_combo$V1, y = PCC_combo$V2, pch = 19, col = PCC_cols_lat, 
     xlab = "PC 1, 4.2% of variation explained", 
     ylab = "PC2, 3.6% of variation explained", 
     #xlim = range(PCC_combo$V1), ylim = range(PCC_combo$V2),
     cex.axis = 1.2, 
     cex.lab = 1.5, 
     cex = 1.5)
     #xaxp  = c(round(range(PCC_combo$V1)[1], 2), round(range(PCC_combo$V1)[2], 2), 2), 
     #yaxp  = c(round(range(PCC_combo$V2)[1], 2), round(range(PCC_combo$V2)[2], 2), 2))
mtext("A", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 0.5)

# Cass, plot of PCA colored by lon
plot(x = ELF_combo$V1, y = ELF_combo$V2, pch = 19, col = ELF_cols_lon, 
     xlab = "PC1, 4.0% of variation explained", 
     ylab = "PC2, 3.8% of variation explained", 
     #xlim = range(ELF_combo$V1), ylim = range(ELF_combo$V2),
     cex.axis = 1.2, 
     cex.lab = 1.5, 
     cex = 1.5)
     #xaxp  = c(round(range(ELF_combo$V1)[1], 2), round(range(ELF_combo$V1)[2], 2), 2), 
     #yaxp  = c(round(range(ELF_combo$V2)[1], 2), round(range(ELF_combo$V2)[2], 2), 2))
mtext("B", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 0.5)

# violin plots of pairwise distance   

# order in rel_dists lists: unrelated_dist, dam_off_dist, sire_off_dist, parent_dist

par(mar = c(5, 4.2, 2.5, 4) + 0.2)             
# vioplot::vioplot(PCC_unrelated_dist, PCC_dam_off_dist, PCC_sire_off_dist, PCC_parent_dist, 
#         col = c("gray90", colors_alt[1], colors_alt[2], colors_alt[4]), 
#         ylab = "pairwise distance (m)", xaxt="n")

boxplot(PCC_unrelated_dist, PCC_dam_off_dist, PCC_sire_off_dist, PCC_parent_dist, 
        col = c("gray90", colors_alt[1], colors_alt[4], colors_alt[7]), 
        ylab = "pairwise distance (m)", xaxt="n", 
        cex.axis = 1.1, 
        cex.lab = 1.5)
points(x = c(1,2,3,4), y = unlist(lapply(X=list(PCC_unrelated_dist, PCC_dam_off_dist, PCC_sire_off_dist, PCC_parent_dist), FUN=mean, na.rm = T)), pch = 8)
legend("topright", legend = c("unrelated", "dam-offspring", "sire-offspring", "parents"), col = c("gray90", colors_alt[1], colors_alt[4], colors_alt[7]), 
       pch = 15, 
       cex = 1.1)
mtext("C", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 0.5)

#par(mar = c(5, 4, 0, 4) + 0.2)              # Additional space for second y-axis
# vioplot::vioplot(ELF_unrelated_dist, ELF_dam_off_dist, ELF_sire_off_dist, ELF_parent_dist,
#         col = c("gray90", colors_alt[1], colors_alt[2], colors_alt[4]), 
#         ylab = "pairwise distance (m)", xaxt="n")
boxplot(ELF_unrelated_dist, ELF_dam_off_dist, ELF_sire_off_dist, ELF_parent_dist,
                 col = c("gray90", colors_alt[1], colors_alt[4], colors_alt[7]), 
                 ylab = "pairwise distance (m)", xaxt="n", 
                 cex.axis = 1.1, 
                 cex.lab = 1.5)

points(x = c(1,2,3,4), y = unlist(lapply(X=list(ELF_unrelated_dist, ELF_dam_off_dist, ELF_sire_off_dist, ELF_parent_dist), FUN=mean, na.rm = T)), pch = 8)
mtext("D", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 0.5)

dev.off()

# scale bar for A and B plots: 

# ELF
# define evenly spaced colors 
elfx2 <- seq(length.out = length(ELF_combo$lon), from = range(ELF_combo$lon, na.rm = T)[1], to = range(ELF_combo$lon, na.rm = T)[2])
elfc2 <- colFunc(elfx2, ELF_pca_colors, nCols = length(elfx2), valRange = range(ELF_combo$lon, na.rm =T))

# ELF_cols_lon and PCC_cols_lat
colRast <- as.raster(matrix(elfc2,nrow=1,ncol=length(ELF_combo$lon)))
plot(NULL, xlim = range(ELF_combo$lon), ylim = c(0,10), axes = F, xlab = "", ylab = "")
axis(side = 1, labels = FALSE, tick = TRUE, main = "latitude", xlim = range(ELF_combo$lon))
xleft = range(ELF_combo$lon)[1]
xright = range(ELF_combo$lon)[2]
ybottom = 0
ytop= 10
rasterImage(colRast,xleft,ybottom,xright,ytop,interpolate=TRUE)
#dev.off()
# find max distance along axis for each site

# ELF
ELF_centroids
# Create a data frame with coordinates for the two points
ELF_max_coords <- data.frame(
  x = c(range(ELF_combo$lon)[1], range(ELF_combo$lon)[2]),
  y = c(4644941, 4644941)
)
# Convert the data frame to an sf object
ELF_max_points <- st_as_sf(ELF_max_coords, coords = c("x", "y"), crs = "+proj=utm +zone=16 +datum=WGS84")
st_distance(ELF_max_points) # 1086 m

# PCC
PCC_centroids
# Create a data frame with coordinates for the two points
PCC_max_coords <- data.frame(
  y = c(range(PCC_combo$lat)[1], range(PCC_combo$lat)[2]),
  x = c(-85.29962, -85.29962)
)
# Convert the data frame to an sf object
PCC_max_points <- st_as_sf(PCC_max_coords, coords = c("x", "y"), crs = "+proj=longlat +datum=WGS84")
st_distance(PCC_max_points) # 1048  m


#PCC, can use same scale bar
pccx2 <- seq(length.out = length(PCC_combo$lon), from = range(PCC_combo$lon, na.rm = T)[1], to = range(PCC_combo$lon, na.rm = T)[2])
pccc2 <- colFunc(pccx2, PCC_pca_colors, nCols = length(pccx2), valRange = range(PCC_combo$lon, na.rm =T))

pdf(file = "../figures/fig2_scale.pdf", width = 12, height = 2)
colRast <- as.raster(matrix(pccc2,nrow=1,ncol=length(PCC_combo$lon)))
plot(NULL, xlim = range(PCC_combo$lon), ylim = c(0,10), axes = F, xlab = "", ylab = "")
axis(side = 1, labels = FALSE, main = "latitude", xlim = range(PCC_combo$lon))
xleft = range(PCC_combo$lon)[1]
xright = range(PCC_combo$lon)[2]
ybottom = 0
ytop= 10
rasterImage(colRast,xleft,ybottom,xright,ytop,interpolate=TRUE)
dev.off()

# for talks

# solid color version 

pdf(file = "../figures/PCA_solid_col.pdf", width = 8, height = 7)
mat <- matrix(c(1, 1, 1, 3, 3, 1, 1, 1, 3, 3, 1, 1, 1, 3, 3, 2, 2, 2, 4, 4, 2, 2, 2, 4, 4, 2, 2, 2, 4, 4), ncol = 6, nrow = 5)
layout(mat)
par(mar = c(4, 4.2, 4, 4))
par(oma=c(0, 0, 0, 0) + 0.1)

# Barry, plot of PCA colored by lat
plot(x = PCC_combo$V1, y = PCC_combo$V2, pch = 19, col = colors[1], 
     xlab = "PC 1", 
     ylab = "PC2", 
     #xlim = range(PCC_combo$V1), ylim = range(PCC_combo$V2),
     cex.axis = 1.2, 
     cex.lab = 1.5)
#xaxp  = c(round(range(PCC_combo$V1)[1], 2), round(range(PCC_combo$V1)[2], 2), 2), 
#yaxp  = c(round(range(PCC_combo$V2)[1], 2), round(range(PCC_combo$V2)[2], 2), 2))
mtext("A", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 0.5)

# Cass, plot of PCA colored by lon
plot(x = ELF_combo$V1, y = ELF_combo$V2, pch = 19, col = colors[2], 
     xlab = "PC1", 
     ylab = "PC2", 
     #xlim = range(ELF_combo$V1), ylim = range(ELF_combo$V2),
     cex.axis = 1.2, 
     cex.lab = 1.5)
#xaxp  = c(round(range(ELF_combo$V1)[1], 2), round(range(ELF_combo$V1)[2], 2), 2), 
#yaxp  = c(round(range(ELF_combo$V2)[1], 2), round(range(ELF_combo$V2)[2], 2), 2))
mtext("B", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 0.5)

# violin plots of pairwise distance   

# order in rel_dists lists: unrelated_dist, dam_off_dist, sire_off_dist, parent_dist

par(mar = c(5, 4.2, 2.5, 4) + 0.2)             
# vioplot::vioplot(PCC_unrelated_dist, PCC_dam_off_dist, PCC_sire_off_dist, PCC_parent_dist, 
#         col = c("gray90", colors_alt[1], colors_alt[2], colors_alt[4]), 
#         ylab = "pairwise distance (m)", xaxt="n")

boxplot(PCC_unrelated_dist, PCC_dam_off_dist, PCC_sire_off_dist, PCC_parent_dist, 
        col = c("gray90", colors_alt[1], colors_alt[2], colors_alt[4]), 
        ylab = "pairwise distance (m)", xaxt="n", 
        cex.axis = 1.1, 
        cex.lab = 1.5)
points(x = c(1,2,3,4), y = unlist(lapply(X=list(PCC_unrelated_dist, PCC_dam_off_dist, PCC_sire_off_dist, PCC_parent_dist), FUN=mean, na.rm = T)), pch = 8)
legend("topright", legend = c("unrelated", "dam-offspring", "sire-offspring", "parents"), col = c("gray90", colors_alt[1], colors_alt[2], colors_alt[4]), 
       pch = 15, 
       cex = 1.1)
mtext("C", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 0.5)

#par(mar = c(5, 4, 0, 4) + 0.2)              # Additional space for second y-axis
# vioplot::vioplot(ELF_unrelated_dist, ELF_dam_off_dist, ELF_sire_off_dist, ELF_parent_dist,
#         col = c("gray90", colors_alt[1], colors_alt[2], colors_alt[4]), 
#         ylab = "pairwise distance (m)", xaxt="n")
boxplot(ELF_unrelated_dist, ELF_dam_off_dist, ELF_sire_off_dist, ELF_parent_dist,
        col = c("gray90", colors_alt[1], colors_alt[2], colors_alt[4]), 
        ylab = "pairwise distance (m)", xaxt="n", 
        cex.axis = 1.1, 
        cex.lab = 1.5)

points(x = c(1,2,3,4), y = unlist(lapply(X=list(ELF_unrelated_dist, ELF_dam_off_dist, ELF_sire_off_dist, ELF_parent_dist), FUN=mean, na.rm = T)), pch = 8)
mtext("D", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 0.5)

dev.off()

# blank map version


pdf(file = "../figures/maps_by_pc2_blue.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
# Barry, lat/lon colored by PC2
plot(PCC_combo$geometry, pch = 19, col = PCC_cols_PC2, 
     xlim = c(-85.29846, -85.2866), 
     ylim = c(42.53393, 42.54238),
     cex.axis = 1.2, 
     cex.lab = 1.5, 
     cex = 1.5)
mtext("longitude", side = 1, line = 0.5, cex = 1.2)
mtext("latitude", side = 2, line = 0.5, cex = 1.2)
box()
# Cass, lat/lon colored by PC2
plot(ELF_combo$geometry, pch = 19, col = ELF_cols_PC2, 
     xlim = c(581913, 582936.4), 
     ylim = c(4644566, 4645642),
     cex.axis = 1.2, 
     cex.lab = 1.5, 
     cex = 1.5)
mtext("longitude", side = 1, line = 0.5, cex = 1.2)
mtext("latitude", side = 2, line = 0.5, cex = 1.2)
box()
dev.off()


# box plots
palette_box = "Hokusai3" #"OKeeffe2"
box_colors <- met.brewer(palette_box, n = 3)

pdf(file = "../figures/rel_dists_talks_blue.pdf", width = 12, height = 6)
par(mar = c(5, 4.5, 2.5, 4) + 0.2)
par(mfrow = c(1,2))
# vioplot::vioplot(PCC_unrelated_dist, PCC_dam_off_dist, PCC_sire_off_dist, PCC_parent_dist, 
#         col = c("gray90", colors_alt[1], colors_alt[2], colors_alt[4]), 
#         ylab = "pairwise distance (m)", xaxt="n")

boxplot(PCC_unrelated_dist, PCC_dam_off_dist, PCC_sire_off_dist, PCC_parent_dist, 
        col = c("gray90", box_colors), 
        ylab = "pairwise geographic distance (m)", xaxt="n", 
        cex.axis = 1.1, 
        cex.lab = 1.5)
points(x = c(1,2,3,4), y = unlist(lapply(X=list(PCC_unrelated_dist, PCC_dam_off_dist, PCC_sire_off_dist, PCC_parent_dist), FUN=mean, na.rm = T)), pch = 8)
legend("topright", legend = c("unrelated", "dam-offspring", "sire-offspring", "parents"), col = c("gray90", box_colors), 
       pch = 15, 
       cex = 1.1)

boxplot(ELF_unrelated_dist, ELF_dam_off_dist, ELF_sire_off_dist, ELF_parent_dist,
        col = c("gray90", box_colors), 
        ylab = "pairwise geographic distance (m)", xaxt="n", 
        cex.axis = 1.1, 
        cex.lab = 1.5)

points(x = c(1,2,3,4), y = unlist(lapply(X=list(ELF_unrelated_dist, ELF_dam_off_dist, ELF_sire_off_dist, ELF_parent_dist), FUN=mean, na.rm = T)), pch = 8)
dev.off()

# Kolmogorov-Smirnov Tests
# use bonferroni corrections! 6 tests total

# Are related distributed different from unrelated distribution? 
# unrelated_dist vs dam_off_dist
PCC_unrelated_dam <- ks.test(x = unique(PCC_unrelated_dist), y = unique(PCC_dam_off_dist), alternative = "two.sided", exact = TRUE)
PCC_unrelated_dam$p.value*6

# unrelated_dist vs sire_off_dist 
PCC_unrelated_sire <- ks.test(x = unique(PCC_unrelated_dist), y = unique(PCC_sire_off_dist), alternative = "two.sided", exact = TRUE)
PCC_unrelated_sire$p.value*6

# unrelated_dist vs parent_dist 
PCC_unrelated_parent <- ks.test(x = unique(PCC_unrelated_dist), y = unique(PCC_parent_dist), alternative = "two.sided", exact = TRUE)
PCC_unrelated_parent$p.value*6

# unrelated_dist vs dam_off_dist
# use bonferroni corrections! 12 tests total
ELF_unrelated_dam <- ks.test(x = unique(ELF_unrelated_dist), y = unique(ELF_dam_off_dist), alternative = "two.sided", exact = TRUE)
ELF_unrelated_dam$p.value*6

# unrelated_dist vs sire_off_dist
ELF_unrelated_sire <- ks.test(x = unique(ELF_unrelated_dist), y = unique(ELF_sire_off_dist), alternative = "two.sided", exact = TRUE)
ELF_unrelated_sire$p.value*6

# unrelated_dist vs parent_dist
ELF_unrelated_parent <- ks.test(x = unique(ELF_unrelated_dist), y = unique(ELF_parent_dist), alternative = "two.sided", exact = TRUE)
ELF_unrelated_parent$p.value*6


# dam_off_dist vs sire_off_dist
PCC_dam_sire <- ks.test(x = unique(PCC_unrelated_dist), y = unique(PCC_sire_off_dist), alternative = "two.sided", exact = TRUE)
PCC_dam_sire$p.value

# dam_off_dist vs parent_dist
PCC_dam_parent <- ks.test(x = unique(PCC_dam_off_dist), y = unique(PCC_parent_dist), alternative = "two.sided", exact = TRUE)
PCC_dam_parent$p.value

# sire_off_dist vs parent_dist
PCC_sire_parent <- ks.test(x = unique(PCC_sire_off_dist), y = unique(PCC_parent_dist), alternative = "two.sided", exact = TRUE)
PCC_sire_parent$p.value

# dam_off_dist vs sire_off_dist
ELF_dam_sire <- ks.test(x = unique(ELF_unrelated_dist), y = unique(ELF_sire_off_dist), alternative = "two.sided", exact = TRUE)
ELF_dam_sire$p.value

# dam_off_dist vs parent_dist
ELF_dam_parent <- ks.test(x = unique(ELF_dam_off_dist), y = unique(ELF_parent_dist), alternative = "two.sided", exact = TRUE)
ELF_dam_parent$p.value

# sire_off_dist vs parent_dist
ELF_sire_parent <- ks.test(x = unique(ELF_sire_off_dist), y = unique(ELF_parent_dist), alternative = "two.sided", exact = TRUE)
ELF_sire_parent$p.value

# logistic regression: does distance predict relatedness (parent/offspring)? 
# PCC_rel_dists  # unrelated_dist, dam_off_dist, sire_off_dist, parent_dist

# unrelated <- cbind.data.frame(rep(0, length(PCC_unrelated_dist)), PCC_unrelated_dist)
# colnames(unrelated) <- c("related", "dist")
# dam_dist <- cbind.data.frame(rep(1, length(PCC_dam_off_dist)), PCC_dam_off_dist)
# colnames(dam_dist) <- c("related", "dist")
# sire_dist <- cbind.data.frame(rep(1, length(PCC_sire_off_dist)), PCC_sire_off_dist)
# colnames(sire_dist) <- c("related", "dist")
# 
# dists <- rbind.data.frame(unrelated, dam_dist, sire_dist)
# 
# plot(related~dist, data = dists)
# with(dists, lines(lowess(related~dist)))
# 
# rel_mod <- glmmTMB(related~dist, data = dists, family = "binomial")
# summary(rel_mod)
# 
# plot(allEffects(rel_mod))
# 
# # -fixef(rel_mod)$cond[1]/fixef(rel_mod)$cond[2]
# 
# fixef(rel_mod)$cond[2]
# 
# newdata = with(dists, data.frame(x = seq(min(dist), max(dist), len = 100)))
# Xmat = model.matrix(~x, data = newdata)
# coefs = fixef(rel_mod)$cond
# fit = as.vector(coefs %*% t(Xmat))
# se = sqrt(diag(Xmat %*% vcov(rel_mod)$cond %*% t(Xmat)))
# q = qnorm(0.975)
# newdata = cbind(newdata, fit = binomial()$linkinv(fit), lower = binomial()$linkinv(fit -
#                                                                                      q * se), upper = binomial()$linkinv(fit + q * se))
# ggplot(data = newdata, aes(y = fit, x = x)) + geom_point(data = dists,
#                                                          aes(y = related, x = dist)) + geom_ribbon(aes(ymin = lower, ymax = upper),
#                                                                                    fill = "blue", alpha = 0.3) + geom_line() + theme_classic()

# palette(PCC_map_colors)
# PCC_box <- c(left = -85.29791, bottom = 42.53295, right = -85.28568, top = 42.54254)
# plot(PCC_combo$lat~PCC_combo$lon, col = as.factor(PCC_combo$V1), pch = 19, main = "", 
#      axes = TRUE, xaxt = "n", yaxt = "n", 
#      xlim = PCC_box[c(1,3)], ylim = PCC_box[c(2,4)])
# mtext("A", side = 3, adj = 0, outer = FALSE, cex = 2, line = 0.5)
# mtext("Barry", side = 3, adj = 0.5, outer = FALSE, cex = 2, line = 0.5)
# 
# palette(ELF_map_colors)
# plot(ELF_combo$lat~ELF_combo$lon, col = as.factor(ELF_combo$V1), pch = 19, main = "", 
#      axes = TRUE, xaxt = "n", yaxt = "n")
# # mtext("Cass", side = 3, font = 2, line = 1, cex = 1.5)
# mtext("B", side = 3, adj = 0, outer = FALSE, cex = 2, line = 0.5)
# mtext("Cass", side = 3, adj = 0.5, outer = FALSE, cex = 2, line = 0.5)
# 
# hist(PCC_rel_dists[[1]], breaks = 50, xlim = c(0,1200), xlab = "pairwise distance [m]", main = "", col = "gray90", 
#      ylab = "Frequency")
# par(new = TRUE)                             # Add new plot
# hist(PCC_rel_dists[[2]], breaks = 20, xlim = c(0,1200), ylim = c(0,30), col = alpha(colors[4], alpha = 0.75), axes = FALSE, xlab = "", ylab = "", main = "")
# hist(PCC_rel_dists[[3]], breaks = 10, add = T, col = alpha(colors[5], alpha = 0.75), axes = FALSE, xlab = "", ylab = "", main = "")
# hist(PCC_rel_dists[[4]], breaks = 20, add = T, col = alpha(colors[1], alpha = 0.75), axes = FALSE, xlab = "", ylab = "", main = "")
# axis(side = 4, at = pretty(c(0,30)))      # Add second axis
# #corners = par("usr")
# # mtext("Related Frequency", side = 4, line = 3, cex = 0.75)             # Add second axis label

# hist(ELF_rel_dists[[1]], breaks = 50, xlim = c(0,1200), xlab = "pairwise distance [m]", main = "", col = "gray90", 
#      ylab = "")
# par(new = TRUE)                             # Add new plot
# hist(ELF_rel_dists[[2]], breaks = 20, xlim = c(0,1200), ylim = c(0,100), col = alpha(colors[4], alpha = 0.75), axes = FALSE, xlab = "", ylab = "", main = "")
# hist(ELF_rel_dists[[3]], add = T, breaks = 20, col = alpha(colors[5], alpha = 0.75), axes = FALSE, xlab = "", ylab = "", main = "")
# hist(ELF_rel_dists[[4]], breaks = 20, add = T, col = alpha(colors[1], alpha = 0.75), axes = FALSE, xlab = "", ylab = "", main = "")
# axis(side = 4, at = pretty(c(0,100)))      # Add second axis
# #mtext("Related Frequency", side = 4, line = 3)             # Add second axis label


# ------------------------------------------------------------------------------------------------------------

### Figure 3: Inbreeding depression 
# ------------------------------------------------------------------------------------------------------------

palette = "Hiroshige" #"OKeeffe2"
colors <- met.brewer(palette, n = 10)

# cols <- unlist(lapply(RO_mod_data$site, function(x){if(x=="PCC"){return(colors[2])}else{return(colors[6])}})) # cols by site

# plot probability of excess zero 
load(file = "../inbreeding_models/zi_Fgrm_plot_data.Robj", verbose = TRUE)
emm_results_zi_df <- zi_plot[[1]]
offspring_binary <- zi_plot[[2]]
data_flt <- zi_plot [[3]]

plot(response~orig_Fgrm, data = emm_results_zi_df, type = "l", col = "black", 
     xlab = "Fgrm", ylab = "Probability of excess zero", 
     cex.axis = 1.2, 
     cex.lab = 1.75, 
     xlim = range(emm_results_zi_df$orig_Fgrm), 
     ylim = c(0,1))
points(offspring_binary~data_flt$Fgrm, pch = 19, col = alpha(colors[2], alpha = 0.1))
polygon(x=c(emm_results_zi_df$orig_Fgrm, rev(emm_results_zi_df$orig_Fgrm)), 
        y = c(emm_results_zi_df$asymp.UCL, rev(emm_results_zi_df$asymp.LCL)), 
        col = alpha("gray", alpha = 0.5), border = NA)

# load data to plot survival model 
###Apparent survival effects of inbreeding (Frgm) on ELF and PCCI adult controlling for SVL (52) -0.0946 to 0.314829
dat2 <- read.csv("../beta_coefficients3_15June2024/survival_Fgrm_ELF&PCCI_52svl_15June24.csv")
surv_est <- read.csv("../beta_coefficients3_15June2024/Individual_estimates_SVL52.csv")

plot(NULL, xlim = c(-0.1, 0.35), ylim = c(0, 1), 
     xlab = "Fgrm", ylab = "Apparent survival")
polygon(x=c(subset(dat2, Site == "ELF")$Fgrm, rev(subset(dat2, Site == "ELF")$Fgrm)), y = c(subset(dat2, Site == "ELF")$UCI, rev(subset(dat2, Site == "ELF")$LCI)), 
        col = alpha(colors[7], alpha = 0.5), border = NA)
polygon(x=c(subset(dat2, Site == "PCCI")$Fgrm, rev(subset(dat2, Site == "PCCI")$Fgrm)), y = c(subset(dat2, Site == "PCCI")$UCI, rev(subset(dat2, Site == "PCCI")$LCI)), 
        col = alpha("gray", alpha = 0.5), border = NA)
points(Phi~Fgrm, data = subset(surv_est, Site_name == "PCCI"), pch = 21, bg = alpha("gray", 0.5), col = alpha("gray", 1))
points(Phi~Fgrm, data = subset(surv_est, Site_name == "ELF"), pch = 21, bg = alpha(colors[7], 0.5), col = alpha(colors[7], 1))

lines(Phi~Fgrm, data = subset(dat2, Site == "ELF"), lty = 2)
lines(Phi~Fgrm, data = subset(dat2, Site == "PCCI"))

# plot together
tiff(file="../figures/fig_3.tiff", width=12, height=6, units="in", res=400)
par(mar = c(5, 4.2, 2.5, 4) + 0.2)             
par(mfrow=c(1,2))
# survival 
plot(NULL, xlim = c(-0.1, 0.35), ylim = c(0, 0.97), 
     xlab = "Fgrm", ylab = "Apparent survival", 
     cex.axis = 1.2,
     cex.lab = 1.75)
polygon(x=c(subset(dat2, Site == "ELF")$Fgrm, rev(subset(dat2, Site == "ELF")$Fgrm)), y = c(subset(dat2, Site == "ELF")$UCI, rev(subset(dat2, Site == "ELF")$LCI)), 
        col = alpha(colors[7], alpha = 0.5), border = NA)
polygon(x=c(subset(dat2, Site == "PCCI")$Fgrm, rev(subset(dat2, Site == "PCCI")$Fgrm)), y = c(subset(dat2, Site == "PCCI")$UCI, rev(subset(dat2, Site == "PCCI")$LCI)), 
        col = alpha("gray", alpha = 0.5), border = NA)
# points(Phi~Fgrm, data = subset(surv_est, Site_name == "PCCI"), pch = 21, bg = alpha("gray", 0.5), col = alpha("gray", 1))
# points(Phi~Fgrm, data = subset(surv_est, Site_name == "ELF"), pch = 21, bg = alpha(colors[7], 0.5), col = alpha(colors[7], 1))
axis(side = 1, tck= 0.03, at = subset(surv_est, Site_name == "PCCI")$Fgrm, labels = FALSE, col = "black")
axis(side = 3, tck = 0.03, at = subset(surv_est, Site_name == "ELF")$Fgrm, labels = FALSE, col = colors[7])

lines(Phi~Fgrm, data = subset(dat2, Site == "ELF"))
lines(Phi~Fgrm, data = subset(dat2, Site == "PCCI"), lty = 2)
mtext("A", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 1.0)
legend(x= -0.118, y = 0.2, legend = c("Barry", "Cass"), lty = c(2,1), col = "black", fill = c("gray", colors[7]))
box()
# reproductive success 
plot(response~orig_Fgrm, data = emm_results_zi_df, type = "l", col = "black", 
     xlab = "Fgrm", ylab = "Probability of excess zero", 
     cex.axis = 1.2, 
     cex.lab = 1.75, 
     xlim = range(emm_results_zi_df$orig_Fgrm), 
     ylim = c(0,1))
points(offspring_binary~data_flt$Fgrm, pch = 21, bg = alpha(colors[2], alpha =0.25), col = alpha(colors[2], alpha = 1))
polygon(x=c(emm_results_zi_df$orig_Fgrm, rev(emm_results_zi_df$orig_Fgrm)), 
        y = c(emm_results_zi_df$asymp.UCL, rev(emm_results_zi_df$asymp.LCL)), 
        col = alpha("gray", alpha = 0.5), border = NA)
mtext("B", side = 3, adj = 0, outer = FALSE, cex = 1.5, line = 1.0)

dev.off()

# compare "inbred" and "not inbred" individuals
# Individuals above the 80th quantile of inbreeding in their population are 
# X times less likely to have offspring and are Y % lower apparent annual survival. 
head(data)
# get Fgrm value that represents the 95th quantile of inbreeding for each site
PCC_quant <- quantile(subset(data, site == "PCC", select = Fgrm)$Fgrm, probs = 0.95) # 0.05989676
ELF_quant <- quantile(subset(data, site == "ELF", select = Fgrm)$Fgrm, probs = 0.95) # 0.0879195
all_quant <- quantile(subset(data, select = Fgrm)$Fgrm, probs = 0.95) # 0.07974094

# survival 
mean(subset(surv_est, Site_name == "PCCI")[which(subset(surv_est, Site_name == "PCCI")$Fgrm >= PCC_quant), "Phi"])
# Average survival of individuals above the 95th quantile of inbreeding is 71% at average adult size 

mean(subset(surv_est, Site_name == "PCCI")[which(subset(surv_est, Site_name == "PCCI")$Fgrm < PCC_quant), "Phi"])
# Average survival of individuals below the 95th quantile of inbreeding is 80.4% at average adult size

mean(subset(surv_est, Site_name == "ELF")[which(subset(surv_est, Site_name == "ELF")$Fgrm >= ELF_quant), "Phi"])
# Average survival of individuals above the 95th quantile of inbreeding is 54.8% at average adult size 

mean(subset(surv_est, Site_name == "ELF")[which(subset(surv_est, Site_name == "ELF")$Fgrm < ELF_quant), "Phi"])
# Average survival of individuals below the 95th quantile of inbreeding is 68.7% at average adult size

100* (mean(subset(surv_est, Site_name == "ELF")[which(subset(surv_est, Site_name == "ELF")$Fgrm < ELF_quant), "Phi"])
 - mean(subset(surv_est, Site_name == "ELF")[which(subset(surv_est, Site_name == "ELF")$Fgrm >= ELF_quant), "Phi"]))

100* (mean(subset(surv_est, Site_name == "PCCI")[which(subset(surv_est, Site_name == "PCCI")$Fgrm < PCC_quant), "Phi"])
      - mean(subset(surv_est, Site_name == "PCCI")[which(subset(surv_est, Site_name == "PCCI")$Fgrm >= PCC_quant), "Phi"]))

(13.79791  + 9.462528) / 2 # 11.63022
# Predicted apparent annual survival of the most inbred individuals in the population decreases by 13.8%/9.46% compared 
# to the most inbred (upper 95th quantile)


# what about compared to Fgrm = 0? 


# reproduction 
# use emm_results_zi_obs_df
100* (mean(subset(surv_est, Site_name == "PCCI")[which(subset(surv_est, Site_name == "PCCI")$Fgrm < PCC_quant), "Phi"])
      - mean(subset(surv_est, Site_name == "PCCI")[which(subset(surv_est, Site_name == "PCCI")$Fgrm >= PCC_quant), "Phi"]))

mean(emm_results_zi_obs_df[which(emm_results_zi_obs_df$orig_Fgrm >= all_quant), "response"])
# 90.7%  probability of excess zero 

mean(emm_results_zi_obs_df[which(emm_results_zi_obs_df$orig_Fgrm < all_quant), "response"])
# 77.2% probability of excess zero 

100 * (mean(emm_results_zi_obs_df[which(emm_results_zi_obs_df$orig_Fgrm >= all_quant), "response"]) - mean(emm_results_zi_obs_df[which(emm_results_zi_obs_df$orig_Fgrm < all_quant), "response"]))


# the probability of having any offspring is 13.46132% less in the upper 95th quantile of inbreeding compared to less inbred individuals 

# what about compared to Fgrm = 0? 


# talk version 
# survival 
# plot together
tiff(file="../figures/survival.tiff", width=6, height=6, units="in", res=400)
par(mar = c(5, 4.2, 2.5, 4) + 0.2)             
plot(NULL, xlim = c(-0.1, 0.35), ylim = c(0, 0.97), 
     xlab = "Fgrm", ylab = "Apparent survival", 
     cex.axis = 1.2,
     cex.lab = 1.75)
polygon(x=c(subset(dat2, Site == "ELF")$Fgrm, rev(subset(dat2, Site == "ELF")$Fgrm)), y = c(subset(dat2, Site == "ELF")$UCI, rev(subset(dat2, Site == "ELF")$LCI)), 
        col = alpha(colors[7], alpha = 0.5), border = NA)
polygon(x=c(subset(dat2, Site == "PCCI")$Fgrm, rev(subset(dat2, Site == "PCCI")$Fgrm)), y = c(subset(dat2, Site == "PCCI")$UCI, rev(subset(dat2, Site == "PCCI")$LCI)), 
        col = alpha("gray", alpha = 0.5), border = NA)
lines(Phi~Fgrm, data = subset(dat2, Site == "ELF"))
lines(Phi~Fgrm, data = subset(dat2, Site == "PCCI"), lty = 2)
legend("bottomleft", legend = c("Barry", "Cass"), lty = c(2,1), col = "black", fill = c("gray", colors[7]))

dev.off()


# reproductive success

tiff(file="../figures/zi_reproductive_output_lines.tiff", width=6, height=6, units="in", res=400)
par(mar = c(5, 4.2, 2.5, 4) + 0.2)             
plot(response~orig_Fgrm, data = emm_results_zi_df, type = "l", col = "black", 
     xlab = "Fgrm", ylab = "Probability of pedigree offspring", 
     cex.axis = 1.2, 
     cex.lab = 1.75, 
     xlim = range(emm_results_zi_df$orig_Fgrm), 
     ylim = c(0,1))
points(offspring_binary~data_flt$Fgrm, pch = 19, col = alpha(colors[2], alpha = 0.1))
polygon(x=c(emm_results_zi_df$orig_Fgrm, rev(emm_results_zi_df$orig_Fgrm)), 
        y = c(emm_results_zi_df$asymp.UCL, rev(emm_results_zi_df$asymp.LCL)), 
        col = alpha("gray", alpha = 0.5), border = NA)
abline(v = 0.00, lty = 2, col = "black")
abline(v = 0.05, lty = 2, col = "black")
dev.off()

# ggplot versions
# survival_plot = ggplot(dat2, aes(x=Fgrm, y=Phi, ymax=UCI, ymin=LCI, fill=Site, linetype=Site)) +
#   geom_ribbon(alpha=0.5) +
#   geom_line() +  
#   ylab("Apparent survival") + xlab("Fgrm") +
#   scale_x_continuous(breaks=seq(-0.10, 0.3, by=0.05), limits=c(-0.0946,0.3149))+
#   scale_y_continuous(breaks=seq(0.0, 1, by=0.1), limits=c(0.0, 1))+
#   theme_classic() + theme(text = element_text(size=14), axis.text = element_text(colour='black'),
#                           #legend.position.inside = c(0.9, 0.9),
#                           axis.line.x = element_line(colour = 'black', size=0, linetype='solid'),
#                           axis.line.y = element_line(colour = 'black', size=0, linetype='solid'), 
#                           axis.text.x=element_text(angle=0,hjust=0.5, family = "Helvetica"), 
#                           axis.text.y=element_text(family = "Helvetica"), 
#                           axis.title.x = element_text(family = "Helvetica"), 
#                           axis.title.y = element_text(family = "Helvetica"), 
#                           panel.border = element_rect(colour = "black", fill=NA, linewidth=0.75)) + 
#   theme(legend.position="none") +
#   scale_fill_manual(values = c("turquoise3","gray")) 
# 
# survival_plot
# 
# # plot 
# Fgrm_plot = ggplot(data = Fgrm_effects, aes(x=Fgrm, y=fit)) +
#   geom_point(data = RO_mod_data, aes(x=Fgrm, y=offspring, fill = sex), shape = 21, stroke = 0, color = "black", size = 3, alpha = 0.5) +
#   geom_ribbon(aes(ymax=lower, ymin=upper), alpha=0.25) +
#   geom_line() + 
#   ylab("Number of offspring") + xlab("Fgrm") + 
#   scale_x_continuous(breaks=seq(-0.10, 0.3, by=0.05), limits=c(-0.0946,0.315))+
#   scale_y_continuous(breaks=seq(0.0, 20, by=5), limits=c(0.0, 20)) +
#   theme_classic() + theme(text = element_text(size=14),axis.text = element_text(colour='black'),
#                           #legend.position.inside = c(0.9, 0.9),
#                           axis.line.x = element_line(colour = 'black', size=0, linetype='solid'),
#                           axis.line.y = element_line(colour = 'black', size=0, linetype='solid'), 
#                           axis.text.x=element_text(angle=0,hjust=0.5, family = "Helvetica"), 
#                           axis.text.y=element_text(family = "Helvetica"), 
#                           axis.title.x = element_text(family = "Helvetica"), 
#                           axis.title.y = element_text(family = "Helvetica"), 
#                           panel.border = element_rect(colour = "black", fill=NA, linewidth=0.75)) + 
#   theme(legend.position="none") +
#   scale_fill_manual(values = c("1" = colors[2] ,"2" = colors[4]))
# 
# Fgrm_plot
# 


# pdf(file = "../figures/fig_2.pdf", width = 8, height = 4.5)
# par(mfrow=c(1,2))
# plot(data_flt$offspring~ data_flt$Fgrm, pch = 19, col = alpha("gray19", alpha = 0.5), 
#      xlab = "Fgrm", 
#      ylab = "number of offspring")
# abline(real_mod, col = site_palette[4], lwd = 1.5)
# mtext("A", side = 3, adj = 0, outer = FALSE, cex = 2, line = 0.5)
# 
# hist(permut_slopes, main = NULL, xlab = "slope of offspring~inbreeding", 
#      xlim = c(-10, 15), breaks = 50, col = "gray80", border = "gray80")
# mtext("10,000 permutations")
# mtext("B", side = 3, adj = 0, outer = FALSE, cex = 2, line = 0.5)
# 
# abline(v = quants, col = "black", lty = 2)
# abline(v = real_coeff, col = colors[4], lwd =1.5)
# dev.off()

# ------------------------------------------------------------------------------------------------------------

### Figure 3/4: Survival and inbreeding models 
# ------------------------------------------------------------------------------------------------------------
data_flt <- data[-which(is.na(data$sex)),]
data_flt$sex <- as.factor(data_flt$sex) 
data_flt$site <- as.factor(data_flt$site)
data_flt <- data_flt[-which(is.na(data_flt$yearsContribute2Ped)),]
data_flt <- data_flt[-which(data_flt$yearsContribute2Ped == 0),]

# subset data
data_flt <- subset(data_flt, select = c(offspring, Fgrm, sex, site, PC1, PC2, PC3, PC4, PC5, PC6, yearsContribute2Ped))
# remove outliers in yearsContribute2Ped
data_flt <- subset(data_flt, subset = yearsContribute2Ped <= 16)

save(data_flt, file = "../inbreeding_models/data_for_model_troubleshooting_06272024.RObj")
# visualize predictors
GGally::ggpairs(data_flt)
plot(offspring~Fgrm, data = data_flt, pch = 19, col = alpha("black", alpha = 0.25))

# updating models to use years contributing offset

cprior <- data.frame(prior = rep("normal(0,3)", 2),
                     class = rep("beta", 2),
                     coef = c("(Intercept)", ""))

mod1 <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6 + offset(log(yearsContribute2Ped)),
                ziformula = ~., 
                data = data_flt, 
                family = nbinom2)
summary(mod1)
plot(DHARMa::simulateResiduals(mod1))
diagnose(mod1) # Unusually large coefficients, Unusually large Z-statistics
DHARMa::testDispersion(mod1)


################################################################################
data_old <- subset(data, InferBirthYear < 2012) # 447 individuals
dim(subset(data_old, site == "PCC")) # 108 PCC
dim(subset(data_old, site == "ELF")) # 339 ELF
GGally::ggpairs(data_old)
plot(offspring~Fgrm, data = data_old, pch = 19, col = alpha("black", alpha = 0.25))

save(data_old, file = "data_for_cstat.Robj")
################################################################################
# silly age approximation 
data_test <- data_flt
data_test <- subset(data_test, !is.na(InferBirthYear))
data_test$dyear <- data_test$InferBirthYear + 10
for(i in 1:nrow(data_test)){
  if(data_test$dyear[i] > 2020){
    data_test$dyear[i] <- 2020
  }
}
data_test$age <- data_test$dyear - data_test$InferBirthYear

test_mod1 <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                ziformula = ~.,
                offset = age,
                data = data_test, 
                family = nbinom2)

test_mod2 <- glmmTMB(offspring ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                ziformula = ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                offset = age,
                data = data_test, 
                family = nbinom2)

test_mod3 <- glmmTMB(offspring ~ sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                ziformula = ~ sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                offset = age,
                data = data_test, 
                family = nbinom2)


test_mod_list <- list(test_mod1, test_mod2, test_mod3)
AICcmodavg::aictab(test_mod_list)

################################################################################

# Poisson vs conway maxwell poisson, likelihood ratio test using truncated data 

mod1alt <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                   ziformula = ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                   data = data_old, 
                   family = compois)
summary(mod1alt)
plot(DHARMa::simulateResiduals(mod1alt))
diagnose(mod1alt)
DHARMa::testDispersion(mod1alt)

mod1 <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                ziformula = ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                data = data_old, 
                family = poisson)
summary(mod1)
plot(DHARMa::simulateResiduals(mod1))
diagnose(mod1)
DHARMa::testDispersion(mod1)

lmtest::lrtest(mod1, mod1alt)

# conclusion: use conway maxwell poisson! 

# mod1: offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6
# mod2: offspring ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6
# mod3: offspring ~ sex + site + PC2 + PC3 + PC4 + PC5 + PC6

# sex as fixed effect, Fgrm * sex --> what sex is it? what exactly does interaction term mean

mod1 <- glmmTMB(offspring ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                ziformula = ~ Fgrm + sex + Fgrm*sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                data = data_old, 
                family = compois)

summary(mod1)
plot(DHARMa::simulateResiduals(mod1))
diagnose(mod1)

mod2 <- glmmTMB(offspring ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                ziformula = ~ Fgrm + sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                data = data_old, 
                family = compois)
summary(mod2)
plot(DHARMa::simulateResiduals(mod2))
diagnose(mod2)

mod3 <- glmmTMB(offspring ~ sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                ziformula = ~ sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                data = data_old, 
                family = compois)
summary(mod3)
plot(DHARMa::simulateResiduals(mod3))

mod_list <- list(mod1, mod2, mod3)
AICcmodavg::aictab(mod_list)

#bbmle::AICctab(mod1,mod2, mod3)

# permutations
fixef_mod2 <- fixef(mod2)
real_coeff <- fixef_mod2$zi["Fgrm"]
real_coeff_count <- fixef_mod2$cond["Fgrm"]

# permute Fgrm values and re run
permuts <- 100
permut_slopes_zi <- vector(length = permuts)
permut_slopes_count <- vector(length = permuts)

data_old_permut <- data_old
for(i in 1:permuts){
  data_old_permut$F_shuffle <- sample(data_old_permut$Fgrm)
  shuffle_mod <- glmmTMB(offspring ~ F_shuffle + sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                         ziformula = ~ F_shuffle + sex + site + PC2 + PC3 + PC4 + PC5 + PC6, 
                         data = data_old_permut, 
                         family = compois)
  
  permut_slopes_zi[i] <- fixef(shuffle_mod)$zi["F_shuffle"]
  permut_slopes_count[i] <- fixef(shuffle_mod)$cond["F_shuffle"]
}

quants_zi <- quantile(permut_slopes_zi, probs = c(0.05, 0.95))
quants_count <- quantile(permut_slopes_count, probs = c(0.05, 0.95))

# predict POIS regression 
ggpredict(mod2, terms = "Fgrm", type = "zero_inflated")
plot(ggpredict(mod2, terms = c("Fgrm", "site"), type = "fixed"))

plot(ggpredict(mod2, terms = "Fgrm", type = 'zero_inflated'))
# newdat <- data.frame(Fgrm = seq(-0.2, 0.4, 0.1))
# newdat$predPOIS <- predict(mod1, newdata = newdat, type = "response")

par(mfrow=c(1,2))
plot(data_old$offspring~ data_old$Fgrm, pch = 19, col = alpha("gray19", alpha = 0.5), 
     xlab = "Fgrm", 
     ylab = "number of offspring")
lines(predPOIS ~ Fgrm, newdat, col = site_palette[4], lwd = 2)

hist(permut_slopes_count, main = NULL, xlab = "effect of Fgrm on inbreeding in count model", 
     xlim = c(-10, 15), col = "gray80", border = "gray80")
mtext("100 permutations")
abline(v = quants_count, col = "black", lty = 2)
abline(v = real_coeff_count, col = colors[4], lwd =1.5)

hist(permut_slopes_zi, main = NULL, xlab = "effect of Fgrm on inbreeding in zi model", 
     xlim = c(-10, 15), col = "gray80", border = "gray80")
mtext("100 permutations")
abline(v = quants_zi, col = "black", lty = 2)
abline(v = real_coeff, col = colors[4], lwd =1.5)




mod1_old <- glmmTMB(offspring ~ Fgrm * sex + (1|site) + sex, 
                ziformula = ~ Fgrm * sex + (1|site) + sex,
                data = data_old, 
                family = compois)
summary(mod1_old)
AIC(mod1_old) # 1312.724
diagnose(mod1_old)

# mod2_old <- glmmTMB(offspring ~ Fgrm + (1|sex) + (1|site) + (1|sex)*Fgrm + PC1,
#                 ziformula = ~ Fgrm + (1|sex) + (1|site) + (1|sex)*Fgrm + PC1,
#                 data = data_old, 
#                 family = compois)
# summary(mod2_old)
# AIC(mod2_old) # 1313.428
# 
# mod3_old <- glmmTMB(offspring ~ Fgrm + (1|sex) + (1|site)  + PC1,
#                 ziformula = ~ Fgrm + (1|sex) + (1|site) + PC1,
#                 data = data_old, 
#                 family = compois)
# summary(mod3_old)
# AIC(mod3_old) # 1314.443

mod4_old <- glmmTMB(offspring ~ Fgrm + (1|sex) + (1|site),
                    ziformula = ~ Fgrm + (1|sex) + (1|site),
                    data = data_old, 
                    family = compois)
summary(mod4_old)
AIC(mod4_old) # 1313.67

mod5_old <- glmmTMB(offspring ~ (1|sex) + (1|site),
                    ziformula = ~ (1|sex) + (1|site),
                    data = data_old, 
                    family = compois)
summary(mod5_old)
AIC(mod5_old) # 1316.129

mod_old.AIC = AIC(mod1_old,mod4_old,mod5_old)

mod_old.AIC[order(mod_old.AIC$AIC),]

mod_old.AIC[order(mod_old.AIC$AIC),"AIC"][1] - mod_old.AIC[order(mod_old.AIC$AIC),"AIC"]

AICcmodavg::aictab(mod1_old,mod4_old,mod5_old)

subset()

plot(data$SVL_by~data$BirthYear, ylim = c(2005, 2020))
abline(a = 0, b = 1)
# ------------------------------------------------------------------------------------------------------------

### Supplemental: correlation between pedigree and genomic inbreeding 
# ------------------------------------------------------------------------------------------------------------

# Correlation between birth year and inferred birth year
# ------------------------------------------------------------------------------------------------------------
png(file = "../figures/supplemental_by_vs_svlby.png", width = 8, height = 8, units = "in", res = 500)
plot(x = data$BirthYear, y = data$estBirthYear, pch = 19, 
     col = alpha("black", alpha = 0.15), 
     xlab = "Known birth year", 
     ylab = "Inferred birth year", 
     ylim = c(2005, 2020))
abline(a = 0, b = 1)
dev.off()

# ------------------------------------------------------------------------------------------------------------

# Error with SVL
# ------------------------------------------------------------------------------------------------------------
PCC_by_est_max <- as.data.frame(matrix(nrow = nrow(PCC_pedigree_flt), ncol = 3))
colnames(PCC_by_est_max) <- c("ID", "SVL", "InferBirthYear")
PCC_by_est_max$ID <- PCC_pedigree_flt$id
PCC_by_est_max <- merge(x = PCC_by_est_max, y = PCC_results[[2]], all.x = TRUE, by = "ID")
# re-do inference with larger sizes
for(i in 1:length(PCC_by_est_max$ID)){
  #tryCatch(
  #{
  ind = PCC_by_est_max$ID[i]
  no_caps <- dim(PCC_exp_metaData[[ind]])[1]
  if(no_caps > 1){
    # chose latest capture year... 
    focal_year <- PCC_exp_metaData[[ind]][which(PCC_exp_metaData[[ind]]$year == max(PCC_exp_metaData[[ind]]$year)),"year"]
    SVL <- PCC_exp_metaData[[ind]][which(PCC_exp_metaData[[ind]]$year == max(PCC_exp_metaData[[ind]]$year)),"SVL"]
    
    # if the snake was captured more than once in the focal_year... 
    if(length(SVL) > 1){
      SVL <- mean(as.numeric(SVL), na.rm = T) # take the mean SVL value
    }
    if(length(focal_year) > 1){
      focal_year <- focal_year[1] # and make sure length(focal_year) == 1
    }
    # if there isn't an SVL measurement for that year... 
    if(is.na(as.numeric(subset(PCC_exp_metaData[[ind]], year == focal_year, select = SVL)$SVL))){
      # chose the largest SVL that isn't NA, and redefine SVL and focal_year
      SVL <- PCC_exp_metaData[[ind]][which(PCC_exp_metaData[[ind]]$SVL == max(as.numeric(PCC_exp_metaData[[ind]]$SVL), na.rm = T)),"SVL"]
      focal_year <-  PCC_exp_metaData[[ind]][which(PCC_exp_metaData[[ind]]$SVL == max(as.numeric(PCC_exp_metaData[[ind]]$SVL), na.rm = T)),"year"]
    }
  }else if (no_caps ==1){
    SVL <- as.numeric(PCC_exp_metaData[[ind]]$SVL)
    focal_year <- as.numeric(PCC_exp_metaData[[ind]]$year)
  }
  if(!is.na(SVL)){
    PCC_by_est_max[i, "InferBirthYear"] <- infer_birth_year(ind, "barry", as.numeric(SVL), focal_year, subset(PCC_meta, ID == ind, select = Sex)$Sex)
    PCC_by_est_max[i, "SVL"] <- SVL
  }else{
    PCC_by_est_max[i, "InferBirthYear"] <- NA
    PCC_by_est_max[i, "SVL"] <- SVL
  }
  #}, 
  #warning = function(e){
  #  cat("warning occured at interation", i, ": ", conditionMessage(e), "\n")
  #}
  #)
}  # warnings checked and can be ignored in this case

ELF_by_est_max <- as.data.frame(matrix(nrow = nrow(ELF_pedigree_flt2), ncol = 3))
colnames(ELF_by_est_max) <- c("ID", "SVL", "InferBirthYear")
ELF_by_est_max$ID <- ELF_pedigree_flt2$id
ELF_by_est_max <- merge(x = ELF_by_est_max, y = ELF_results[[2]], all.x = TRUE, by = "ID")

for(i in 1:length(ELF_by_est_max$ID)){
  #for(i in 1:141){
  
  #tryCatch(
  #{
  ind = ELF_by_est_max$ID[i]
  no_caps <- dim(ELF_exp_metaData[[ind]])[1]
  if(no_caps > 1){
    # chose first capture year... 
    focal_year <- ELF_exp_metaData[[ind]][which(ELF_exp_metaData[[ind]]$year == max(ELF_exp_metaData[[ind]]$year)),"year"]
    SVL <- ELF_exp_metaData[[ind]][which(ELF_exp_metaData[[ind]]$year == max(ELF_exp_metaData[[ind]]$year)),"SVL"]
    
    # if the snake was captured more than once in the focal_year... 
    if(length(SVL) > 1){
      SVL <- mean(as.numeric(SVL), na.rm = T) # take the mean SVL value
    }
    if(length(focal_year) > 1){
      focal_year <- focal_year[1] # and make sure length(focal_year) == 1
    }
    # if there isn't an SVL measurement for that year... 
    if(is.na(as.numeric(SVL))){
      # chose the smallest SVL that isn't NA, and redefine SVL and focal_year
      SVL <- ELF_exp_metaData[[ind]][which(as.numeric(ELF_exp_metaData[[ind]]$SVL) == max(as.numeric(ELF_exp_metaData[[ind]]$SVL), na.rm = T)),"SVL"]
      focal_year <-  ELF_exp_metaData[[ind]][which(as.numeric(ELF_exp_metaData[[ind]]$SVL) == max(as.numeric(ELF_exp_metaData[[ind]]$SVL), na.rm = T)),"year"]
      if(length(SVL) == 0){
        SVL <- NA
      }
    }
  }else if (no_caps ==1){
    SVL <- as.numeric(ELF_exp_metaData[[ind]]$SVL)
    focal_year <- as.numeric(ELF_exp_metaData[[ind]]$year)
  }
  if(!is.na(SVL) & !is.nan(SVL)){
    ELF_by_est_max[i, "InferBirthYear"] <- infer_birth_year(ind, "cass", as.numeric(SVL), focal_year, subset(ELF_meta, ID == ind, select = Sex)$Sex)
    ELF_by_est_max[i, "SVL"] <- SVL
  }else{
    ELF_by_est_max[i, "InferBirthYear"] <- NA
    ELF_by_est_max[i, "SVL"] <- SVL
  }
  #}, 
  #warning = function(e){
  #  cat("warning occured at interation", i, ": ", conditionMessage(e), "\n")
  #}
  #)
}  # warnings checked and can be ignored in this case

by_est_max <- rbind.data.frame(PCC_by_est_max, ELF_by_est_max)

#pdf(file = "../figures/supplemental_by_infer_error.pdf", width = 8, height = 8)
png(file = "../figures/supplemental_by_infer_error.png", width = 8, height = 8, units = "in", res = 500)

plot(x = by_est_max$SVL, y = (as.numeric(by_est_max$BirthYear) - by_est_max$InferBirthYear), 
     pch = 19, 
     col = alpha("black", alpha = 0.25), 
     xlab = "SVL", 
     ylab = "Error in inferred birth year")
abline(h = 0, lty = 2)
abline(h = -2, lty = 3)
abline(h = 2, lty = 3)
abline(v = 50, lty = 2, col = "red")
dev.off()



# ------------------------------------------------------------------------------------------------------------

# Inbreeding at each site
# ------------------------------------------------------------------------------------------------------------

pdf(file = "../figures/supplemental_Fgrm_dist_lines.pdf", width = 8, height = 4)
par(mfrow=c(1,2))
hist(subset(data, site == "PCC", select = Fgrm)$Fgrm, 
     xlab = "Fgrm", main = "Barry")
abline(v = 0.05, col = colors[1], lty =2, lwd = 2)
hist(subset(data, site == "ELF", select = Fgrm)$Fgrm, 
     xlab = "Fgrm", main = "Cass")
abline(v = 0.05, col = colors[1], lty =2, lwd = 2)
dev.off()

pdf(file = "../figures/supplemental_Fgrm_dist.pdf", width = 8, height = 5)
par(mfrow=c(1,2))
hist(subset(data, site == "PCC", select = Fgrm)$Fgrm, 
     xlab = "Fgrm", main = "Barry")
hist(subset(data, site == "ELF", select = Fgrm)$Fgrm, 
     xlab = "Fgrm", main = "Cass")
dev.off()

# ------------------------------------------------------------------------------------------------------------

# Joint PCA
# ------------------------------------------------------------------------------------------------------------

pca_colors <- data$site
pca_colors[pca_colors == "ELF"] <- colors[2]
pca_colors[pca_colors == "PCC"] <- colors[1]

#pdf(file = "../figures/supplemental_joint_pca.pdf", width = 8, height = 8)
png(file = "../figures/supplemental_joint_pca.png", width = 8, height = 8, units = "in", res = 500)

plot(x = data$PC1, y = data$PC2, 
     xlab = paste0("PC1, ", round(joint_PCA[[2]][1]*100, 1), "% of variation explained"), 
     ylab = paste0("PC2, ", round(joint_PCA[[2]][2]*100, 1), "% of variation explained"), 
     pch = 19, 
     col = alpha(pca_colors, alpha = 0.5))
legend("bottomright", legend = c("Barry", "Cass"), pch = 19, col = unique(pca_colors))
dev.off()

# additional version that contains outlier individual for supplement

load("../vcf_filtering/joint_pca_w_outlier.RObj", verbose = T) # list(bb_pca, var_explained)

pca_colors <- grepl("ELF", rownames(joint_PCA_outlier[[1]]))
pca_colors[pca_colors == TRUE] <- colors[2]
pca_colors[pca_colors == FALSE] <- colors[1]

png(file = "../figures/supplemental_joint_pca.png", width = 8, height = 8, units = "in", res = 500)
plot(x = joint_PCA_outlier[[1]][,1], y = joint_PCA_outlier[[1]][,2], 
     xlab = paste0("PC1, ", round(joint_PCA_outlier[[2]][1]*100, 1), "% of variation explained"), 
     ylab = paste0("PC2, ", round(joint_PCA_outlier[[2]][2]*100, 1), "% of variation explained"), 
     pch = 19, 
     col = alpha(pca_colors, alpha = 0.5))

points(x = joint_PCA_outlier[[1]][which(rownames(joint_PCA_outlier[[1]]) == "ELF_269"),1], 
       y = joint_PCA_outlier[[1]][which(rownames(joint_PCA_outlier[[1]]) == "ELF_269"),2], 
       col = "red")
legend("bottomright", legend = c("Barry", "Cass"), pch = 19, col = unique(pca_colors))
dev.off()


# ------------------------------------------------------------------------------------------------------------

# genetic IBD
# ------------------------------------------------------------------------------------------------------------
dim(PCC_centroid_dist) # dim: 256 256
dim(PCC_PWP) # dim: 261 261

# filter out individual w/o spatial information
PCC_PWP_flt <- PCC_PWP[match(colnames(PCC_centroid_dist), colnames(PCC_PWP)),match(colnames(PCC_centroid_dist), colnames(PCC_PWP))]

dim(ELF_centroid_dist) # dim: 621 621
dim(ELF_PWP) # dim: 777 777

# filter out individual w/o spatial information
ELF_PWP_flt <- ELF_PWP[match(colnames(ELF_centroid_dist), colnames(ELF_PWP)),match(colnames(ELF_centroid_dist), colnames(ELF_PWP))]
# what about with no missing data? 
ELF_PWP_nMD_flt <- ELF_PWP_nMD[match(colnames(ELF_centroid_dist), colnames(ELF_PWP)),match(colnames(ELF_centroid_dist), colnames(ELF_PWP))]

tiff(file = "../figures/supplemental_genetic_IBD.tiff", width = 8, height = 4, units = "in", res = 480)
par(mfrow=c(1,2))
scatter.smooth(PCC_PWP_flt[upper.tri(PCC_PWP_flt, diag = FALSE)] ~ PCC_centroid_dist[upper.tri(PCC_centroid_dist, diag = FALSE)], span = 2/3, degree = 2, pch = 19, col = alpha(colors[1], alpha = 0.15), 
               xlab = "Pairwise distance between capture centroids", ylab = "Pairwise pi at polymorphic sites", main = "Barry")
scatter.smooth(ELF_PWP_flt[upper.tri(ELF_PWP_flt, diag = FALSE)] ~ ELF_centroid_dist[upper.tri(ELF_centroid_dist, diag = FALSE)], span = 2/3, degree = 2, pch = 19, col = alpha(colors[2], alpha = 0.15), 
               xlab = "Pairwise distance between capture centroids", ylab = "Pairwise pi at polymorphic sites", main = "Cass")

dev.off()


plot(ELF_PWP_flt[upper.tri(ELF_PWP_flt, diag = FALSE)] ~ ELF_centroid_dist[upper.tri(ELF_centroid_dist, diag = FALSE)], pch = 19, col = alpha(colors[2], alpha = 0.15), 
               xlab = "Pairwise distance between capture centroids", ylab = "Pairwise pi at polymorphic sites", main = "Cass")

points(PWP~geo_dist, data = ELF_unique_rents_rel)

# no missing data, use ELF_PWP_nMD_flt
plot(ELF_PWP_nMD_flt[upper.tri(ELF_PWP_nMD_flt, diag = FALSE)] ~ ELF_centroid_dist[upper.tri(ELF_centroid_dist, diag = FALSE)], pch = 19, col = alpha(colors[2], alpha = 0.15), 
     xlab = "Pairwise distance between capture centroids", ylab = "Pairwise pi at polymorphic sites, no MD", main = "Cass")

points(PWP_nMD~geo_dist, data = ELF_unique_rents_rel)

# ------------------------------------------------------------------------------------------------------------

# compare to Fgrm 
# ------------------------------------------------------------------------------------------------------------
mod <- lm(Fgrm~ped_F, data = data[which(data$Fgrm < 0.8),])
summary(mod)

#pdf(file = "../figures/supplemental_pedF_vs_Fgrm.pdf", width = 8, height = 8)
png(file = "../figures/supplemental_pedF_vs_Fgrm.png", width = 7, height = 7, unit = "in", res = 500)
plot(x = data[which(data$Fgrm < 0.8),]$ped_F, y = data[which(data$Fgrm < 0.8),]$Fgrm, 
     xlab = "pedigree inbreeding", ylab = "Fgrm", 
     col = alpha("black", alpha = 0.25), pch = 19)
abline(mod, col = "red", lty = 2)
abline(h = 0, col = "black", lty = 2)
dev.off()


cor(cbind.data.frame(data$Fgrm, data$ped_F))
cor.test(data$Fgrm, data$ped_F)
# ------------------------------------------------------------------------------------------------------------


# ELF_269 offspring
# ------------------------------------------------------------------------------------------------------------
het_quant <- quantile(subset(data, site == "ELF", select = het)$het, probs = c(0.75, 0.8, 0.85, 0.9, 0.95))

pdf(file = "../figures/supplemental_cass_het.pdf", width = 8, height = 8)

hist(subset(data, site == "ELF", select = het)$het, 
     xlab = "heterozygosity", 
     main = NULL)
abline(v = subset(data, id == "ELF_473", select = het)$het, lty = 1, col = "black")
abline(v = subset(data, id == "ELF_534", select = het)$het, lty = 1, col = "black")
abline(v = subset(data, id == "ELF_666", select = het)$het, lty = 1, col = "black")
abline(v = subset(data, id == "ELF_774", select = het)$het, lty = 1, col = "black")
abline(v = subset(data, id == "ELF_828", select = het)$het, lty = 1, col = "black")

abline(v = het_quant[c(3)], col = "red", lty = 2)

legend("topleft", legend = "85th quantile", lty = 2, col = "red")
dev.off()

# ------------------------------------------------------------------------------------------------------------

# Pedigree confidence estimates
# ------------------------------------------------------------------------------------------------------------
# ELF.ped.conf
# PCC.ped.conf

pdf(file = "../figures/supplemental_ped_errors.pdf", width = 8, height = 10)
par(mfcol=c(4,3))

### mismatches
barplot(height = PCC.ped.conf[["PedErrors"]][,3,2], 
        main = "Barry sire mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes",         
        col = colors[1], 
        las = 2)

barplot(height = PCC.ped.conf[["PedErrors"]][,3,1], 
        main = "Barry dam mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes", 
        col = colors[1], 
        las = 2)

barplot(height = ELF.ped.conf[["PedErrors"]][,3,2], 
        main = "Cass sire mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes", 
        col = colors[2], 
        las = 2)

barplot(height = ELF.ped.conf[["PedErrors"]][,3,1], 
        main = "Cass dam mismatches", 
        ylab = "Proportion mismatch", 
        xlab = "individual + parent codes", 
        col = colors[2], 
        las = 2)
### false negatives
barplot(height = PCC.ped.conf[["PedErrors"]][,1,2], 
        main = "Barry sire false negatives", 
        ylab = "False negative rate", 
        xlab = "individual + parent codes",         
        col = colors[1], 
        las = 2)

barplot(height = PCC.ped.conf[["PedErrors"]][,1,1], 
        main = "Barry dam false negatives", 
        ylab = "False negative rate", 
        xlab = "individual + parent codes", 
        col = colors[1], 
        las = 2)

barplot(height = ELF.ped.conf[["PedErrors"]][,1,2], 
        main = "Cass sire false negatives", 
        ylab = "False negative rate", 
        xlab = "individual + parent codes", 
        col = colors[2], 
        las = 2)

barplot(height = ELF.ped.conf[["PedErrors"]][,1,1], 
        main = "Cass dam false negatives", 
        ylab = "False negative rate", 
        xlab = "individual + parent codes", 
        col = colors[2], 
        las = 2)
### 
barplot(height = PCC.ped.conf[["PedErrors"]][,2,2], 
        main = "Barry Sire false positives", 
        ylab = "False positive rate", 
        xlab = "individual + parent codes",         
        col = colors[1], 
        las = 2)

barplot(height = PCC.ped.conf[["PedErrors"]][,2,1], 
        main = "Barry dam false positives", 
        ylab = "False positive rate", 
        xlab = "individual + parent codes", 
        col = colors[1], 
        las = 2)
barplot(height = ELF.ped.conf[["PedErrors"]][,2,2], 
        main = "Cass sire false positives", 
        ylab = "False positive rate", 
        xlab = "individual + parent codes", 
        col = colors[2], 
        las = 2, 
        yaxt = "n")
axis(2)
barplot(height = ELF.ped.conf[["PedErrors"]][,2,1], 
        main = "Cass dam false positives", 
        ylab = "False positive rate", 
        xlab = "individual + parent codes", 
        col = colors[2], 
        las = 2, 
        yaxt = "n")
axis(2)

### 
dev.off()

write.table(PCC.ped.conf[["ConfProb"]], file = "../figures/PCC_conf_prob.txt")
write.table(ELF.ped.conf[["ConfProb"]], file = "../figures/ELF_conf_prob.txt")

# ------------------------------------------------------------------------------------------------------------

# Fgrm in count model
# ------------------------------------------------------------------------------------------------------------
# load data to plot reproductive output model 
load("../inbreeding_models/Fgrm_plot.Robj", verbose = T)
Fgrm_effects <- plot_list[[1]] # output of allEffect()
RO_mod_data <- plot_list[[2]] # "data_flt", raw unscaled data used for inbreeding models 

tiff(file="../figures/Supplemental_countFgrm.tiff", width=5, height=5, units="in", res=400)
# number of offspring
# cols <- unlist(lapply(RO_mod_data$site, function(x){if(x=="PCC"){return(colors[2])}else{return(colors[6])}}))
plot(offspring~Fgrm, data = RO_mod_data, pch = 19, col = alpha(colors[2], alpha = 0.5), 
     xlab = "Fgrm", ylab = "Number of pedigree offspring", 
     cex.axis = 1, 
     cex.lab = 1)
lines(fit~Fgrm, data = Fgrm_effects)
polygon(x=c(Fgrm_effects$Fgrm, rev(Fgrm_effects$Fgrm)), y = c(Fgrm_effects$upper, rev(Fgrm_effects$lower)), 
        col = alpha("gray", alpha = 0.5), border = NA)
dev.off()

# for talks

tiff(file="../figures/raw_reproductive_output.tiff", width=6, height=6, units="in", res=400)
par(mar = c(5, 4.2, 2.5, 4) + 0.2)             
plot(offspring~Fgrm, data = RO_mod_data, pch = 19, col = alpha(colors[2], alpha = 0.5), 
     xlab = "Fgrm", ylab = "Number of pedigree offspring", 
     cex.axis = 1.2, 
     cex.lab = 1.75)
dev.off()


tiff(file="../figures/count_reproductive_output.tiff", width=6, height=6, units="in", res=400)
par(mar = c(5, 4.2, 2.5, 4) + 0.2)             
plot(offspring~Fgrm, data = RO_mod_data, pch = 19, col = alpha(colors[2], alpha = 0.5), 
     xlab = "Fgrm", ylab = "Number of pedigree offspring", 
     cex.axis = 1.2, 
     cex.lab = 1.75)
lines(fit~Fgrm, data = Fgrm_effects)
polygon(x=c(Fgrm_effects$Fgrm, rev(Fgrm_effects$Fgrm)), y = c(Fgrm_effects$upper, rev(Fgrm_effects$lower)), 
        col = alpha("gray", alpha = 0.5), border = NA)

dev.off()

# ------------------------------------------------------------------------------------------------------------

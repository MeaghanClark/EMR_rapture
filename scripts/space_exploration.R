
## load libraries 
library(sf)
library(MetBrewer)
library(maps)
library(sequoia)
library(stringr)

## define color palette 
palette = "OKeeffe2"
met.brewer(palette, n=264, type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = 264)

# date <- "05262023"
date <- "04132024"

# load functions and data 
# ------------------------------------------------------------------------------------------------------------

source("./pedigree_funcs.R") 
  # key functions: get_coords() and plot_snakes()

find_distance <- function(pairs_mat, dist_mat){
  pairs_dist <- vector(length = nrow(pairs_mat))
  for(i in 1:nrow(pairs_mat)){
    pairs_dist[i] <- dist_mat[pairs_mat[i,][1],pairs_mat[i,][2]]
  }
  return(pairs_dist)
}

load(paste0("../pedigree_reconstruction/ELF_pedigree_results_", date, ".Robj"), verbose = T)
load(paste0("../pedigree_reconstruction/PCC_pedigree_results_", date, ".Robj"), verbose = T)
PCC_pedigree <- PCC_results[[4]]$Pedigree
PCC_pedigree$id <- str_remove(PCC_pedigree$id, "_P[0-9]")
PCC_pedigree$dam <- str_remove(PCC_pedigree$dam, "_P[0-9]")
PCC_pedigree$sire <- str_remove(PCC_pedigree$sire, "_P[0-9]")

ELF_pedigree <- ELF_results[[4]]$Pedigree
ELF_pedigree$id <- str_remove(ELF_pedigree$id, "_P[0-9]")
ELF_pedigree$dam <- str_remove(ELF_pedigree$dam, "_P[0-9]")
ELF_pedigree$sire <- str_remove(ELF_pedigree$sire, "_P[0-9]")

load(paste0("../pedigree_exploration/ELF_expanded_metaData.Robj"), verbose = T) # ELF_exp_metaData
load(paste0("../pedigree_exploration/PCC_expanded_metaData.Robj"), verbose = T) # PCC_exp_metaData

# loading coordinates from saved obj, but code to generate is below
load(paste0("../pedigree_exploration/coords",date,".Robj"), verbose = TRUE)
PCC_coords <- coords[[1]]
ELF_coords <- coords[[2]]

# PWP data 
load("../vcf_filtering/pwp_02082024.Robj", verbose = T)
PCC_PWP <- PWP[[1]]
colnames(PCC_PWP) <- str_remove(colnames(PCC_PWP), "_P[0-9]") # fix tech dupe names
rownames(PCC_PWP) <- str_remove(rownames(PCC_PWP), "_P[0-9]")

ELF_PWP <- PWP[[2]]
colnames(ELF_PWP) <- str_remove(colnames(ELF_PWP), "_P[0-9]") # fix tech dupe names
rownames(ELF_PWP) <- str_remove(rownames(ELF_PWP), "_P[0-9]")

# ------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------
### 0.5: Better site map
# ------------------------------------------------------------------------------------------------------------

pdf(file = "../pedigree_exploration/site_map.pdf", width = 12, height = 6)
maps::map(database= "state", regions = "michigan")
maps::map(database= "state", regions = c("wisconsin", "illinois", "ohio", "indiana"), add = TRUE, col = "gray10")
maps::map(database= "county", regions = "michigan,barry", add = TRUE, col = colors[1], fill = T)
maps::map(database= "county", regions = "michigan,cass", add = TRUE, col = colors[2], fill = T )
dev.off()

# ------------------------------------------------------------------------------------------------------------

### 1: Get coordinates for all individuals in pedigrees -----------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
PCC_coords <- get_coords(exp_metaData =PCC_exp_metaData, target_inds =PCC_pedigree$id, easting_col ="DD.easting", northing_col ="DD.northing", doy = FALSE)

# filter coords
PCC_coords_flt <- na.omit(PCC_coords)
# replace easting coords that aren't negative (assuming those are errors)
PCC_coords_flt$easting[PCC_coords_flt$easting > 0] <- PCC_coords_flt$easting[PCC_coords_flt$easting > 0]*-1

PCC_coords_flt <- PCC_coords_flt[-which(PCC_coords_flt$easting < -85.3583),] # PCC_268, got to be incorrect coord
PCC_coords_flt <- PCC_coords_flt[-which(PCC_coords_flt$northing < 42.518),] # PCC_303, got to be incorrect coord

ELF_coords <- get_coords(exp_metaData =ELF_exp_metaData, target_inds =ELF_pedigree$id, easting_col ="UTM.easting", northing_col ="UTM.northing", doy = TRUE)
colnames(ELF_coords) <- c("ID", "easting", "northing", "age", "year", "doy")
ELF_coords_flt <- na.omit(ELF_coords)
ELF_coords_flt <- ELF_coords_flt[-which(ELF_coords_flt$easting > 583500),] # ELF_361, coord is ecolab... 

# hist(ELF_coords_flt$easting)
# hist(ELF_coords_flt$northing)

coords <- list(PCC_coords_flt, ELF_coords_flt)
names(coords) <- c("PCC_coords_flt", "ELF_coords_flt")
save(coords, file = paste0("../pedigree_exploration/coords", date, ".Robj"))

# ------------------------------------------------------------------------------------------------------------

### 2: Generate pairwise distance matrix -----------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
#PCC_coords <- PCC_coords_flt
#ELF_coords <- ELF_coords_flt

# use PCC_coords, ELF_coords

PCC_inds <- unique(PCC_coords$ID) # 260

ELF_inds <- unique(ELF_coords$ID) # 777

# make spatial object? 

PCC_sf <- st_as_sf(PCC_coords, coords = c("easting", "northing"))
PCC_sf <- st_set_crs(PCC_sf, "+proj=longlat +datum=WGS84")
plot(PCC_sf)

ELF_sf <- st_as_sf(ELF_coords, coords = c("easting", "northing"))
ELF_sf <- st_set_crs(ELF_sf, "+proj=utm +zone=16 +datum=WGS84")
plot(ELF_sf)

# find centroid of each snake's capture locations 

# covert to multipoints when appropriate 
PCC_multi_sf <- aggregate(PCC_sf["geometry"], by = list(PCC_sf$ID), FUN = function(x) {
  st_multipoint(as.matrix(x))
})

ELF_multi_sf <- aggregate(ELF_sf["geometry"], by = list(ELF_sf$ID), FUN = function(x) {
  st_multipoint(as.matrix(x))
})

# Calculate centroids for each group
PCC_centroids <- st_centroid(PCC_multi_sf)
ELF_centroids <- st_centroid(ELF_multi_sf)

# Plot the GPS points, this isn't very informative
group_palette <- met.brewer("Archambault", n = length(unique(PCC_sf$ID)))  # Add more colors if needed
plot(PCC_sf$geometry[which(PCC_sf$ID == "PCC_257")], pch = 16, col = group_palette, main = "GPS Points")
plot(PCC_centroids$geometry[which(PCC_centroids$Group.1 == "PCC_257")], pch = 17, col = "black", add = TRUE)

group_palette <- met.brewer("Archambault", n = length(unique(ELF_sf$ID)))  # Add more colors if needed
plot(ELF_sf$geometry[which(ELF_sf$ID == "ELF_533")], pch = 16, col = group_palette, main = "GPS Points")
plot(ELF_centroids$geometry[which(ELF_centroids$Group.1 == "ELF_533")], pch = 17, col = "black", add = TRUE)

# sapply(ELF_sf$ID[duplicated(ELF_sf$ID)], FUN = function(x) {sum(ELF_sf$ID == x)}) # find individuals with multiple recaptures

# IN PROGRESS record average distance between the centroid and capture locations
# PCC_sf$geometry[which(PCC_sf$ID == "PCC_66")]

# results in new columns for coords objects: 2 columns for lat and long of centroid, 1 column for average distance from centroid

# calculate pairwise distance between snakes (using centroid locations)
PCC_centroid_dist <- st_distance(PCC_centroids) # dim: 253, 253
colnames(PCC_centroid_dist) <- PCC_centroids$Group.1
rownames(PCC_centroid_dist) <- PCC_centroids$Group.1

ELF_centroid_dist <- st_distance(ELF_centroids)
colnames(ELF_centroid_dist) <- ELF_centroids$Group.1
rownames(ELF_centroid_dist) <- ELF_centroids$Group.1
# ------------------------------------------------------------------------------------------------------------
# save centroids 
centroids <- list(PCC_centroids, ELF_centroids)
save(centroids, file = "../space_exploration/centroids.Robj")

### 2: PWP vs geo dist -----------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
site_palette <- met.brewer("Archambault", n = 6)

# PCC ----------------------------------------------------------------------------
dim(PCC_centroid_dist) # dim: 250 250
dim(PCC_PWP) # dim: 261 261

# filter out individual w/o spatial information
PCC_PWP_flt <- PCC_PWP[match(colnames(PCC_centroid_dist), colnames(PCC_PWP)),match(colnames(PCC_centroid_dist), colnames(PCC_PWP))]

scatter.smooth(PCC_PWP_flt[upper.tri(PCC_PWP_flt, diag = FALSE)] ~ PCC_centroid_dist[upper.tri(PCC_centroid_dist, diag = FALSE)], span = 2/3, degree = 2, pch = 19, col = alpha(site_palette[1], alpha = 0.15), 
               xlab = "Pairwise distance between capture centroids", ylab = "Pairwise pi at polymorphic sites", main = "Barry")
# ELF ----------------------------------------------------------------------------
dim(ELF_centroid_dist) # dim: 629 629
dim(ELF_PWP) # dim: 786 786

# filter out individual w/o spatial information
ELF_PWP_flt <- ELF_PWP[match(colnames(ELF_centroid_dist), colnames(ELF_PWP)),match(colnames(ELF_centroid_dist), colnames(ELF_PWP))]


scatter.smooth(ELF_PWP_flt[upper.tri(ELF_PWP_flt, diag = FALSE)] ~ ELF_centroid_dist[upper.tri(ELF_centroid_dist, diag = FALSE)], span = 2/3, degree = 2, pch = 19, col = alpha(colors[2], alpha = 0.15), 
               xlab = "Pairwise distance between capture centroids", ylab = "Pairwise pi at polymorphic sites", main = "Cass")

pdf(paste0("../space_exploration/pi_vs_dist_", date, ".pdf"), width = 12, height = 6)
par(mfrow=c(1,2))
scatter.smooth(PCC_PWP_flt[upper.tri(PCC_PWP_flt, diag = FALSE)] ~ PCC_centroid_dist[upper.tri(PCC_centroid_dist, diag = FALSE)], span = 2/3, degree = 2, pch = 19, col = alpha(site_palette[1], alpha = 0.15), 
               xlab = "Pairwise distance between capture centroids", ylab = "Pairwise pi at polymorphic sites", main = "Barry")
scatter.smooth(ELF_PWP_flt[upper.tri(ELF_PWP_flt, diag = FALSE)] ~ ELF_centroid_dist[upper.tri(ELF_centroid_dist, diag = FALSE)], span = 2/3, degree = 2, pch = 19, col = alpha(site_palette[2], alpha = 0.15), 
               xlab = "Pairwise distance between capture centroids", ylab = "Pairwise pi at polymorphic sites", main = "Cass")
dev.off()

pdf("../space_exploration/pi_vs_dist_poster.pdf", width = 13, height = 6.5)
par(mfrow=c(1,2))
cex_val = 2
scatter.smooth(PCC_PWP_flt[upper.tri(PCC_PWP_flt, diag = FALSE)] ~ PCC_centroid_dist[upper.tri(PCC_centroid_dist, diag = FALSE)], span = 2/3, degree = 2, pch = 16, col = alpha(site_palette[1], alpha = 0.15), 
               cex = 1, cex.lab = cex_val, cex.axis = cex_val, cex.main =cex_val, 
               xlab = "Pairwise distance between capture centroids", ylab = "Pairwise pi at polymorphic sites", main = "Barry")
scatter.smooth(ELF_PWP_flt[upper.tri(ELF_PWP_flt, diag = FALSE)] ~ ELF_centroid_dist[upper.tri(ELF_centroid_dist, diag = FALSE)], span = 2/3, degree = 2, pch = 16, col = alpha(site_palette[2], alpha = 0.15), 
               cex = 1, cex.lab = cex_val, cex.axis = cex_val, cex.main =cex_val, 
               xlab = "Pairwise distance between capture centroids", ylab = "Pairwise pi at polymorphic sites", main = "Cass")
dev.off()

# ------------------------------------------------------------------------------------------------------------

### 1: Are offspring found near their parents? -----------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# PCC ------------------------------------------------------------------------------------------------------------

# use PCC_centroid_dist
PCC_inds <- rownames(PCC_centroid_dist) # use centroid dist rowname order for simplicity

# distance between  individuals 
PCC_relM <- GetRelM(Pedigree = PCC_pedigree, GenBack = 2, Return = "Matrix", patmat = TRUE)
PCC_relM <- PCC_relM[grepl("PCC", colnames(PCC_relM)),grepl("PCC", colnames(PCC_relM))] # remove dummy individuals
PCC_relM_flt <- PCC_relM[match(dimnames(PCC_centroid_dist)[[1]], dimnames(PCC_relM)[[1]]),match(dimnames(PCC_centroid_dist)[[1]], dimnames(PCC_relM)[[1]])]

# calculate dam-offspring distances
dam_off_pairs <- which(PCC_relM_flt == "M", arr.ind = TRUE)
dam_off_dist <- vector(length = nrow(dam_off_pairs))
for(i in 1:nrow(dam_off_pairs)){
  dam_off_dist[i] <- PCC_centroid_dist[dam_off_pairs[i,][1],dam_off_pairs[i,][2]]
}

# calculate sire-offspring distances
sire_off_pairs <- which(PCC_relM_flt == "P", arr.ind = TRUE)
sire_off_dist <- vector(length = nrow(sire_off_pairs))
for(i in 1:nrow(sire_off_pairs)){
  sire_off_dist[i] <-  PCC_centroid_dist[sire_off_pairs[i,][1],sire_off_pairs[i,][2]]
}

# calculate parental distances
PCC_parents <- na.omit(PCC_pedigree[,c("dam", "sire")])
PCC_unique_rents <- unique(PCC_parents)
PCC_unique_rents <- PCC_unique_rents[grepl("PCC_", PCC_unique_rents[,1]) & grepl("PCC_", PCC_unique_rents[,2]),] # no dummy individuals

parent_dist <- vector(length = nrow(PCC_unique_rents))
for(i in 1:nrow(PCC_unique_rents)){
  tryCatch(
    {
      parent_dist[i] <-  PCC_centroid_dist[which(rownames(PCC_centroid_dist) == PCC_unique_rents[i,1]),
                                           which(colnames(PCC_centroid_dist) == PCC_unique_rents[i,2])]
    }, 
    error = function(e){
      cat("Error occured at interation", i, ": ", conditionMessage(e), "\n")
    }
  )
}

parent_dist[c(3, 11, 24)] <- NA

# calculate unrelated distances 

# make unsymmetrical 
PCC_relM_flt_asym <- PCC_relM_flt
PCC_relM_flt_asym[lower.tri(PCC_relM_flt_asym)] <- NA
unrel_pairs <- which(PCC_relM_flt_asym == "U", arr.ind = TRUE)

unrelated_dist <- vector(length = nrow(unrel_pairs))

for(i in 1:nrow(unrel_pairs)){
  unrelated_dist[i] <-  PCC_centroid_dist[unrel_pairs[i,][1],unrel_pairs[i,][2]]
}

mean(unrelated_dist)
mean(dam_off_dist)
mean(sire_off_dist)
mean(parent_dist, na.rm=T)

t.test(x = dam_off_dist, y = unrelated_dist)
t.test(x = sire_off_dist, y = unrelated_dist)

pdf(file = "../space_exploration/PCC_dist_comp.pdf", width = 6, height = 6)
hist(sample(unrelated_dist, size = 500, replace = TRUE), xlim = c(0,1200), xlab = "pairwise distance", main = "Barry", col = "white")
hist(dam_off_dist, add = T, col = alpha(site_palette[4], alpha = 0.4))
hist(sire_off_dist, add = T, col = alpha(site_palette[5], alpha = 0.4))
hist(parent_dist, add = T, col = alpha(site_palette[1], alpha = 0.4))
legend("topright", legend = c("unrelated, subsampled 500", "dam-offspring", "sire-offspring"), col = c("gray", site_palette[4], site_palette[5]), pch = 15)
dev.off()

pdf(file = "../space_exploration/PCC_dist_comp_alt.pdf", width = 10, height = 4)
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
hist(unrelated_dist, breaks = 50, xlim = c(0,1200), xlab = "pairwise distance [m]", main = "Barry", col = "gray90", 
     ylab = "Unrelated Frequency")
par(new = TRUE)                             # Add new plot
hist(dam_off_dist, breaks = 10, xlim = c(0,1200), ylim = c(0,40), col = alpha(site_palette[4], alpha = 0.5), axes = FALSE, xlab = "", ylab = "", main = "")
hist(sire_off_dist, add = T, col = alpha(site_palette[5], alpha = 0.5), axes = FALSE, xlab = "", ylab = "", main = "")
hist(parent_dist, breaks = 10, add = T, col = alpha(site_palette[1], alpha = 0.5), axes = FALSE, xlab = "", ylab = "", main = "")
axis(side = 4, at = pretty(c(0,40)))      # Add second axis
mtext("Related Frequency", side = 4, line = 3)             # Add second axis label
legend("topright", legend = c("unrelated", "dam-offspring", "sire-offspring", "parents"), col = c("gray90", site_palette[4], site_palette[5], site_palette[1]), pch = 15)
dev.off()


PCC_rel_dists <- list(unrelated_dist, dam_off_dist, sire_off_dist, parent_dist)
# ELF ------------------------------------------------------------------------------------------------------------

# use ELF_centroid_dist
ELF_inds <- rownames(ELF_centroid_dist) # use centroid dist rowname order for simplicity

# distance between  individuals 
ELF_relM <- GetRelM(Pedigree = ELF_pedigree, GenBack = 2, Return = "Matrix", patmat = TRUE)
ELF_relM <- ELF_relM[grepl("ELF", colnames(ELF_relM)),grepl("ELF", colnames(ELF_relM))] # remove dummy individuals
ELF_relM_flt <- ELF_relM[match(dimnames(ELF_centroid_dist)[[1]], dimnames(ELF_relM)[[1]]),match(dimnames(ELF_centroid_dist)[[1]], dimnames(ELF_relM)[[1]])]

# calculate dam-offspring distances
dam_off_pairs <- which(ELF_relM_flt == "M", arr.ind = TRUE)
dam_off_dist <- vector(length = nrow(dam_off_pairs))
for(i in 1:nrow(dam_off_pairs)){
  dam_off_dist[i] <- ELF_centroid_dist[dam_off_pairs[i,][1],dam_off_pairs[i,][2]]
}

# calculate sire-offspring distances
sire_off_pairs <- which(ELF_relM_flt == "P", arr.ind = TRUE)
sire_off_dist <- vector(length = nrow(sire_off_pairs))
for(i in 1:nrow(sire_off_pairs)){
  sire_off_dist[i] <-  ELF_centroid_dist[sire_off_pairs[i,][1],sire_off_pairs[i,][2]]
}

# calculate parental distances
ELF_parents <- na.omit(ELF_pedigree[,c("dam", "sire")])
ELF_unique_rents <- unique(ELF_parents)
ELF_unique_rents <- ELF_unique_rents[grepl("ELF_", ELF_unique_rents[,1]) & grepl("ELF_", ELF_unique_rents[,2]),] # no dummy individuals

parent_dist <- vector(length = nrow(ELF_unique_rents))
for(i in 1:nrow(ELF_unique_rents)){
  if(length(which(rownames(ELF_centroid_dist) == ELF_unique_rents[i,1])) != 0 & length(which(rownames(ELF_centroid_dist) == ELF_unique_rents[i,2])) != 0){
    parent_dist[i] <-  ELF_centroid_dist[which(rownames(ELF_centroid_dist) == ELF_unique_rents[i,1]),
                                         which(colnames(ELF_centroid_dist) == ELF_unique_rents[i,2])]
  }else{parent_dist[i] <- NA}
}

parent_dist <- na.omit(parent_dist)
# no spatial locations for 584

# calculate unrelated distances 

# make unsymmetrical 
ELF_relM_flt_asym <- ELF_relM_flt
ELF_relM_flt_asym[lower.tri(ELF_relM_flt_asym)] <- NA
unrel_pairs <- which(ELF_relM_flt_asym == "U", arr.ind = TRUE)

unrelated_dist <- vector(length = nrow(unrel_pairs))

for(i in 1:nrow(unrel_pairs)){
  unrelated_dist[i] <-  ELF_centroid_dist[unrel_pairs[i,][1],unrel_pairs[i,][2]]
}

mean(unrelated_dist)
mean(dam_off_dist)
mean(sire_off_dist)
mean(parent_dist)

t.test(x = dam_off_dist, y = unrelated_dist)
t.test(x = sire_off_dist, y = unrelated_dist)

pdf(file = "../space_exploration/ELF_dist_comp.pdf", width = 6, height = 6)
hist(sample(unrelated_dist, size = 3000, replace = TRUE), xlab = "pairwise distance", main = "Cass")
hist(dam_off_dist, add = T, col = alpha(site_palette[4], alpha = 0.5))
hist(sire_off_dist, add = T, col = alpha(site_palette[5], alpha = 0.5))
hist(parent_dist, add = T, col = alpha(site_palette[1], alpha = 0.5))
legend("topright", legend = c("unrelated, subsampled 3000", "dam-offspring", "sire-offspring"), col = c("gray", site_palette[4], site_palette[5]), pch = 15)
dev.off()

pdf(file = "../space_exploration/ELF_dist_comp_alt.pdf", width = 10, height = 4)
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
hist(unrelated_dist, breaks = 50, xlim = c(0,1200), xlab = "pairwise distance [m]", main = "Cass", col = "gray90", 
     ylab = "Unrelated Frequency")
par(new = TRUE)                             # Add new plot
hist(dam_off_dist, breaks = 25, xlim = c(0,1200), ylim = c(0,175), col = alpha(site_palette[4], alpha = 0.5), axes = FALSE, xlab = "", ylab = "", main = "")
hist(sire_off_dist, breaks = 25, add = T, col = alpha(site_palette[5], alpha = 0.5), axes = FALSE, xlab = "", ylab = "", main = "")
hist(parent_dist, breaks = 10, add = T, col = alpha(site_palette[1], alpha = 0.5), axes = FALSE, xlab = "", ylab = "", main = "")
axis(side = 4, at = pretty(c(0,175)))      # Add second axis
mtext("Related Frequency", side = 4, line = 3)             # Add second axis label
legend("topright", legend = c("unrelated", "dam-offspring", "sire-offspring", "parents"), col = c("gray90", site_palette[4], site_palette[5], site_palette[1]), pch = 15)
dev.off()

ELF_rel_dists <- list(unrelated_dist, dam_off_dist, sire_off_dist, parent_dist)

rel_dists <- list(PCC_rel_dists, ELF_rel_dists)
save(rel_dists, file = "../space_exploration/rel_dists_04132024.Robj")
# ------------------------------------------------------------------------------------------------------------

### Meiosis distance -----------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
## pairwise distance matrices 
# PCC_centroid_dist
# ELF_centroid_dist

## pairwise rel matrices 
# PCC_relM_flt
# ELF_relM_flt

# PCC -----------------------------------------------------------------------------
# calculate one meiosis pairs -- M, P
PCC_one_meiosis_pairs <- which(PCC_relM_flt == "M" | PCC_relM_flt == "P", arr.ind = TRUE)
PCC_one_meiosis <- find_distance(PCC_one_meiosis_pairs, PCC_centroid_dist)

# calculate two meiosis pairs -- GP, half-siblings, siblings (could also be considered 4... but in this case we're tracking distance from one parent, not the set)
# FS, MHS, PHS -- listed twice in symmetrical matrix 
# GP/GO -- listed once as "GP" once as "GO" in symmetrical matrix
# use PCC_relM_flt_asym
PCC_two_meiosis_pairs <- which(PCC_relM_flt_asym == "FS" | PCC_relM_flt_asym == "MHS" | PCC_relM_flt_asym == "PHS" | PCC_relM_flt_asym == "MGM" | PCC_relM_flt_asym == "MGF" | PCC_relM_flt_asym == "PGM" | PCC_relM_flt_asym == "PGF", arr.ind = TRUE)
PCC_two_meiosis <- find_distance(PCC_two_meiosis_pairs, PCC_centroid_dist)

# three meiosis pairs -- FA, FN, HA, HN
PCC_three_meiosis_pairs <- which(PCC_relM_flt == "FA" | PCC_relM_flt == "HA", arr.ind = TRUE)
PCC_three_meiosis <- find_distance(PCC_three_meiosis_pairs, PCC_centroid_dist)

# four meiosis pairs -- FC1, DFC1
PCC_four_meiosis_pairs <- which(PCC_relM_flt_asym == "FC1" | PCC_relM_flt_asym == "DFC1", arr.ind = TRUE)
PCC_four_meiosis <- find_distance(PCC_four_meiosis_pairs, PCC_centroid_dist)

# unrealted pairs 
PCC_unrel_pairs <- which(PCC_relM_flt_asym == "U", arr.ind = TRUE)
PCC_unrelated <- find_distance(PCC_unrel_pairs, PCC_centroid_dist)

PCC_meiosis_data <- cbind.data.frame(c(rep(1, length(PCC_one_meiosis)), rep(2, length(PCC_two_meiosis)), rep(3, length(PCC_three_meiosis)), rep(4, length(PCC_four_meiosis)), rep(5, length(PCC_unrelated))), c(PCC_one_meiosis, PCC_two_meiosis, PCC_three_meiosis, PCC_four_meiosis, PCC_unrelated))
colnames(PCC_meiosis_data) <- c("meioses", "dist")

boxplot(dist~meioses, data = PCC_meiosis_data, names = c("one", "two", "three", "four", "unrelated"), 
        xlab = "meiosis events", ylab = "distance [m]", 
        col = met.brewer("Archambault", n = 9)[c(2,4,5,6,7,8)])
abline(h = median(PCC_unrelated), col = met.brewer("Archambault", n = 9)[9], lty = 2, lwd = 3)

# ELF ------------------------------------------------------------------------------------------------------------
# calculate one meiosis pairs -- M, P
ELF_one_meiosis_pairs <- which(ELF_relM_flt == "M" | ELF_relM_flt == "P", arr.ind = TRUE)
ELF_one_meiosis <- find_distance(ELF_one_meiosis_pairs, ELF_centroid_dist)

# calculate two meiosis pairs -- GP, half-siblings, siblings (could also be considered 4... but in this case we're tracking distance from one parent, not the set)
# FS, MHS, PHS -- listed twice in symmetrical matrix 
# GP/GO -- listed once as "GP" once as "GO" in symmetrical matrix
# use ELF_relM_flt_asym
ELF_two_meiosis_pairs <- which(ELF_relM_flt_asym == "FS" | ELF_relM_flt_asym == "MHS" | ELF_relM_flt_asym == "PHS" | ELF_relM_flt_asym == "MGM" | ELF_relM_flt_asym == "MGF" | ELF_relM_flt_asym == "PGM" | ELF_relM_flt_asym == "PGF", arr.ind = TRUE)
ELF_two_meiosis <- find_distance(ELF_two_meiosis_pairs, ELF_centroid_dist)

# three meiosis pairs -- FA, FN, HA, HN
ELF_three_meiosis_pairs <- which(ELF_relM_flt == "FA" | ELF_relM_flt == "HA", arr.ind = TRUE)
ELF_three_meiosis <- find_distance(ELF_three_meiosis_pairs, ELF_centroid_dist)

# four meiosis pairs -- FC1, DFC1
ELF_four_meiosis_pairs <- which(ELF_relM_flt_asym == "FC1" | ELF_relM_flt_asym == "DFC1", arr.ind = TRUE)
ELF_four_meiosis <- find_distance(ELF_four_meiosis_pairs, ELF_centroid_dist)

# unrealted pairs 
ELF_unrel_pairs <- which(ELF_relM_flt_asym == "U", arr.ind = TRUE)
ELF_unrelated <- find_distance(ELF_unrel_pairs, ELF_centroid_dist)

ELF_meiosis_data <- cbind.data.frame(c(rep(1, length(ELF_one_meiosis)), rep(2, length(ELF_two_meiosis)), rep(3, length(ELF_three_meiosis)), rep(4, length(ELF_four_meiosis)), rep(5, length(ELF_unrelated))), c(ELF_one_meiosis, ELF_two_meiosis, ELF_three_meiosis, ELF_four_meiosis, ELF_unrelated))
colnames(ELF_meiosis_data) <- c("meioses", "dist")

boxplot(dist~meioses, data = ELF_meiosis_data, names = c("one", "two", "three", "four", "unrelated"), 
        xlab = "meiosis events", ylab = "distance [m]", 
        col = met.brewer("Archambault", n = 9)[c(2,4,5,6,7,8)])
abline(h = median(ELF_unrelated), col = met.brewer("Archambault", n = 9)[9], lty = 2, lwd = 3)


## figure 
pdf(file = "../space_exploration/meioses.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
boxplot(dist~meioses, data = PCC_meiosis_data, names = c("one", "two", "three", "four", "unrelated"), 
        xlab = "meiosis events", ylab = "distance [m]", 
        col = met.brewer("Archambault", n = 9)[c(2,4,5,6,7,8)], 
        main = "Barry County")
abline(h = median(PCC_unrelated), col = met.brewer("Archambault", n = 9)[9], lty = 2, lwd = 3)

boxplot(dist~meioses, data = ELF_meiosis_data, names = c("one", "two", "three", "four", "unrelated"), 
        xlab = "meiosis events", ylab = "distance [m]", 
        col = met.brewer("Archambault", n = 9)[c(2,4,5,6,7,8)], 
        main = "Cass County")
abline(h = median(ELF_unrelated), col = met.brewer("Archambault", n = 9)[9], lty = 2, lwd = 3)

dev.off()

# ------------------------------------------------------------------------------------------------------------

### Distribution of distances ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
### PCC----------------------------------------------------------------------------
# find distances between capture events 
PCC_recap_multi <- PCC_multi_sf$geometry[st_geometry_type(PCC_multi_sf) == "MULTIPOINT"] # 31 features
PCC_recal_ids <- PCC_multi_sf$Group.1[st_geometry_type(PCC_multi_sf) == "MULTIPOINT"] 

PCC_self_dist <- vector(mode = "list", length = length(PCC_recap_multi))
for(i in 1:length(PCC_recap_multi)){
  dist_mat <- st_distance(st_cast(PCC_recap_multi[i], to = "POINT"))
  PCC_self_dist[[i]] <- dist_mat[upper.tri(dist_mat, diag = FALSE)]
}

hist(unlist(PCC_self_dist))

### ELF ----------------------------------------------------------------------------
# find distances between capture events 
ELF_recap_multi <- ELF_multi_sf$geometry[st_geometry_type(ELF_multi_sf) == "MULTIPOINT"] # 31 features
ELF_recal_ids <- ELF_multi_sf$Group.1[st_geometry_type(ELF_multi_sf) == "MULTIPOINT"] 

ELF_self_dist <- vector(mode = "list", length = length(ELF_recap_multi))
for(i in 1:length(ELF_recap_multi)){
  dist_mat <- st_distance(st_cast(ELF_recap_multi[i], to = "POINT"))
  ELF_self_dist[[i]] <- dist_mat[upper.tri(dist_mat, diag = FALSE)]
}

hist(unlist(ELF_self_dist))

### add to meiosis figure? figure
PCC_meiosis_data_wS <- cbind.data.frame(c(rep(0,length(unlist(PCC_self_dist))), PCC_meiosis_data$meioses),
                 c(unlist(PCC_self_dist), PCC_meiosis_data$dist))
colnames(PCC_meiosis_data_wS) <- c("meioses", "dist")
ELF_meiosis_data_wS <- cbind.data.frame(c(rep(0,length(unlist(ELF_self_dist))), ELF_meiosis_data$meioses),
                                        c(unlist(ELF_self_dist), ELF_meiosis_data$dist))
colnames(ELF_meiosis_data_wS) <- c("meioses", "dist")

pdf(file = "../space_exploration/meioses_withSelf.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
boxplot(dist~meioses, data = PCC_meiosis_data_wS, names = c("recaptures", "one", "two", "three", "four", "unrelated"), 
        xlab = "meiosis events", ylab = "distance [m]", 
        col = met.brewer("Archambault", n = 9)[c(2,4,5,6,7,8)])
abline(h = median(PCC_unrelated), col = met.brewer("Archambault", n = 9)[9], lty = 2, lwd = 3)

boxplot(dist~meioses, data = ELF_meiosis_data_wS, names = c("recaptures", "one", "two", "three", "four", "unrelated"), 
        xlab = "meiosis events", ylab = "distance [m]", 
        col = met.brewer("Archambault", n = 9)[c(2,4,5,6,7,8)])
abline(h = median(ELF_unrelated), col = met.brewer("Archambault", n = 9)[9], lty = 2, lwd = 3)

dev.off()

# for poster
pdf(file = "../space_exploration/meioses_withSelf_poster.pdf", width = 13, height = 6.5)
par(mfrow=c(1,2))
cex_val = 2
boxplot(dist~meioses, data = PCC_meiosis_data_wS, names = c("recaptures", "one", "two", "three", "four", "unrelated"), 
        xlab = "meiosis events", ylab = "distance [m]", 
        col = met.brewer("Archambault", n = 9)[c(2,4,5,6,7,8)], 
        main = "Barry County", 
        cex.lab = cex_val, cex.main = cex_val, cex.axis = cex_val)
abline(h = median(PCC_unrelated), col = met.brewer("Archambault", n = 9)[9], lty = 2, lwd = 3)

boxplot(dist~meioses, data = ELF_meiosis_data_wS, names = c("recaptures", "one", "two", "three", "four", "unrelated"), 
        xlab = "meiosis events", ylab = "distance [m]", 
        col = met.brewer("Archambault", n = 9)[c(2,4,5,6,7,8)],
        main = "Cass County", 
        cex.lab = cex_val, cex.main = cex_val, cex.axis = cex_val)
abline(h = median(ELF_unrelated), col = met.brewer("Archambault", n = 9)[9], lty = 2, lwd = 3)

dev.off()


# self distances
pdf(file = "../space_exploration/recap_dist.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
vioplot::vioplot(unlist(PCC_self_dist), col = met.brewer("Archambault", n = 6)[1], names = "Barry County", ylab = "distance [m]")
vioplot::vioplot(unlist(ELF_self_dist), col = met.brewer("Archambault", n = 6)[2], names = "Cass County", ylab = "distance[m]")
dev.off()

# ------------------------------------------------------------------------------------------------------------

## PCA in space------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------
# load PCA scores
# load("../pedigree_exploration/PCA_scores_05262023.Robj", verbose = TRUE)
load("../vcf_filtering/PCA_loadings_04132024.Robj", verbose = TRUE)

# list, each object are the PCA scores for each pop from filter_vcf_popgen.R
pca_PCC <- pca[[2]]
pca_ELF <- pca[[1]]

# fix duplicate ids 
rownames(pca_PCC)[which(grepl("_P", rownames(pca_PCC)))][1] <- "PCC_300"
rownames(pca_PCC)[which(grepl("_P", rownames(pca_PCC)))][1] <- "PCC_302"
rownames(pca_PCC)[which(grepl("_P", rownames(pca_PCC)))][1] <- "PCC_62"

rownames(pca_ELF)[which(grepl("_P", rownames(pca_ELF)))] <- "ELF_544"

# load metadata (with coords)
# "extra" metadata for analyses
load(paste0("../pedigree_exploration/ELF_expanded_metaData.Robj"), verbose = T) # ELF_exp_metaData
load(paste0("../pedigree_exploration/PCC_expanded_metaData.Robj"), verbose = T) # PCC_exp_metaData

# Get centroid lat long
PCC_cent_df <- cbind.data.frame(PCC_centroids$Group.1, st_coordinates(PCC_centroids))

ELF_cent_df <- cbind.data.frame(ELF_centroids$Group.1, st_coordinates(ELF_centroids))

### PCC -------------------------------------------------------------------------------------------------------------
spgeo<- SpatialPoints(PCC_coords[c("easting", "northing")], proj4string=CRS("+proj=longlat +datum=WGS84")) 
box = c(left = -85.30487, bottom = 42.53295, right = -85.28, top = 42.54254)

dat <- SpatialPointsDataFrame(coords = spgeo, data = PCC_coords[,c(1,4,5)])
dat.df<-data.frame(dat)
ordered_scores_PCC_PC1 <- as.factor(pca_PCC[,1][match(dat.df$ID, rownames(pca_PCC))])
ordered_scores_PCC_PC2 <- as.factor(pca_PCC[,2][match(dat.df$ID, rownames(pca_PCC))])

dat.df <- cbind.data.frame(dat.df, ordered_scores_PCC_PC1, ordered_scores_PCC_PC2)

colnames(dat.df) <- c("ID", "age", "year", "long", "lat", "optional", "PC1", "PC2")

# ggmap::register_google("AIzaSyCJDzIzh11qvcFK3ltqGP1DgpTaf8OkGGo")
map <- get_map(location = box, maptype = "satellite", source = "google")
#map <- get_map(location = box, maptype = "terrain-background", source = "stamen")

## define color palette 
palette = "OKeeffe2"
met.brewer(palette, n=length(unique(ordered_scores_PCC_PC1)), type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = length(unique(ordered_scores_PCC_PC1)))

### Make map transparent-- run for hiding geographic locations! 
map_attributes <- attributes(map)
map_transparent <- matrix(adjustcolor(map, 
                                    alpha.f = 0.0), 
                                    nrow = nrow(map))
attributes(map_transparent) <- map_attributes
###

PCC_PC1 <-ggmap(map_transparent) + geom_point(aes(x = long, y = lat, color = PC1), 
                            data = dat.df, size = 2.5) + 
  #geom_path(aes(x = long, y = lat, color = ordered_scores_PCC_PC1), 
  #          data = dat.df) +
  theme_nothing() + #theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title="PCC PC1") + ylab("") + xlab("") + 
  scale_color_manual(values=met.brewer(palette, n = length(unique(dat.df$PC1))))

for(i in 1:length(unique(dat.df$ID))){
  if(nrow(dat.df[which(dat.df$ID == dat.df$ID[i]),]) > 1){
    PCC_PC1 <- PCC_PC1 + geom_path(aes(x = long, y = lat, color = PC1), 
                     data = dat.df[which(dat.df$ID == dat.df$ID[i]),])
  }
}

PCC_PC1

PCC_PC2 <-ggmap(map_transparent) + geom_point(aes(x = long, y = lat, color = PC2), 
                                  data = dat.df, size = 2.5) + 
  #geom_path(aes(x = long, y = lat, color = ordered_scores_PCC_PC2), 
  #          data = dat.df) +
  theme_nothing() + #theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title="PCC PC2") + ylab("") + xlab("") + 
  scale_color_manual(values=met.brewer(palette, n = length(unique(dat.df$PC2))))

for(i in 1:length(unique(dat.df$ID))){
  if(nrow(dat.df[which(dat.df$ID == dat.df$ID[i]),]) > 1){
    PCC_PC2 <- PCC_PC2 + geom_path(aes(x = long, y = lat, color = PC2), 
                                   data = dat.df[which(dat.df$ID == dat.df$ID[i]),])
  }
}

PCC_PC2


pdf(file = "../pedigree_exploration/PCC_pca_in_space_blank_updated.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
PCC_PC1
PCC_PC2
dev.off()
### 

# with centroid locations
box = c(left = -85.30487, bottom = 42.53295, right = -85.28, top = 42.54254)

dat <- SpatialPointsDataFrame(coords = PCC_cent_df[,2:3], data = PCC_cent_df)
dat.df<-data.frame(dat)
colnames(dat.df) <- c("ID", "long", "lat", "rep", "rep", "optional")

ordered_scores_PCC_PC1_cent <- as.factor(pca_PCC[,1][match(dat.df$ID, rownames(pca_PCC))])
ordered_scores_PCC_PC2_cent <- as.factor(pca_PCC[,2][match(dat.df$ID, rownames(pca_PCC))])

dat.df <- cbind.data.frame(dat.df, ordered_scores_PCC_PC1_cent)

dat.df$ordered_scores_PCC_PC1 <- as.factor(dat.df$ordered_scores_PCC_PC1_cent)

dat.df <- dat.df[,c(1,2,3,7)]
colnames(dat.df) <- c("ID", "long", "lat", "PC1")

map <- get_map(location = box, maptype = "satellite", source = "google")

### Make map transparent-- run for hiding geographic locations! 
map_attributes <- attributes(map)
map_transparent <- matrix(adjustcolor(map, 
                          alpha.f = 0.0), 
                          nrow = nrow(map))
attributes(map_transparent) <- map_attributes
###

## define color palette 
palette = "OKeeffe2"
met.brewer(palette, n=length(unique(ordered_scores_PCC_PC1_cent)), type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = length(unique(ordered_scores_PCC_PC1_cent)))

PCC_PC1_cent <-ggmap(map_transparent) + geom_point(aes(x = long, y = lat, color = PC1), 
                                            data = dat.df, size = 2.0) + 
  theme_nothing() + #theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title="") + ylab("") + xlab("") + 
  scale_color_manual(values=met.brewer(palette, n = length(unique(dat.df$PC1))))

pdf(file = "../space_exploration/PCC_PC1_centroid_blank.pdf", width = 6, height = 6)
PCC_PC1_cent
dev.off()

pdf(file = "../space_exploration/PCC_blank_map.pdf", width = 6, height = 6)
ggmap(map) +   theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(title="") + ylab("") + xlab("") 
dev.off()

pdf(file = "../space_exploration/PCC_offspring_in_space_updated_palette.pdf", width = 6, height = 6)
colors
dev.off()

### ELF -------------------------------------------------------------------------------------------------------------

sputm <- SpatialPoints(ELF_coords[c("easting", "northing")], proj4string=CRS("+proj=utm +zone=16 +datum=WGS84")) 
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
box = c(left = -86.01221, bottom = 41.94555, right = -85.98966, top = 41.95972)

dat <- SpatialPointsDataFrame(coords = spgeo, data = ELF_coords[,c(1,4,5)])
dat.df<-data.frame(dat)
ordered_scores_ELF_PC1 <- as.factor(pca_ELF[,1][match(dat.df$ID, rownames(pca_ELF))])
ordered_scores_ELF_PC2 <- as.factor(pca_ELF[,2][match(dat.df$ID, rownames(pca_ELF))])


dat.df <- cbind.data.frame(dat.df, ordered_scores_ELF_PC1, ordered_scores_ELF_PC2)

colnames(dat.df) <- c("ID", "age", "year", "long", "lat", "optional", "PC1", "PC2")

map <- get_map(location = box, maptype = "satellite", source = "google")

### Make map transparent-- run for hiding geographic locations! 
map_attributes <- attributes(map)
map_transparent <- matrix(adjustcolor(map, 
                                      alpha.f = 0.0), 
                          nrow = nrow(map))
attributes(map_transparent) <- map_attributes
###

## define color palette 
palette = "OKeeffe2"
met.brewer(palette, n=length(unique(ordered_scores_ELF_PC1)), type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = length(unique(ordered_scores_ELF_PC1)))

ELF_PC1 <-ggmap(map_transparent) + geom_point(aes(x = long, y = lat, color = PC1), 
                                  data = dat.df, size = 2.5) + 
  #geom_path(aes(x = long, y = lat, color = ordered_scores_ELF_PC1), 
  #          data = dat.df) +
  theme_nothing() + #theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title="ELF PC1") + ylab("") + xlab("") + 
  scale_color_manual(values=met.brewer(palette, n = length(unique(dat.df$PC1))))

for(i in 1:length(unique(dat.df$ID))){
  if(nrow(dat.df[which(dat.df$ID == dat.df$ID[i]),]) > 1){
    ELF_PC1 <- ELF_PC1 + geom_path(aes(x = long, y = lat, color = PC1), 
                                   data = dat.df[which(dat.df$ID == dat.df$ID[i]),])
  }
}

ELF_PC1

ELF_PC2 <-ggmap(map) + geom_point(aes(x = long, y = lat, color = PC2), 
                                  data = dat.df, size = 2.5) + 
  #geom_path(aes(x = long, y = lat, color = ordered_scores_ELF_PC2), 
  #          data = dat.df) +
  theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title="ELF PC2") + ylab("") + xlab("") + 
  scale_color_manual(values=met.brewer(palette, n = length(unique(dat.df$PC2))))

for(i in 1:length(unique(dat.df$ID))){
  if(nrow(dat.df[which(dat.df$ID == dat.df$ID[i]),]) > 1){
    ELF_PC2 <- ELF_PC2 + geom_path(aes(x = long, y = lat, color = PC2), 
                                   data = dat.df[which(dat.df$ID == dat.df$ID[i]),])
  }
}

ELF_PC2


pdf(file = "../pedigree_exploration/ELF_pca_in_space_updated.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
ELF_PC1
ELF_PC2
dev.off()

### 
# with centroid locations
sputm <- SpatialPoints(ELF_cent_df[c("X", "Y")], proj4string=CRS("+proj=utm +zone=16 +datum=WGS84")) 
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
box = c(left = -86.01221, bottom = 41.94555, right = -85.98966, top = 41.95972)

dat <- SpatialPointsDataFrame(coords = spgeo, data = ELF_cent_df)
dat.df<-data.frame(dat)
colnames(dat.df) <- c("ID", "easting", "northing", "long", "lat", "optional")

ordered_scores_ELF_PC1_cent <- as.factor(pca_ELF[,1][match(dat.df$ID, rownames(pca_ELF))])
ordered_scores_ELF_PC2_cent <- as.factor(pca_ELF[,2][match(dat.df$ID, rownames(pca_ELF))])

dat.df <- cbind.data.frame(dat.df, ordered_scores_ELF_PC1_cent)

dat.df$ordered_scores_ELF_PC1 <- as.factor(dat.df$ordered_scores_ELF_PC1)

dat.df <- dat.df[,c(1,4,5,7)]
colnames(dat.df) <- c("ID", "long", "lat", "PC1")

map <- get_map(location = box, maptype = "satellite", source = "google")

### Make map transparent-- run for hiding geographic locations! 
map_attributes <- attributes(map)
map_transparent <- matrix(adjustcolor(map, 
                          alpha.f = 0.0), 
                          nrow = nrow(map))
attributes(map_transparent) <- map_attributes
###

## define color palette 
palette = "OKeeffe2"
met.brewer(palette, n=length(unique(ordered_scores_ELF_PC1)), type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = length(unique(ordered_scores_ELF_PC1)))

ELF_PC1_cent <-ggmap(map) + geom_point(aes(x = long, y = lat, color = PC1), 
                                       data = dat.df, size = 1.5) + 
  #theme_nothing() + 
  theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title="") + ylab("") + xlab("") + 
  scale_color_manual(values=met.brewer(palette, n = length(unique(dat.df$PC1))))

pdf(file = "../space_exploration/ELF_PC1_centroid.pdf", width = 6, height = 6)
ELF_PC1_cent
dev.off()

pdf(file = "../space_exploration/ELF_blank_map.pdf", width = 6, height = 6)
ggmap(map) +   theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(title="") + ylab("") + xlab("") 
dev.off()
# --------------------------------------------------------------------------------------------------------


## Map colored by number of offspring ------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------
# get number of offspring
PCC_off <- offspring_per_ind(PCC_pedigree, plot = F)
ELF_off <- offspring_per_ind(ELF_pedigree_flt, plot = F)

# format offspring data and combine
PCC_moms <- PCC_off[[1]]
PCC_moms <- PCC_moms[which(grepl("_", names(PCC_moms)))]
PCC_dads <- PCC_off[[2]]
PCC_dads <- PCC_dads[which(grepl("_", names(PCC_dads)))]
PCC_offspring_data <- cbind.data.frame(c(names(PCC_moms), names(PCC_dads)), c(PCC_moms, PCC_dads), c(rep(1, length(PCC_moms)), rep(2, length(PCC_dads))))
colnames(PCC_offspring_data) <- c("id", "offspring", "sex")

ELF_moms <- ELF_off[[1]]
ELF_moms <- ELF_moms[which(grepl("_", names(ELF_moms)))]
ELF_dads <- ELF_off[[2]]
ELF_dads <- ELF_dads[which(grepl("_", names(ELF_dads)))]
ELF_offspring_data <- cbind.data.frame(c(names(ELF_moms), names(ELF_dads)), c(ELF_moms, ELF_dads), c(rep(1, length(ELF_moms)), rep(2, length(ELF_dads))))
colnames(ELF_offspring_data) <- c("id", "offspring", "sex")

# --------------------------------------------------------------------------------------------
# Get centroid lat long
PCC_cent_df <- cbind.data.frame(PCC_centroids$Group.1, st_coordinates(PCC_centroids))

ELF_cent_df <- cbind.data.frame(ELF_centroids$Group.1, st_coordinates(ELF_centroids))


### PCC -------------------------------------------------------------------------------------------------------------
#spgeo<- SpatialPoints(PCC_coords[c("easting", "northing")], proj4string=CRS("+proj=longlat +datum=WGS84")) 
box = c(left = -85.30487, bottom = 42.53295, right = -85.28, top = 42.54254)

dat <- SpatialPointsDataFrame(coords = PCC_cent_df[,2:3], data = PCC_cent_df)
dat.df<-data.frame(dat)
colnames(dat.df) <- c("ID", "long", "lat", "rep", "rep", "optional")

ordered_offspring <- PCC_offspring_data[match(dat.df$ID, PCC_offspring_data$id),"offspring"]
ordered_offspring[which(is.na(ordered_offspring))] <- 0

dat.df <- cbind.data.frame(dat.df, ordered_offspring)

dat.df$ordered_offspring <- as.factor(dat.df$ordered_offspring)

dat.df <- dat.df[,c(1,2,3,7)]
colnames(dat.df) <- c("ID", "long", "lat", "offspring")
dat.df <- dat.df[which(dat.df$offspring != 0),]

map <- get_map(location = box, maptype = "satellite", source = "google")

## define color palette 
palette = "OKeeffe2"
met.brewer(palette, n=length(unique(ordered_offspring)), type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = length(unique(ordered_offspring)))

PCC_offspring_map <-ggmap(map) + geom_point(aes(x = long, y = lat, color = offspring), 
                                  data = dat.df, size = 2.0) + 
  theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title="") + ylab("") + xlab("") + 
  scale_color_manual(values=met.brewer(palette, n = length(unique(dat.df$offspring))))

pdf(file = "../space_exploration/PCC_offspring_in_space_updated.pdf", width = 6, height = 6)
PCC_offspring_map
dev.off()

pdf(file = "../space_exploration/PCC_offspring_in_space_updated_palette.pdf", width = 6, height = 6)
colors
dev.off()

### ELF -------------------------------------------------------------------------------------------------------------
sputm <- SpatialPoints(ELF_cent_df[c("X", "Y")], proj4string=CRS("+proj=utm +zone=16 +datum=WGS84")) 
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
box = c(left = -86.01221, bottom = 41.94555, right = -85.98966, top = 41.95972)
#box = c(left = -86.01182, bottom = 41.94555, right = -85.9899, top = 41.95842)


dat <- SpatialPointsDataFrame(coords = spgeo, data = ELF_cent_df)
dat.df<-data.frame(dat)
colnames(dat.df) <- c("ID", "easting", "northing", "long", "lat", "optional")

ordered_offspring <- ELF_offspring_data[match(dat.df$ID, ELF_offspring_data$id),"offspring"]
ordered_offspring[which(is.na(ordered_offspring))] <- 0

dat.df <- cbind.data.frame(dat.df, ordered_offspring)

dat.df$ordered_offspring <- as.factor(dat.df$ordered_offspring)

dat.df <- dat.df[,c(1,4,5,7)]
colnames(dat.df) <- c("ID", "long", "lat", "offspring")

dat.df <- dat.df[which(dat.df$offspring != 0),]


map <- get_map(location = box, maptype = "satellite", source = "google")

## define color palette 
palette = "OKeeffe2"
met.brewer(palette, n=length(unique(ordered_offspring)), type="continuous")
MetBrewer::colorblind.friendly(palette)
colors <- met.brewer(palette, n = length(unique(ordered_offspring)))

ELF_offspring_map <-ggmap(map) + geom_point(aes(x = long, y = lat, color = offspring), 
                                            data = dat.df, size = 1.5) + 
  theme(legend.position = "none", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  labs(title="") + ylab("") + xlab("") + 
  scale_color_manual(values=met.brewer(palette, n = length(unique(dat.df$offspring))))

ELF_offspring_map

pdf(file = "../space_exploration/ELF_offspring_in_space_updated.pdf", width = 6, height = 6)
ELF_offspring_map
dev.off()

pdf(file = "../space_exploration/ELF_offspring_in_space_updated_palette.pdf", width = 6, height = 6)
colors
dev.off()

# --------------------------------------------------------------------------------------------

## Space in PCA ------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# run PCA code first! also fix this code so you don't have to do that! 
ELF_coords$easting

ELF_pca_inds <- sep_pop$ELF@ind.names
ELF_inds_w_by[which(ELF_inds_w_by$id == "ELF_544_P10"),"id"] <- "ELF_544_P5"

color_vec <- rep("gray", length(ELF_pca_inds))

ELF_coords[match(ELF_pca_inds, ELF_coords$ID),]
color_vec[match(ELF_coords$ID, ELF_pca_inds)] <- met.brewer("Greek", n = 491)

plot(x = pca_ELF$scores[, 1], y = pca_ELF$scores[, 2],
     col = alpha(color_vec, 0.7), 
     cex = 1, pch = 19, cex.lab = 1.5,
     xlab = "PC1, 3.5% variation", ylab = "PC2, 3.0% variation")




# --------------------------------------------------------------------------------------------------------
## code graveyard
# get parent offspring information

# PCC_PO <- sapply(PCC_inds, FUN = function(x){
#   tryCatch(
#     GetDescendants(x, PCC_pedigree), error=function(e){}
#   )
#   })
# 
## make vector of parent-offspring distances 
#parent_offspring_distances <- vector()
#for(i in 1:length(PCC_PO)){
#  if(!is.null(PCC_PO[i])){
#    desc_list <- PCC_PO[i]
#    par_centroid <- PCC_centroids[which(PCC_centroids$Group.1 == desc_list[[1]]$id),]
#    off <- desc_list[[1]]$offspring[grepl("PCC", desc_list[[1]]$offspring)]
#    # correct technical replicate names
#    off <- str_remove(off, "_P[0-9]")
#    
#    off_centroid <- lapply(off, FUN = function(x) {
#      PCC_centroids[which(PCC_centroids$Group.1 == x),]
#    })
#    off_dist <- unlist(lapply(off_centroid, FUN = function(x) {
#      st_distance(x = par_centroid, y = x)
#      }))
#    parent_offspring_distances <- c(parent_offspring_distances, off_dist)
#   }
# }

# retain only individuals with spatial data

distances <- matrix(nrow = length(PCC_inds), ncol = 3)
colnames(distances) <- c("dam", "sire", "avg_unrelated")
unrelated_distances <- vector()
for(i in 1:length(PCC_inds)){
  # get centroid of focal individual
  focal_centroid <- PCC_centroids[which(PCC_centroids$Group.1 == PCC_inds[i]),]
  
  # make list of relatives
  rel_vec <- na.omit(PCC_relM[which(rownames(PCC_relM) == PCC_inds[i], arr.ind = TRUE),])
  
  # parents
  dam <- names(rel_vec[rel_vec == "M"])
  sire <- names(rel_vec[rel_vec == "P"])
  
  if(length(dam) > 0){
    dam_dist <- st_distance(x = focal_centroid, y = PCC_centroids[which(PCC_centroids$Group.1 == dam),])
  }else { dam_dist <- NA }
  if(length(sire) > 0){
    sire_dist <- st_distance(x = focal_centroid, y = PCC_centroids[which(PCC_centroids$Group.1 == sire),])
  }else { sire_dist <- NA }
  
  not_rel <- names(rel_vec[rel_vec == "U"])
  # correct technical replicate names
  not_rel <- str_remove(not_rel, "_P[0-9]")
  
  # make list of dists from focal ind to unrelated individuals 
  not_rel_dist <- unlist(lapply(not_rel, FUN = function(x) {
    st_distance(x = focal_centroid, y = PCC_centroids[which(PCC_centroids$Group.1 == x),])
  }))
  distances[i,"dam"] <- dam_dist
  distances[i,"sire"] <- sire_dist
  distances[i,"avg_unrelated"] <- mean(not_rel_dist)
  unrelated_distance <- c(unrelated_distances, not_rel_dist)
}

distances <- as.data.frame(distances)




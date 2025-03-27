# sample code to plot a pedigree from Sequoia using kinship2 "plot.pedigree"
# M. Clark, 08/13/2024

library(sequoia)
library(kinship2)
library(MetBrewer)

colFunc <- function (x, cols, nCols, valRange){
  if (is.null(valRange)){
    valRange <- c(min(x), max(x))
  }
  cols <- (grDevices::colorRampPalette(cols))(nCols)[findInterval(x, seq(valRange[1], valRange[2], length.out = nCols))]
  return(cols)
}

# define color palette
ped_colors_alt <- met.brewer("Hiroshige", n = 20)

# make pedigree data frame with just desired individuals 
# ELF_pedigree_flt_nD is in the same format as a pedigree dataframe output from Sequoia
ELF_ex_ped1 <- subset(ELF_pedigree_flt_nD, id %in% ped1_inds)

# polish pedigree
ELF_ex_ped1_polish <- sequoia::PedPolish(ELF_ex_ped1, DropNonSNPd=FALSE,
                                         FillParents = TRUE)
# add sex information 
# be aware that if any individuals have unknown sex, they will be assigned "male" here
ELF_ex_ped1_polish$Sex <- ifelse(ELF_meta[match(ELF_ex_ped1_polish$id, ELF_meta$ID),"Sex"] == 1, yes = "female", no = "male")

# fill in sex info for dummy individuals
ELF_ex_ped1_polish[grepl("F0", ELF_ex_ped1_polish$id), "Sex"] <- "female" # females
ELF_ex_ped1_polish[grepl("M0", ELF_ex_ped1_polish$id), "Sex"] <- "male" # males

# Fix the sex of parents, add parents that are missing from the pedigree
ELF_ex_ped1_fix <- with(ELF_ex_ped1_polish, kinship2::fixParents(id=id, dadid=sire, momid=dam, sex=Sex))

# get kinship2 pedigree object
ELF_ex_ped1_k <- with(ELF_ex_ped1_fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))

# get Fgrm in the same order as individuals in the pedigree object
# ELF_fgrm is a dataframe with one column of Fgrm values, and rownames are the individual name 
ELF_ex_ped1_Fgrm <- ELF_fgrm[match(ELF_ex_ped1_k$id, rownames(ELF_fgrm)),1]

# generate genotype matrix to distinguish genotyped and non-genotyped individuals
genotyped <- !is.na(ELF_ex_ped1_Fgrm)
not_geno <- is.na(ELF_ex_ped1_Fgrm)
ELF_ped1_genotyped <- vector(length=length(ELF_ex_ped1_k$id))
ELF_ped1_genotyped[genotyped] <- 1
ELF_ped1_genotyped[not_geno] <- 0

# define colors for each individual based on range of color values
ELF_ped1_cols <- colFunc(ELF_ex_ped1_Fgrm, ped_colors_alt, nCols = length(all_fgrm), valRange = range(all_fgrm, na.rm =T))
ELF_ped1_cols[is.na(ELF_ex_ped1_Fgrm)] <- "gray55" # individuals who were not genotyped are gray

# make pedigree plot 
# "ELF_ped1_cols" is a vector that lists the desired color for each individual/symbol 
plot.pedigree(ELF_ex_ped1_k, id = rep("", length(ELF_ex_ped1_k$id)), symbolsize=1, 
              density = c(-1, 35, 65, 20), 
              affected = ELF_ped1_genotyped, col = ELF_ped1_cols, 
              width = 100, 
              align = c(6,2))

# create scale bar
# define evenly spaced colors for scale bar
x2 <- seq(length.out = length(all_fgrm), from = range(all_fgrm, na.rm = T)[1], to = range(all_fgrm, na.rm = T)[2])
c2 <- colFunc(x2, ped_colors_alt, nCols = length(x2), valRange = range(all_fgrm, na.rm =T))

# start plotting
colRast <- as.raster(matrix(c2,nrow=1,ncol=length(c2)))
plot(NULL, xlim = range(all_fgrm), ylim = c(0,10), axes = F, xlab = "", ylab = "")
axis(side = 1, at = c(-0.1, 0, 0.1, 0.2, 0.3, 0.315), labels = TRUE, main = "Fgrm", xlim = range(all_fgrm))
xleft = range(all_fgrm)[1]
xright = range(all_fgrm)[2]
ybottom = 0
ytop= 10
rasterImage(colRast,xleft,ybottom,xright,ytop,interpolate=TRUE)

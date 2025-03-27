## File name: Fgrm_functions.R## Purpose: Functions to calculate individual Fgrm from SNP data based on equation for FhatIII in Yang et al. 2011
## M. I. Clark, July 2023

get_pop_AF <- function(t_gt){
	# This function calculates allele frequencies for each locus in t_gt  # t_gt is a genotype matrix where each row is an individual and each column is a SNP. Genotypes are encoded as 0/1/2 or NA for missing data  p <- vector(length = ncol(t_gt))  for(i in 1:ncol(t_gt)){    p[i] <- sum(t_gt[,i], na.rm = TRUE) / (2 * sum(!is.na(t_gt[,i])))  }  return(p)}calc_Fgrm <- function(t_gt){
	# This function calculates individual Fgrm 
    # t_gt is a genotype matrix where each row is an individual and each column is a SNP. Genotypes are encoded as 0/1/2 or NA for missing data    # filter out monomorphic sites   t_gt_poly <- t_gt[,which(!colSums(t_gt, na.rm = TRUE)==0)]    # calculate population-wide allele frequencies   p <- get_pop_AF(t_gt_poly)
  
  # calculate Fgrm  Fgrm <- vector(length = nrow(t_gt_poly))  for(j in 1:nrow(t_gt_poly)){    Fgrm_i <- vector(length = ncol(t_gt_poly))    for(i in 1:ncol(t_gt_poly)){      x <- t_gt_poly[j,i]      Fgrm_i[i] <- (x^2 - (1 + 2*p[i])*x + 2*(p[i]^2)) / (2*p[i])*(1-p[i]) # sites have to be polymorphic?     }     Fgrm[j] <- sum(Fgrm_i, na.rm = TRUE) / sum(!is.na(Fgrm_i))  }  return(Fgrm)}


calcFi <- function(x, p){
	    (x^2 - (1 + 2*p)*x + 2*(p^2)) / (2*p)*(1-p) # sites have to be polymorphic? 
}


x <- vector(length = length(seq(0,1,0.05)))
j = 1
for(i in seq(0,1,0.05)){
	x[j] <- calcFi(0,i)
	j <- j + 1
}

y <- vector(length = length(seq(0,1,0.05)))
j = 1
for(i in seq(0,1,0.05)){
	y[j] <- calcFi(1,i)
	j <- j + 1
}

z <- vector(length = length(seq(0,1,0.05)))
j = 1
for(i in seq(0,1,0.05)){
	z[j] <- calcFi(2,i)
	j <- j + 1
}

plot(NULL, xlim = c(0,1), ylim = c(-1,2), xlab = "allele frequency of reference allele in population", ylab = "Fgrm for single SNP in one individual")
lines(x = seq(0,1,0.05), y = x, col = "red")
lines(x = seq(0,1,0.05), y = y, col = "blue")
lines(x = seq(0,1,0.05), y = z, col = "green")
abline(v = 0.5, col = "gray", lty = 2)
legend("topright", legend = c("homozygous for alt allele", "heterozygous", "homozygous for ref allele"), lty = 1, col = c("red", "blue", "green"))




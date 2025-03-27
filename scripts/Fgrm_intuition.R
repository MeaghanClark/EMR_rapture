calcFgrm_i <- function(x, p){
	    (x^2 - (1 + 2*p)*x + 2*(p^2)) / (2*p)*(1-p)  
}


x <- vector(length = length(seq(0,1,0.05)))
j = 1
for(i in seq(0,1,0.05)){
	x[j] <- calcFgrm_i(0,i)
	j <- j + 1
}

y <- vector(length = length(seq(0,1,0.05)))
j = 1
for(i in seq(0,1,0.05)){
	y[j] <- calcFgrm_i(1,i)
	j <- j + 1
}

z <- vector(length = length(seq(0,1,0.05)))
j = 1
for(i in seq(0,1,0.05)){
	z[j] <- calcFgrm_i(2,i)
	j <- j + 1
}

plot(NULL, xlim = c(0,1), ylim = c(-1,1), xlab = "allele frequency of reference allele in population", ylab = "Fgrm for single SNP in one individual")
lines(x = seq(0,1,0.05), y = x, col = "red")
lines(x = seq(0,1,0.05), y = y, col = "blue")
lines(x = seq(0,1,0.05), y = z, col = "green")
abline(v = 0.5, col = "gray", lty = 2)
legend("topright", legend = c("homozygous for alt allele", "heterozygous", "homozygous for ref allele"), lty = 1, col = c("red", "blue", "green"))

abline(h = 0, col = "gray", lty = 2)

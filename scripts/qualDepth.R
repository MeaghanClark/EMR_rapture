#!/usr/bin/env Rscript

flist = args[1] #"/mnt/research/Fitz_Lab/projects/massasauga/EMR_WGS/variants/bamstats/EMR_bamstats_depth.txt"
depth <- read.table(flist[1],head=TRUE)

save(depth, file = args[2]) # "/mnt/research/Fitz_Lab/projects/massasauga/EMR_WGS/variants/bamstats/EMR_bamstats_depth.Robj"

outfile_hist = arg[3]

pdf(file = outfile_hist, height = 6, width = 6)
# total depth hist: ----------------------------------------------------------
depth <- as.numeric(df[,1])
hist(depth, main = NULL,  xlab = colnames(df)[1])
abline(v = median(depth), col = "black", lwd = 2)
mtext(median(depth), at = median(depth))

abline(v = (median(depth)*1.5), col = "red", lwd = 2)
mtext((median(depth)*1.5), at = (median(depth)*1.5), col = "red")

abline(v = (median(depth)*0.5), col = "red", lwd = 2)
mtext((median(depth)*0.5), at = (median(depth)*0.5), col = "red")

abline(v = (median(depth)*0.75), col = "blue", lwd = 2)
mtext((median(depth)*0.75), at = (median(depth)*0.75), col = "blue", line = 1.5)

abline(v = (median(depth)*1.25), col = "blue", lwd = 2)
mtext((median(depth)*1.25), at = (median(depth)*1.25), col = "blue", line = 1.5)

abline(v = (median(depth)*2), col = "green", lwd = 2)
mtext((median(depth)*2), at = (median(depth)*2), col = "green")

legend("topright", c("+/- 0.5", "+/- 0.25", "x2"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

### zoom 1 ----------------------------------------------------------
hist(depth, main = NULL, xlim = c(0, 1e6), xlab = colnames(df)[1], breaks = 100)
abline(v = median(depth), col = "black", lwd = 2)
mtext(median(depth), at = median(depth))

abline(v = (median(depth)*1.5), col = "red", lwd = 2)
mtext((median(depth)*1.5), at = (median(depth)*1.5), col = "red")

abline(v = (median(depth)*0.5), col = "red", lwd = 2)
mtext((median(depth)*0.5), at = (median(depth)*0.5), col = "red")

abline(v = (median(depth)*0.75), col = "blue", lwd = 2)
mtext((median(depth)*0.75), at = (median(depth)*0.75), col = "blue")

abline(v = (median(depth)*1.25), col = "blue", lwd = 2)
mtext((median(depth)*1.25), at = (median(depth)*1.25), col = "blue")

abline(v = (median(depth)*2), col = "green", lwd = 2)
mtext((median(depth)*2), at = (median(depth)*2), col = "green")

legend("topright", c("+/- 0.5", "+/- 0.25", "x2"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

hist(depth, main = NULL, xlim = c(0, 1e6), xlab = colnames(df)[1], breaks = 1000)
abline(v = median(depth), col = "black", lwd = 2)
mtext(median(depth), at = median(depth))

abline(v = (median(depth)*1.5), col = "red", lwd = 2)
mtext((median(depth)*1.5), at = (median(depth)*1.5), col = "red")

abline(v = (median(depth)*0.5), col = "red", lwd = 2)
mtext((median(depth)*0.5), at = (median(depth)*0.5), col = "red")

abline(v = (median(depth)*0.75), col = "blue", lwd = 2)
mtext((median(depth)*0.75), at = (median(depth)*0.75), col = "blue")

abline(v = (median(depth)*1.25), col = "blue", lwd = 2)
mtext((median(depth)*1.25), at = (median(depth)*1.25), col = "blue")

abline(v = (median(depth)*2), col = "green", lwd = 2)
mtext((median(depth)*2), at = (median(depth)*2), col = "green")

legend("topright", c("+/- 0.5", "+/- 0.25", "x2"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

hist(depth, main = NULL, xlim = c(0, 1e6), xlab = colnames(df)[1], breaks = 10000)
abline(v = median(depth), col = "black", lwd = 2)
mtext(median(depth), at = median(depth))

abline(v = (median(depth)*1.5), col = "red", lwd = 2)
mtext((median(depth)*1.5), at = (median(depth)*1.5), col = "red")

abline(v = (median(depth)*0.5), col = "red", lwd = 2)
mtext((median(depth)*0.5), at = (median(depth)*0.5), col = "red")

abline(v = (median(depth)*0.75), col = "blue", lwd = 2)
mtext((median(depth)*0.75), at = (median(depth)*0.75), col = "blue")

abline(v = (median(depth)*1.25), col = "blue", lwd = 2)
mtext((median(depth)*1.25), at = (median(depth)*1.25), col = "blue")

abline(v = (median(depth)*2), col = "green", lwd = 2)
mtext((median(depth)*2), at = (median(depth)*2), col = "green")

legend("topright", c("+/- 0.5", "+/- 0.25", "x2"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

# zoom 2: ----------------------------------------------------------
hist(depth, main = NULL, xlim = c(0, 10000), xlab = colnames(df)[1], breaks = 100)
abline(v = median(depth), col = "black", lwd = 2)
mtext(median(depth), at = median(depth))

abline(v = (median(depth)*1.5), col = "red", lwd = 2)
mtext((median(depth)*1.5), at = (median(depth)*1.5), col = "red")

abline(v = (median(depth)*0.5), col = "red", lwd = 2)
mtext((median(depth)*0.5), at = (median(depth)*0.5), col = "red")

abline(v = (median(depth)*0.75), col = "blue", lwd = 2)
mtext((median(depth)*0.75), at = (median(depth)*0.75), col = "blue")

abline(v = (median(depth)*1.25), col = "blue", lwd = 2)
mtext((median(depth)*1.25), at = (median(depth)*1.25), col = "blue")

abline(v = (median(depth)*2), col = "green", lwd = 2)
mtext((median(depth)*2), at = (median(depth)*2), col = "green")

legend("topright", c("+/- 0.5", "+/- 0.25", "x2"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

hist(depth, main = NULL, xlim = c(0, 10000), xlab = colnames(df)[1], breaks = 1000)
abline(v = median(depth), col = "black", lwd = 2)
mtext(median(depth), at = median(depth))

abline(v = (median(depth)*1.5), col = "red", lwd = 2)
mtext((median(depth)*1.5), at = (median(depth)*1.5), col = "red")

abline(v = (median(depth)*0.5), col = "red", lwd = 2)
mtext((median(depth)*0.5), at = (median(depth)*0.5), col = "red")

abline(v = (median(depth)*0.75), col = "blue", lwd = 2)
mtext((median(depth)*0.75), at = (median(depth)*0.75), col = "blue")

abline(v = (median(depth)*1.25), col = "blue", lwd = 2)
mtext((median(depth)*1.25), at = (median(depth)*1.25), col = "blue")

abline(v = (median(depth)*2), col = "green", lwd = 2)
mtext((median(depth)*2), at = (median(depth)*2), col = "green")

legend("topright", c("+/- 0.5", "+/- 0.25", "x2"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

hist(depth, main = NULL, xlim = c(0, 10000), xlab = colnames(df)[1], breaks = 10000)
abline(v = median(depth), col = "black", lwd = 2)
mtext(median(depth), at = median(depth))

abline(v = (median(depth)*1.5), col = "red", lwd = 2)
mtext((median(depth)*1.5), at = (median(depth)*1.5), col = "red")

abline(v = (median(depth)*0.5), col = "red", lwd = 2)
mtext((median(depth)*0.5), at = (median(depth)*0.5), col = "red")

abline(v = (median(depth)*0.75), col = "blue", lwd = 2)
mtext((median(depth)*0.75), at = (median(depth)*0.75), col = "blue")

abline(v = (median(depth)*1.25), col = "blue", lwd = 2)
mtext((median(depth)*1.25), at = (median(depth)*1.25), col = "blue")

abline(v = (median(depth)*2), col = "green", lwd = 2)
mtext((median(depth)*2), at = (median(depth)*2), col = "green")

legend("topright", c("+/- 0.5", "+/- 0.25", "x2"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

# zoom 3: ----------------------------------------------------------
hist(depth, main = NULL, xlim = c(0, 6000), xlab = colnames(df)[1], breaks = 100)
abline(v = median(depth), col = "black", lwd = 2)
mtext(median(depth), at = median(depth))

abline(v = (median(depth)*1.5), col = "red", lwd = 2)
mtext((median(depth)*1.5), at = (median(depth)*1.5), col = "red")

abline(v = (median(depth)*0.5), col = "red", lwd = 2)
mtext((median(depth)*0.5), at = (median(depth)*0.5), col = "red")

abline(v = (median(depth)*0.75), col = "blue", lwd = 2)
mtext((median(depth)*0.75), at = (median(depth)*0.75), col = "blue")

abline(v = (median(depth)*1.25), col = "blue", lwd = 2)
mtext((median(depth)*1.25), at = (median(depth)*1.25), col = "blue")

abline(v = (median(depth)*2), col = "green", lwd = 2)
mtext((median(depth)*2), at = (median(depth)*2), col = "green")

legend("topright", c("+/- 0.5", "+/- 0.25", "x2"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

hist(depth, main = NULL, xlim = c(0, 6000), xlab = colnames(df)[1], breaks = 1000)
abline(v = median(depth), col = "black", lwd = 2)
mtext(median(depth), at = median(depth))

abline(v = (median(depth)*1.5), col = "red", lwd = 2)
mtext((median(depth)*1.5), at = (median(depth)*1.5), col = "red")

abline(v = (median(depth)*0.5), col = "red", lwd = 2)
mtext((median(depth)*0.5), at = (median(depth)*0.5), col = "red")

abline(v = (median(depth)*0.75), col = "blue", lwd = 2)
mtext((median(depth)*0.75), at = (median(depth)*0.75), col = "blue")

abline(v = (median(depth)*1.25), col = "blue", lwd = 2)
mtext((median(depth)*1.25), at = (median(depth)*1.25), col = "blue")

abline(v = (median(depth)*2), col = "green", lwd = 2)
mtext((median(depth)*2), at = (median(depth)*2), col = "green")

legend("topright", c("+/- 0.5", "+/- 0.25", "x2"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

hist(depth, main = NULL, xlim = c(0, 6000), xlab = colnames(df)[1], breaks = 10000)
abline(v = median(depth), col = "black", lwd = 2)
mtext(median(depth), at = median(depth))

abline(v = (median(depth)*1.5), col = "red", lwd = 2)
mtext((median(depth)*1.5), at = (median(depth)*1.5), col = "red")

abline(v = (median(depth)*0.5), col = "red", lwd = 2)
mtext((median(depth)*0.5), at = (median(depth)*0.5), col = "red")

abline(v = (median(depth)*0.75), col = "blue", lwd = 2)
mtext((median(depth)*0.75), at = (median(depth)*0.75), col = "blue")

abline(v = (median(depth)*1.25), col = "blue", lwd = 2)
mtext((median(depth)*1.25), at = (median(depth)*1.25), col = "blue")

abline(v = (median(depth)*2), col = "green", lwd = 2)
mtext((median(depth)*2), at = (median(depth)*2), col = "green")

legend("topright", c("+/- 0.5", "+/- 0.25", "x2"), col = c("red", "blue", "green"), lty = 1, lwd = 2)

dev.off()









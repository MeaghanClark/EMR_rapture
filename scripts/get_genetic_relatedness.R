library(adegenet)
# install dev version of dartRverse 
# devtools::install_github("green-striped-gecko/dartR.base@dev")
library(dartR.base)
library(SNPRelate)

### load genotype data 
# import filtered genotype data from filter_vcf_popgen.R

load(file = "../pedigree_reconstruction/PCC_filt_gt_ped.Robj")
load(file = "../pedigree_reconstruction/ELF_filt_gt_ped.Robj")

# drop ELF_269
ELF_filt_gt_ped <- ELF_filt_gt_ped[,-which(colnames(ELF_filt_gt_ped) == "ELF_269")]

### convert to SNPRelate format

  # convert to genlight object
ELF_genlight <- new("genlight", gen = t(ELF_filt_gt_ped), ind.names = colnames(ELF_filt_gt_ped))
PCC_genlight <- new("genlight", gen = t(PCC_filt_gt_ped), ind.names = colnames(PCC_filt_gt_ped))

  # add dummy information to loc.all slot of genlight objects
    # https://groups.google.com/g/dartr/c/ZGOOKlf2mUs/m/P2PcHyWOBgAJ

ELF_genlight$loc.all <- rep("G/C",nLoc(ELF_genlight))
PCC_genlight$loc.all <- rep("G/C",nLoc(PCC_genlight))

  # add population ID as a factor to genlight object
ELF_genlight@pop <- as.factor(rep("ELF",dim(ELF_genlight)[1]))
PCC_genlight@pop <- as.factor(rep("PCC",dim(PCC_genlight)[1]))

  # convert genlight object to plink
    # must download binary file of plink 1.9 per help(gl2plink) details 
gl2plink(x = ELF_genlight, outfile = "ELF_plink", outpath = "../vcf_stats/", bed.files = TRUE, 
         plink.bin.path = "/Users/meaghan/Desktop/EMR_rapture/scripts/plink_mac_20241022")

gl2plink(x = PCC_genlight, outfile = "PCC_plink", outpath = "../vcf_stats/", bed.files = TRUE, 
         plink.bin.path = "/Users/meaghan/Desktop/EMR_rapture/scripts/plink_mac_20241022")

  # covert plink to snprelate using snpgdsBED2GDS

ELF.bed.fn <- "/Users/meaghan/Desktop/EMR_rapture/vcf_stats/ELF_plink.bed"
ELF.fam.fn <- "/Users/meaghan/Desktop/EMR_rapture/vcf_stats/ELF_plink.fam"
ELF.bim.fn <- "/Users/meaghan/Desktop/EMR_rapture/vcf_stats/ELF_plink.bim"
snpgdsBED2GDS(ELF.bed.fn, ELF.fam.fn, ELF.bim.fn, "/Users/meaghan/Desktop/EMR_rapture/vcf_stats/ELF.gds")
snpgdsSummary("/Users/meaghan/Desktop/EMR_rapture/vcf_stats/ELF.gds")

PCC.bed.fn <- "/Users/meaghan/Desktop/EMR_rapture/vcf_stats/PCC_plink.bed"
PCC.fam.fn <- "/Users/meaghan/Desktop/EMR_rapture/vcf_stats/PCC_plink.fam"
PCC.bim.fn <- "/Users/meaghan/Desktop/EMR_rapture/vcf_stats/PCC_plink.bim"
snpgdsBED2GDS(PCC.bed.fn, PCC.fam.fn, PCC.bim.fn, "/Users/meaghan/Desktop/EMR_rapture/vcf_stats/PCC.gds")
snpgdsSummary("/Users/meaghan/Desktop/EMR_rapture/vcf_stats/PCC.gds")


### calculate relatedness with snpgdsIBDKING

# ELF
# open genofile 
(genofile <- snpgdsOpen("/Users/meaghan/Desktop/EMR_rapture/vcf_stats/ELF.gds"))

ibd.robust <- snpgdsIBDKING(genofile)

ELF_ibd <- snpgdsIBDSelection(ibd.robust)
head(ELF_ibd)

# close genofile 
closefn.gds(genofile)

# PCC

(genofile <- snpgdsOpen("/Users/meaghan/Desktop/EMR_rapture/vcf_stats/PCC.gds"))

ibd.robust <- snpgdsIBDKING(genofile)

PCC_ibd <- snpgdsIBDSelection(ibd.robust)
head(PCC_ibd)

closefn.gds(genofile)

#  save kinship data frames

save(ELF_ibd, file = "/Users/meaghan/Desktop/EMR_rapture/vcf_stats/ELF_ibd.Robj")
save(PCC_ibd, file = "/Users/meaghan/Desktop/EMR_rapture/vcf_stats/PCC_ibd.Robj")


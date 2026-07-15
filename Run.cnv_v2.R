# Detecting and filtering copy number variants in radseq data

library(rCNV)
library(adegenet)

# homedir
setwd("/Volumes/CatDisk/SDZooRAD/2024/rCNV")

## Get attributes of samples  ---------------------------------------------------------------
att <- read.table("/Volumes/CatDisk/SDZooRAD/2024/Rworking/elatAttributes.txt", header = T, sep = "\t", colClasses = "character", quote = "")


# import vcf data from main stacks output - no filters
vcf.file <- "/Volumes/CatDisk/SDZooRAD/2024/stacks/noReps/mac3R7r7nuc/populations.snps.vcf"
vcf <- readVCF(vcf.file)

# genpop data to compare loci
h <- "/Volumes/CatDisk/SDZooRAD/2024/Rworking/data/allSamplesNeutHaps.gen"

haps <- read.genepop(h)

# calculate missing data
missing.data <- get.miss(vcf)

miss.samples <- get.miss(vcf, type = "samples")
miss.sam.filter <- which(miss.samples$perSample$f_miss>0.5)+9 # column indices
length(miss.sam.filter) # individuals removed
# 30

# filter samples
vcf.samp.filt <- data.frame(vcf$vcf)[,-miss.sam.filter]
# create new object
vcf.50 <- list(vcf=vcf.samp.filt)

# relatedness filters
rel <- relatedness(vcf.50)
which(rel$relatedness_Ajk>0.9)
closeKin <- which(rel$relatedness_Ajk>0.9 & rel$indv1!=rel$indv2)
relFilt <- unique(rel[closeKin,1])
relFiltInd <- which(names(vcf.50$vcf) %in% relFilt)

names(vcf.50$vcf)[relFiltInd]

vcf.rel <- subset(data.frame(vcf.50$vcf), select = -relFiltInd)

vcf.filt <- list(vcf=vcf.rel)

# Allele depths table
ad.tab <- hetTgen(vcf.50, info.type="AD")

# correct GT mismatches
gt.tab <- hetTgen(vcf.50, info.type = "GT")
ad.tab <- ad.correct(ad.tab, gt.table = gt.tab)

# Normalized allele depths table
ad.norm <- cpm.normal(ad.tab, method = "MedR")
#OUTLIERS DETECTED
# Consider removing the samples:
#   Elat112 Elat117

# remove outliers
ad.norm$AD <- ad.norm$AD[,-ad.norm$outliers$column]

# Calc h and Fis
hz <- h.zygosity(vcf.filt)
Fis <- mean(hz$Fis, na.rm = T)

# allele.info
a.info <- allele.info(X=ad.tab, x.norm = ad.norm, Fis = Fis)

# detect deviants
deviants <- dupGet(a.info, Fis = Fis, test = c("z.05", "chi.05"), plot = T, )

# flag CNVs
dups.int <- cnv(a.info, test = c("z.05", "chi.05"), filter = "intersection")

dups.kmeans <- cnv(a.info, test = c("z.05", "chi.05"), filter = "kmeans")

table(deviants$dup.stat)
table(dups.int$dup.stat)
table(dups.kmeans$dup.stat)

# get dups locus-wise
dup.loci.all <- dups.kmeans[which(dups.kmeans$dup.stat=="cnv"),] # all cnv snps
dup.loci <- unique(dup.loci.all$CHROM) # all cnv loci
length(dup.loci)

# Number of paralogs in previously filtered loci
sum(dup.loci %in% locNames(haps))

# save tables
write.csv(x=ad.tab, file = "ad.tab.50.corrected.csv")
write.csv(x=ad.norm$AD, file = "ad.norm.50.cor.noOut.csv")
write.csv(x=gt.tab, file = "gt.tab.csv")
write.csv(x=a.info, file = "a.info.csv")

# save paralogs to file
write(dup.loci, file = "paralogs.txt", ncolumns = 1)



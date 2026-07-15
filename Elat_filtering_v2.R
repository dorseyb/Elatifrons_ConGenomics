
# Description -------------------------------------------------------------
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # # ## ##
# R script to filter data after process_radtags in Stacks found in Dorsey et al 2023 Conservation genomics of Dioon holmgrenii (Zamiaceae) reveals
# a history of range expansion, fragmentation, and isolation
# of populations. Conservation Genetics
#
# 9/9/24
# Modified for use in Conservation Genomics Provides Essential Insights into the Provenance of Ex Situ Collections and their Utility
# for Species Recovery: A Case Study in Encephalartos latifrons
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # # ## ##

# Libraries ------------------------------------------------------------------------------------------------------------------------------

setwd("/Volumes/CatDisk/SDZooRAD/2024/Rworking/Revision")
library(adegenet)
library(hierfstat)
library(qvalue)
library(pcadapt)


# Functions ------------------------------------------------------------------------------------------------------------------------------

# write Structure files
source("/Users/bdorsey/Library/CloudStorage/OneDrive-TheHuntingtonLibrary,ArtMuseum,andBotanicalGardens/Documents/Rfunctions/genind2structure.R")

#  From https://rdrr.io/github/romunov/zvau/src/R/writeGenPop.R
writeGenPop <- function(gi, file.name, comment) {
  
  if (is.list(gi)) {
    # do all genind objects have the same number of loci?
    if (length(unique(sapply(gi, nLoc))) != 1) stop("Number of loci per individual genind object in a list is not equal for all.")
    gi.char <- gi[[1]]
    loc.names <- locNames(gi[[1]])
  } else {
    gi.char <- gi
    loc.names <- locNames(gi)
  }
  
  # Calculate the length of two alleles.
  for (ln in 1:length(locNames(gi.char))){
      lng <- as.character(na.omit(genind2df(gi.char)[, locNames(gi.char)[ln]]))
      if(length(lng)==0){
        next
      }
      else{
        lng <- unique(nchar(lng))
        stopifnot(length(lng) == 1)
        break
      }
  }
  
  cat(paste(comment, "\n"), file = file.name)
  cat(paste(paste(loc.names, collapse = ", "), "\n"), file = file.name, append = TRUE)
  
  if (is.list(gi)) {
    pop.names <- seq_len(length(gi))
  } else {
    pop.names <- popNames(gi)
  }
  
  for (i in pop.names) {
    cat("pop\n", file = file.name, append = TRUE)
    if (is.list(gi)) {
      intm <- gi[[i]]
      loc.names <- locNames(gi[[i]])
    } else {
      intm <- gi[pop(gi) == i, drop = FALSE]
    }
    ind.names <- indNames(intm)
    intm <- genind2df(intm, sep = "")
    intm[is.na(intm)] <- paste(rep("0", lng), collapse = "")
    out <- cbind(names = paste(ind.names, ",", sep = ""), intm[, loc.names])
    write.table(out, file = file.name, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
  }
  
  return(NULL)
}

# Add attributes to genind object - requires csv file with one column ID matching 
# individual names in genind object and other attributes in subsequent columns
addAtts <- function(attfile, genOb){
  ## Get attributes of samples
  # att <- read.table("elatAttributes.txt", header = T, sep = "\t", colClasses = "character", quote = "")
  att <- read.table(attfile, header = T, sep = "\t", colClasses = "character", quote = "")
  # add row names for sorting
  row.names(att) <- att$ID
  # order by individual code in genind obj
  # att <- att[row.names(allSampleNeutHaps@tab),]
  att <- att[row.names(genOb@tab),]
  # check order
  # cbind(att$ID,row.names(genOb@tab))
  
  # add as Strata and print Pop to status (Garden, wild, planted)
  strata(genOb) <- att

  return(genOb)
}

# Directories and Files ------------------------------------------------------------------------------------------------------------------

# working dir for paths
dir <- getwd()

# Stacks directory with output from populations program - genepop and haplotype files
stackdir <- '/Volumes/CatDisk/SDZooRAD/2024/stacks/noReps'

# Import data from Stacks -----------------------------------------------------

#nuclear Loci
nuc.haps.genepop <- paste(stackdir, "mac3R7r7nuc/populations.haps.genepop.gen", sep='/')
nuc.haps.tsv <- paste(stackdir, "mac3R7r7nuc/populations.haplotypes.tsv", sep='/')
nuc.snps.genepop <- paste(stackdir, "mac3R7r7nuc/populations.snps.genepop.gen", sep='/')

## Read in files----------
# haplotypes in genepop format
haps <- read.genepop(file=nuc.haps.genepop)

# snps in genepop format - needed for BayeScan
snps <- read.genepop(file=nuc.snps.genepop)

# raw haplotypes for fineRAD
haps.raw <- read.delim(nuc.haps.tsv, check.names = F)


# set attributes
attFile <- "elatAttributes.txt"
haps <- addAtts(attFile, haps)
snps <- addAtts(attFile, snps)

setPop(haps) <- ~popAbs
setPop(snps) <- ~popAbs

# Count minor alleles ---------------------------------------------------------------------------
mafHaps <- minorAllele(haps)
maf05 <- which(mafHaps < 0.05)
length(maf05)/length(mafHaps)

maf01 <- which(mafHaps < 0.01)
length(maf01)/length(mafHaps)

mafSnps <- minorAllele(snps)
maf05Snps <- which(mafSnps < 0.05)
length(maf05Snps)/length(mafSnps)

maf01Snps <- which(mafSnps < 0.01)
length(maf01Snps)/length(mafSnps)


# pcadapt for non-neutral loci----------------------------------------------------------------------------------------------------

# use plink to convert vcf to bed

# read in bed file for pcadapt
bedFile <- "data/elat.wild.70.snps.bed"
snps70.pcadapt <- read.pcadapt(bedFile, type = "bed")

# run with many pca axes to decide how many to use
ktest1 <- pcadapt(snps70.pcadapt, K=20, min.maf = 0)

# check screeplot to find optimal # pc's
plot(ktest1, option = "screeplot")

# score plot
plot(ktest1, option = "scores")
plot(ktest1, option = "scores", i=3, j=4)
plot(ktest1, option = "scores", i=5, j=6)
plot(ktest1, option = "scores", i=7, j=8)
plot(ktest1, option = "scores", i=9, j=10)

# Choose K=3 for nuc and for org

# set min maf
maf <- 0

# run with K=3
adapt.K3.res <- pcadapt(snps70.pcadapt, K=3, min.maf = maf)
plot(adapt.K3.res, option = "manhattan")
plot(adapt.K3.res, option = "qqplot")
hist(adapt.K3.res$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

# choose outliers using qvalues
qval <- qvalue(adapt.K3.res$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha) # indices of outliers
length(outliers)
# 807 nuc

# select outlier snps
snpsOuts <- snps[,outliers]

# Get outlier names
snpOutLociNames <- locNames(snps[,outliers])
# Remove allele numbers
snpOutLociNames <- gsub("_[0-9]+", "", x=snpOutLociNames)
# Get unique loci
snpOutLociNames <- sort(unique(snpOutLociNames))
# 463 loci

# get names of neutral loci
NeutLoci <- setdiff(locNames(haps), snpOutLociNames)
# 2820 loci

# all samples, neutral loci
allSampleNeutHaps <- haps[, loc = NeutLoci]

# Convert to df
allSamplesNeutHaps.df <- genind2df(allSampleNeutHaps, oneColPerAll = T)
allSamplesNeutHaps.df <- cbind(row.names(allSamplesNeutHaps.df), allSamplesNeutHaps.df)
names(allSamplesNeutHaps.df)[1:2] <- c("Individual", "Population")

# haps under selection
allSampleSelHaps <- haps[,loc = snpOutLociNames]

# loci under selection
allSampleSelSnps <- snps[,outliers]

write.table(NeutLoci, file="neutralLoci20260405.txt")

# Filter paralogs/CNV using rCNV -----------------------------------------------

## Get dups with Run.cnv.R script first-----------------------------------------

## import list of dups, convert to character vector-----------------------------
dups <- read.table("/Volumes/CatDisk/SDZooRAD/2024/rCNV/paralogs.txt")
dups <- as.character(dups$V1)
length(dups)
# 285

## get indices of NON dup names = loci to keep----------------------------------
nodupInd <- which(!(locNames(allSampleNeutHaps) %in% dups))
length(nodupInd)
length(locNames(allSampleNeutHaps)) - length(nodupInd)
# 232

## remove dup loci-----------------------------------------------
allSampleNeutHaps <- allSampleNeutHaps[,loc = nodupInd]


# Calculate missing data ------------------------------------------------------------------------------------------------------------------

indTyped <- propTyped(allSampleNeutHaps, by="ind")
locTyped <- propTyped(allSampleNeutHaps, by="loc")

m <- hist(indTyped)
m <- hist(1-indTyped)

# make a table of loci presence and save
attOrder <- match(names(indTyped), att$ID)
indTypedID <- cbind(att$RubyID[attOrder], indTyped)
colnames(indTypedID) <- c("RubyID","Loci Present")
write.csv(indTypedID, "PercentLociPresent.csv", quote = F)

# remove 2 very low data individuals
lt0.1 <- which(indTyped<0.1)
allSampleNeutHaps <- allSampleNeutHaps[-(lt0.1),]

n <- hist(locTyped)


# Testing effect of missing data ----------------------------------------------------------------------------------------------------

# Ind per "pop" pre missing filter
table(allSampleNeutHaps@pop)
# seed KNBG   WM   GH   TZ    R  BEL RBGA SDBG  MBC   CH  SDZ   UC    L   JC  HBG   FS   TR   ST   VK    B 
#    6   25   15   17   10    5    1    2    2    4   10    2    3    7    4   10   13   20    3    7    1 

preMissingHaps <- allSampleNeutHaps

# previous filtering
# loc70ind10 <- preMissingHaps[,loc=which(locTyped > 0.7)]
# loc70ind10 <- loc70ind10[which(indTyped > 0.1),]


# filter loci to 80% occupancy first and then individuals
length(which(locTyped < 0.8))
# 1650 removed
loc80haps <- preMissingHaps[,loc=which(locTyped >= 0.8)]
setPop(loc80haps) <- ~popAbs

table(loc80haps@pop)
# seed KNBG   WM   GH   TZ    R  BEL RBGA SDBG  MBC   CH  SDZ   UC    L   JC  HBG   FS   TR   ST   VK    B 
#    6   25   15   17   10    5    1    2    2    4   10    2    3    7    4   10   13   20    3    7    1 

loc80IndTyped <- propTyped(loc80haps, by="ind")
l <- hist(loc80IndTyped)

# Individuals > 20%
loc80Ind20haps <- loc80haps[which(loc80IndTyped >= 0.2), , drop=T]
popNames(loc80Ind20haps)
setPop(loc80Ind20haps) <- ~popAbs
table(loc80Ind20haps@pop)
# seed KNBG   WM   GH   TZ    R  BEL RBGA SDBG  MBC   CH  SDZ   UC    L   JC  HBG   FS   TR   ST   VK    B 
#    5   24   15   16   10    5    1    2    2    4   10    2    3    7    4   10   13   20    3    7    1

# Individuals > 50%
length(which(loc80IndTyped < 0.5))
loc80Ind50haps <- loc80haps[which(loc80IndTyped >= 0.5), , drop=T]
popNames(loc80Ind50haps)
setPop(loc80Ind50haps) <- ~popAbs
table(loc80Ind50haps@pop)

indNames(loc80Ind50haps)


# seed KNBG   WM   GH   TZ    R  BEL RBGA SDBG  MBC   CH  SDZ   UC    L   JC  HBG   FS   TR   ST   VK    B 
#    1   24   15   12   10    5    1    2    2    4    3    2    3    7    4   10   13   20    3    7    1 

# Individuals > 70%
loc80Ind70haps <- loc80haps[which(loc80IndTyped >= 0.7), , drop=T]
popNames(loc80Ind70haps)
setPop(loc80Ind70haps) <- ~popAbs
table(loc80Ind70haps@pop)

# seed KNBG   WM   GH   TZ    R  BEL RBGA SDBG  MBC   CH  SDZ   UC    L   JC  HBG   FS   TR   VK    B 
#    1   22   15    9   10    5    1    2    2    4    2    2    3    7    4   10   12   15    7    1 


# number of loci
length(which(locTyped >=0.8))
# 938

# loci and individuals
finalLoci <- locNames(loc80Ind70haps)
Ind20 <- indNames(loc80Ind20haps)
Ind50 <- indNames(loc80Ind50haps)

# Separate wild samples --------------------------------------------------------

# wild samples
setPop(loc80Ind20haps) <- ~Status
hapsInd20Wild <- loc80Ind20haps[pop = c("Wild", "planted"), drop = T]

setPop(loc80Ind50haps) <- ~Status
hapsInd50WildPL <- loc80Ind50haps[pop = c("Wild", "planted"), drop = T]
hapsInd50Wild <- loc80Ind50haps[pop = c("Wild"), drop = T]
indNames(hapsInd50Wild)
pop(loc80Ind50haps)["Elat63"]

# filter Raw ------------------------------------------------------------------

# get the indices of filtered loci in raw dataframe column 1
rawFinalLoci <- which(haps.raw[,1] %in% locNames(hapsInd50Wild))

# Ind > 20%
hapsNeutRawInd20 <- haps.raw[rawFinalLoci,Ind20]
# replace '-' with '' for missing data
hapsNeutRawInd20[hapsNeutRawInd20 == '-'] <- ''

# get wild samples only
finalWildInd20 <- indNames(hapsInd20Wild)
hapsInd20WildRaw <- haps.raw[rawFinalLoci,finalWildInd20]
# replace '-' with '' for missing data
hapsInd20WildRaw[hapsInd20WildRaw == '-'] <- ''
write.table(finalWildInd20, "wildsamplesInd20.txt", quote = F, row.names = F, col.names = F)

# Ind > 50%
hapsNeutRawInd50 <- haps.raw[rawFinalLoci,Ind50]
# replace '-' with '' for missing data
hapsNeutRawInd50[hapsNeutRawInd50 == '-'] <- ''
colnames(hapsNeutRawInd50)

# get wild samples only
finalWildInd50 <- indNames(hapsInd50Wild)
hapsInd50WildRaw <- haps.raw[rawFinalLoci,finalWildInd50]
# replace '-' with '' for missing data
hapsInd50WildRaw[hapsInd50WildRaw == '-'] <- ''

write.table(finalWildInd50, "wildsamplesInd50.txt", quote = F, row.names = F, col.names = F)

# write files ------------------------------------------------------------------------------------------------------------------

## nuc loci--------------------------------------------------------------------------------

# print filtered loci names to file
write(finalLoci, "filt80Loci.txt")

# Genepop files
writeGenPop(allSampleNeutHaps, paste(dir,"data/allSamplesNeutHaps.gen",sep = "/"), "Neutral loci>70, inds>10")

writeGenPop(loc80Ind20haps, paste(dir,"data/loc80Ind20haps.gen",sep = "/"), "Neutral loci>80, inds>20")
writeGenPop(loc80Ind50haps, paste(dir,"data/loc80Ind50haps.gen",sep = "/"), "Neutral loci>80, inds>50")
writeGenPop(loc80Ind70haps, paste(dir,"data/loc80Ind70haps.gen",sep = "/"), "Neutral loci>80, inds>70")

writeGenPop(hapsInd50Wild, paste(dir,"data/loc80Ind50Wildhaps.gen",sep = "/"), "Neutral loci>80, inds>50, wild only")

# write non-paralog loci for phylo
write.table(locNames(loc80Ind50haps), file = "/Volumes/CatDisk/SDZooRAD/2024/phylo/Revision/analyses/iqtree/nonParalogs.txt")

# Write files for fineRADstructure
write.table(hapsNeutRawInd20, file = "hapsInd20_frs.txt", quote = F, sep = "\t", na = "", row.names = F)
write.table(hapsInd20WildRaw, file = "hapsWild20Ind_frs.txt", quote = F, sep = "\t", na = "", row.names = F)

write.table(hapsNeutRawInd50, file = "hapsInd50_frs.txt", quote = F, sep = "\t", na = "", row.names = F)
write.table(hapsInd50WildRaw, file = "hapsWild50Ind_frs.txt", quote = F, sep = "\t", na = "", row.names = F)

# structure files

setPop(loc80Ind50haps) <- ~popAbs
popNames(loc80Ind50haps)

setPop(hapsInd50Wild) <- ~popAbs
popNames(hapsInd50Wild)

genind2structure(loc80Ind50haps, "/Volumes/CatDisk/SDZooRAD/2024/structure/revision/allSamples/allSamplesFilt.structure", pops = T)
genind2structure(hapsInd50Wild, "/Volumes/CatDisk/SDZooRAD/2024/structure/revision/wildSamples/hapsWildNeut.structure", pops = T)

# Use WL to write snp files in Stacks (populations)


## snps for BayesAss------------------------------------------------------------

# file of filtered snps from WL
s <- '/Volumes/CatDisk/SDZooRAD/2024/stacks/noReps/revision/randSnpsWL/ind50/populations.snps.p.snps.gen'

snps <- read.genepop(s)

snps <- addAtts(attFile, snps)

setPop(snps) <- ~popAbs
popNames(snps)

snpsWildBAPops <- snps[pop = c("GH","TR","VK","CH","FS"), drop = T]
popNames(snpsWildBAPops)

# write files
writeGenPop(snpsWildBAPops, "data/snpsWildBApops.gen", "random snps filtered")


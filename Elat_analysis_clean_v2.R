############################################################################################################################################################
# R script to  run analyses found in Dorsey et al Conservation Genomics Provides Essential Insights into the Provenance of Ex Situ Collections and their Utility
# for Species Recovery: A Case Study in Encephalartos latifrons
# 
# Best to run by sections in Rstudio rather than calling the script from the terminal.
# Some options/directories need to be set within sections
#############################################################################################################################################################



# Libraries ---------------------------------------------------------------

{ # this block loads all libraries and functions
  
  setwd("/Volumes/CatDisk/SDZooRAD/2024/Rworking/Revision")
  # devtools::install_github("thierrygosselin/radiator")
  # BiocManager::install("qvalue")
  # install.packages("fields")
  library(adegenet)
  library(ggplot2)  
  library(ggrepel)
  library(gridExtra)
  library(ggnewscale)
  library(fields)
  library(plyr)
  library(rgl)
  library(mmod)
  library(hierfstat)
  library(parallel)
  library(boot)
  library(poppr)
  library(vegan)
  library(pairwiseAdonis)
  
  # default params----
  par.defaults <- par(no.readonly = T)
  par.defaults$mar
  # Functions ---------------------------------------------------------------
  
  source("/Users/bdorsey/Library/CloudStorage/OneDrive-TheHuntingtonLibrary,ArtMuseum,andBotanicalGardens/Documents/Rfunctions/genind2structure.R")
  
  # get Hoban functions
  source("/Volumes/CatDisk/SDZooRAD/2024/Rworking/Hoban/Fa_sample_funcs.R")
  
  # Add attributes
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
    # popNames(genOb)
    # setPop(allSampleNeutHaps) <- ~popAbs
    return(genOb)
  }
  
  
  
  makePlotData <- function(gInd=NULL, dapcObj=NULL) {
    numAx <- dapcObj$n.da
    co <- as.data.frame(dapcObj$ind.coord)
    stratMatch <- match(rownames(dapcObj$ind.coord), gInd$strata$ID)
    co$pop <- gInd$strata$popAbs[stratMatch]
    co$status <- gInd$strata$Status[stratMatch]
    co$frs <- gInd$strata$frsGroup[stratMatch]
    return(co)
  }
  
  newAssignplot <- function (x, only.grp = NULL, subset = NULL, new.pred = NULL, cex.lab = 0.75, pch = 3)
  {
    require(fields)
    if (!inherits(x, "dapc")) 
      stop("x is not a dapc object")
    if (!is.null(new.pred)) {
      n.new <- length(new.pred$assign)
      x$grp <- c(as.character(x$grp), rep("unknown", n.new))
      x$assign <- c(as.character(x$assign), as.character(new.pred$assign))
      x$posterior <- rbind(x$posterior, new.pred$posterior)
    }
    if (!is.null(only.grp)) {
      only.grp <- as.character(only.grp)
      ori.grp <- as.character(x$grp)
      x$grp <- x$grp[only.grp == ori.grp]
      x$assign <- x$assign[only.grp == ori.grp]
      x$posterior <- x$posterior[only.grp == ori.grp, , drop = FALSE]
    }
    else if (!is.null(subset)) {
      x$grp <- x$grp[subset]
      x$assign <- x$assign[subset]
      x$posterior <- x$posterior[subset, , drop = FALSE]
    }
    n.grp <- ncol(x$posterior)
    n.ind <- nrow(x$posterior)
    Z <- t(x$posterior)
    Z <- Z[, ncol(Z):1, drop = FALSE]
    image.plot(x = 1:n.grp, y = seq(0.5, by = 1, le = n.ind), Z, 
               col = rev(heat.colors(100)), yaxt = "n", ylab = "", 
               xaxt = "n", xlab = "Clusters")
    axis(side = 1, at = 1:n.grp, tick = FALSE, labels = colnames(x$posterior))
    # axis(side = 2, at = seq(0.5, by = 1, le = n.ind), labels = rev(rownames(x$posterior)), 
    # las = 1, cex.axis = cex.lab)
    abline(h = 1:n.ind, col = "lightgrey")
    abline(v = seq(0.5, by = 1, le = n.grp))
    box()
    newGrp <- colnames(x$posterior)
    x.real.coord <- rev(match(x$grp, newGrp))
    y.real.coord <- seq(0.5, by = 1, le = n.ind)
    points(x.real.coord, y.real.coord, col = "deepskyblue2", 
           pch = pch)
    return(invisible(match.call()))
  }
  
  medianFunc <- function(data, indices) {
    data <- as.matrix(data)
    dt <- data[indices,]
    med <- apply(dt,2,median, na.rm = T)
    return(med)
  }
  
  medianFunc2 <- function(data, indices) {
    data <- as.matrix(data)
    dt <- data[indices]
    med <- median(dt, na.rm = T)
    return(med)
  }
  
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
    lng <- as.character(na.omit(genind2df(gi.char)[, locNames(gi.char)[1]]))
    lng <- unique(nchar(lng))
    
    stopifnot(length(lng) == 1)
    
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
  
  # plotting function
  plotPC <- function (genind=NULL,pca=NULL,axes=c(1,2,3,4)) {
    pcCo <- as.data.frame(pca$x)
    pcCo$pop <- genind$pop
    pcCo$tax <- genind@strata$Taxon
    pcCo$source <- genind$strata$Status
    pcCo$loc <- genind$strata$Local_Gen
    pcCo$locS <- genind@strata$Local_sp
    pcCo$comp <- genind@strata$comLab
    pcCo$ab <- genind@strata$popAbs
    axNum <- length(axes)
    
    pc12.loc <- ggplot(pcCo) + 
      geom_point(aes(PC1, PC2, color=loc), size = 5) +
      theme(legend.position = "top") 
    # +
    # geom_text_repel(data=pcCo,  x=pcCo$PC1, y=pcCo$PC2, label=pcCo$loc, max.overlaps = 30, size = 3)
    
    pc34.loc <- ggplot(pcCo) + 
      geom_point(aes(PC3, PC4, color=loc), size = 5) +
      theme(legend.position = "top") 
    # +
    #   geom_text_repel(data=pcCo,  x=pcCo$PC3, y=pcCo$PC4, label=pcCo$loc, max.overlaps = 30, size = 3)
    # 
    pc12.tax <- ggplot(pcCo) +
      geom_point(aes(PC1, PC2, color=tax)) +
      theme(legend.position = "top") +
      geom_text_repel(data=pcCo,  x=pcCo$PC1, y=pcCo$PC2, label=pcCo$loc, max.overlaps = 30, size = 3)
    
    pc34.tax <- ggplot(pcCo) + 
      geom_point(aes(PC3, PC4, color=tax)) +
      theme(legend.position = "top") +
      geom_text_repel(data=pcCo,  x=pcCo$PC3, y=pcCo$PC4, label=pcCo$loc, max.overlaps = 30, size = 3)
    
    pc12.gw <- ggplot(pcCo) + 
      geom_point(aes(PC1, PC2, color=source)) +
      theme(legend.position = "top") +
      geom_text_repel(data=pcCo,  x=pcCo$PC1, y=pcCo$PC2, label=pcCo$locS, max.overlaps = 30, size = 3)
    
    pc34.gw <- ggplot(pcCo) + 
      geom_point(aes(PC3, PC4, color=source)) +
      theme(legend.position = "top") +
      geom_text_repel(data=pcCo,  x=pcCo$PC3, y=pcCo$PC4, label=pcCo$locS, max.overlaps = 30, size = 3)
    
    pc12.locS <- ggplot(pcCo) + 
      geom_point(aes(PC1, PC2, color=loc)) +
      theme(legend.position = "top") +
      geom_text_repel(data=pcCo,  x=pcCo$PC1, y=pcCo$PC2, label=pcCo$ab, max.overlaps = 30, size = 3)
    
    pc34.locS <- ggplot(pcCo) + 
      geom_point(aes(PC3, PC4, color=loc)) +
      theme(legend.position = "top") +
      geom_text_repel(data=pcCo,  x=pcCo$PC3, y=pcCo$PC4, label=pcCo$ab, max.overlaps = 30, size = 3)
    
    pc12.comp <- ggplot(pcCo) + 
      geom_point(aes(PC1, PC2, color=source)) +
      # theme(legend.position = "top") +
      geom_text_repel(data=pcCo,  x=pcCo$PC1, y=pcCo$PC2, label=pcCo$comp, max.overlaps = 30, size = 3)
    
    pc34.comp <- ggplot(pcCo) + 
      geom_point(aes(PC3, PC4, color=source)) +
      # theme(legend.position = "top") +
      geom_text_repel(data=pcCo,  x=pcCo$PC3, y=pcCo$PC4, label=pcCo$comp, max.overlaps = 30, size = 3)
    
    return((list(loc12=pc12.loc, tax12=pc12.tax, source12=pc12.gw,loc34=pc34.loc, 
                 tax34=pc34.tax, source34=pc34.gw, locS12=pc12.locS, locS34=pc34.locS,
                 comp12=pc12.comp, comp34=pc34.comp)))
  }
  
  # count private alleles and calculate proportions in each population
  writePrivAlleles <- function(giObj,gpObj,hierObj,outfile) {
    priv.alls.dose <- poppr::private_alleles(giObj)
    # priv.alls.dose[1:4, 1:5]
    #total count of private allele copies per population
    private.count <- rowSums(priv.alls.dose)
    
    #count of alleles that are private to each population
    priv.alls <- poppr::private_alleles(giObj,count.alleles = F)
    # priv.alls[1:4,1:10]
    priv.alls.sums <- rowSums(priv.alls)
    
    # number of distinct alleles present in each pop
    distinctallelesInpops <- apply(gpObj@tab,1,function(l){length(l[l != 0])})
    
    # percentage of alleles present in each pop that are private
    perc.priv <- round((priv.alls.sums/distinctallelesInpops*100),2)
    
    # number of individuals per population
    Indcounts <- summary(giObj@pop)
    
    # possible number of genotypes per population
    all.poss <- Indcounts*nLoc(giObj)
    
    # number of missing genotypes per individual
    n.na.gt <- apply(is.na(hierObj),1,sum)
    
    # number of missing genotypes per pop
    pop.na.gt <- as.numeric(t(rowsum(n.na.gt,hierObj$pop)))
    
    # number of present genotypes
    actual.genos <- all.poss-pop.na.gt
    
    # actual number of alleles (copies) per pop
    total.alls <- actual.genos*2
    
    # percentage of copies that are private alleles
    perc.private.copies <- as.numeric(round((private.count/total.alls)*100,2))
    # as proportion (freq)
    freq.private.copies <- as.numeric(round((private.count/total.alls),4))
    names(perc.private.copies) <- names(perc.priv)
    
    # put into table
    private.allele.df <- data.frame(names(private.count),Indcounts,actual.genos,total.alls,distinctallelesInpops,private.count,perc.private.copies,priv.alls.sums,perc.priv)
    
    pa.df.sort <- arrange(private.allele.df,perc.priv)
    # pa.df.sort
    
    write.csv(pa.df.sort, file = outfile)
    
    return(pa.df.sort)
    
  }
  
  # Calculate diversity stats
  getInDiv <- function(giObj, hierObj) {
    basics <- basic.stats(giObj, diploid = T, digits = 3)
    basics$overall
    
    basics.pop.meds <- matrix( 
      c(
        apply(basics$Fis,2,median, na.rm=T),
        apply(basics$Ho,2,median, na.rm=T),
        apply(basics$Hs,2,median, na.rm=T)
      ),
      nrow = 3, ncol = length(popNames(giObj)), byrow = T, dimnames = list(c("Fis", "Ho", "Hs"), popNames(giObj))
    )

    basics.pop.mean <- matrix(
      c(round(apply(basics$Fis,2,mean, na.rm=T),3),
        round(apply(basics$Ho,2,mean, na.rm=T),3),
        round(apply(basics$Hs,2,mean, na.rm=T),3)
      ),
      nrow = 3, ncol = length(popNames(giObj)), byrow = T, dimnames = list(c("Fis", "Ho", "Hs"), popNames(giObj))
    )
    
    basics.pop.sd <- matrix(
      c(round(apply(basics$Fis,2,sd, na.rm=T),3),
        round(apply(basics$Ho,2,sd, na.rm=T),3),
        round(apply(basics$Hs,2,sd, na.rm=T),3)
      ),
      nrow = 3, ncol = length(popNames(giObj)), byrow = T, dimnames = list(c("Fis", "Ho", "Hs"), popNames(giObj))
    )
    
    basics.pop.mad <- matrix(
      c(round(apply(basics$Fis,2,mad, constant=1, na.rm=T),3),
        round(apply(basics$Ho,2,mad, constant=1, na.rm=T),3),
        round(apply(basics$Hs,2,mad, constant=1, na.rm=T),3)
      ),
      nrow = 3, ncol = length(popNames(giObj)), byrow = T, dimnames = list(c("Fis", "Ho", "Hs"), popNames(giObj))
    )
    
    # allelic richnes
    
    arHier <- allelic.richness(hierObj)
    arMeds <- apply(arHier$Ar,2,median, na.rm = T)
    arMeans <- colMeans(arHier$Ar, na.rm = T)
    arMads <- apply(arHier$Ar,2,mad, constant=1, na.rm = T)
    
    # % polymorphism
    pol <- numeric()
    for (i in popNames(giObj)){
      j <- which(popNames(giObj)==i)
      x <- giObj[pop = i, drop = T]
      p <- sum(isPoly(x, "locus"))/nLoc(x)
      pol[j] <- p
    }
    names(pol) <- popNames(giObj)
  
    return(list("FullBasics"=basics,
                "basics.meds" = basics.pop.meds, 
                "basics.mads" = basics.pop.mad,
                "basics.means" = basics.pop.mean, 
                "basics.sd" = basics.pop.sd, 
                "ar" = arHier$Ar,
                "arMeds" = arMeds, 
                "arMads" = arMads,
                "arMeans" = arMeans,
                "poly" = pol))
    
  }
  
  # get top 5 closest relatives
  getClosestKin <- function(coancSet, numKin=5) {
    tmp <-as.data.frame(t(apply(coancSet, 1, function (x) {
      y <- sort(x, decreasing = T)[1:numKin]
      z <- names(y)
      p <- att[z, "popAbs"]
      return(c(z,p,y))
    })))
    # set column names
    l <- c("ID","Pop","Score")
    names(tmp) <- paste(rep(l, each = numKin), 1:numKin, sep = "")
    # add original population column
    tmp$"Source" <- att[rownames(tmp), "popAbs"]
    tmp$"Match" <- tmp$Source == tmp$Pop1
    return(tmp)
  }
 
  # calculate number and frequency of alleles captured in gardens
  calcAlleleSetFreqs <- function(ghaps,wildAlleles=wAlleles,num_wildAlleles=num_wAlleles){
    
    # Alleles counts and list OG
    gAlleleCounts <- colSums(ghaps@tab, na.rm = T)
    gAlleles <- which(gAlleleCounts>0)
    # gAlleleNamesPos <- gAlleleNames[gAllelesPos]
    num_gAlleles <- length(gAlleles)
    
    # Wild alleles captured in garden
    wINg <- wildAlleles %in% gAlleles
    num_wINg <- sum(wINg)
    num_wNOTg <- sum(!wINg)
    
    freq_wNOTg <- num_wNOTg/num_wildAlleles
    freq_wINg <- round(num_wINg/num_wildAlleles,4)
    
    # Alleles in gardens, not in wild pops
    gINw <- gAlleles %in% wildAlleles
    num_gINw <- sum(gINw)
    num_gNOTw <- sum(!gINw)
    
    freq_gNOTw <- round(num_gNOTw/num_gAlleles,4)
    
    # add to results matrix
    # global_res[num,] <- c(freq_wINg, freq_gNOTw)
    return(c(freq_wINg, freq_gNOTw))
    
  } 
}

# Get data and parse ----------------------------------------------------------

dir <- getwd()

## Read in filtered haps data  ---------------------------------------------------------------
h10 <- read.genepop("/Volumes/CatDisk/SDZooRAD/2024/Rworking/data/allSamplesNeutHaps.gen")
# loc80ind20 <- read.genepop(paste(dir,"data/loc80ind20haps.gen", sep="/"))
loc80ind50 <- read.genepop(paste(dir,"data/loc80ind50haps.gen", sep="/"))
# loc80ind70 <- read.genepop(paste(dir,"data/loc80ind70haps.gen", sep="/"))


## add attributes to genind object  ---------------------------------------------------------------

attFile <- ("elatAttributes.txt")
att <- read.table(attFile, header = T, sep = "\t", colClasses = "character", quote = "")
row.names(att) <- att$ID

# garden plants dropped during filtering
h10 <- addAtts(attFile, h10)
setPop(h10) <- ~Status
table(h10$pop)

droppedGard <- setdiff(indNames(h10[pop="Garden"]),indNames(hapsGard))

att[droppedGard,]

# ID Status Local_Gen Taxon Local_sp popAbs  comLab   RubyID frsGroup
# Elat1   Elat1 Garden  CPGxPoll  Elat  Clifton   seed  seed-1 LF_CL_01  GardLat
# Elat2   Elat2 Garden  CPGxPoll  Elat  Clifton   seed  seed-2 LF_CL_02  GardLat
# Elat4   Elat4 Garden  CPGxPoll  Elat  Clifton   seed  seed-4 LF_CL_05  GardLat
# Elat5   Elat5 Garden  CPGxPoll  Elat  Clifton   seed  seed-5 LF_CL_06  GardLat
# Elat6   Elat6 Garden  CPGxPoll  Elat  Clifton   seed  seed-6 LF_CL_07  GardLat
# Elat22 Elat22 Garden   Trappes  Elat     KNBG   KNBG KNBG-22 LF_KN_20  GardLat

loc80ind20 <- addAtts(attFile, loc80ind20)
loc80ind50 <- addAtts(attFile, loc80ind50)
loc80ind70 <- addAtts(attFile, loc80ind70)
  
## set data set---------------------------------------------------------------
# choose 1
# haps <- loc80ind20
haps <- loc80ind50

## split by G/W into genind objects  ---------------------------------------------------------------

# set pop to gard/wild/planted
setPop(haps) <- ~Status
popNames(haps)

table(haps@pop)

# split into wild and garden
hapsWild <- haps[pop = "Wild", drop = T]
setPop(hapsWild) <- ~popAbs

hapsWildPL <- haps[pop = c("Wild", "planted"), drop = T]
setPop(hapsWildPL) <- ~popAbs
table(hapsWildPL@pop)

hapsGard <- haps[pop = "Garden", drop=T]
table(hapsGard@pop)
setPop(hapsGard) <- ~popAbs

hapsGardPL <- haps[pop = c("Garden", "planted"), drop = T]
setPop(hapsGardPL) <- ~popAbs

# keep all alleles for dapc
hapsWild4dapc <- haps[pop = "Wild"]
setPop(hapsWild4dapc) <- ~popAbs 

hapsGardPL4dapc <- haps[pop = c("Garden", "planted")]
setPop(hapsGardPL4dapc) <- ~popAbs 


# check
popNames(hapsWild)
popNames(hapsWildPL)
popNames(hapsGard)
popNames(hapsGardPL)


# create popmaps for populations stats-----

# change haps object above for each of these:
popMapWildPLInd20 <- cbind(indNames(hapsWildPL), as.character(hapsWildPL@strata$popAbs))
write.table(popMapWildPLInd20, file = "ind20PopMap.txt", sep = "\t", quote = F, row.names = F, col.names = F)

popMapWildPLInd50 <- cbind(indNames(hapsWildPL), as.character(hapsWildPL@strata$popAbs))
write.table(popMapWildPLInd50, file = "ind50PopMap.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# Use populations program to calculate Pi.

# PCA --------------------------------------------------------------------
# Exploratory only
# set data partition

# dat <- hapsWildPL
dat <- hapsWild

# run pca
pcaRes <- prcomp(tab(dat, freq = T, NA.method = "mean"))
summary(pcaRes)
# plot results
plotList <- plotPC(dat, pcaRes)

names(plotList)

plotList$loc12
plotList$loc34
plotList$source12
plotList$source34
plotList$tax12
plotList$locS12
plotList$locS34
plotList$comp12
plotList$comp34


# DAPC --------------------------------------------------------------------

## set data partition -----------------------------------------------------
dat <- hapsWild4dapc
popNames(dat)

dat <- hapsWildPL
popNames(dat)


# remove Strowan - for checking the influence of these plants
# dat <- hapsWildPL[pop=c("VK", "WM", "TR", "GH", "CH", "FS", "B" ), drop = T]
# popNames(dat)
# 
# dat <- hapsWildPL[pop=c("WM", "TR", "GH", "CH", "FS", "ST", "B" ), drop = T]
# popNames(dat)

## a.score optimization ------------------------------
apcVec <- numeric()
da <- dapc(dat, n.da = 3, n.pca = 80)
for (i in 1:10) {
  temp <- optim.a.score(da, n.sim=100, plot = F)
  # names(temp)
  apcVec[i] <- temp$best
}
mean(apcVec)
# hapsPL ind20 = 14
# hapsPL ind50 = 5
# hapsPL ind70 = 10

apc <- round(mean(apcVec),0)
sd(apcVec)
apc <- 5
## by "population" --------------------------------------------------------------------
# = whatever is set in the genind obj

haps.da.pop <- dapc(dat, n.pca=apc, var.loadings=T)#, n.da = 3
# choose 2 LDs


## assign garden plants --------------------------------------------------------------------
# need wild clusters in dapc
# both partitions must have the same loci and alleles
assign.gard <- predict.dapc(haps.da.pop, newdata = hapsGardPL4dapc)

assign.gard$assign
assign.gard$posterior

# Plot DAPC ---------------------------------------------------------------

# partition
set <- dat

# dapc
da <- haps.da.pop

# df of coords and labels
co <- makePlotData(set, da)

# group centroids
cents <- as.data.frame(ddply(co, ~pop, summarize, mean1=mean(LD1), mean2=mean(LD2)))
cents <- as.data.frame(ddply(co, ~pop, summarize, mean1=mean(LD1), mean2=mean(LD2), mean3=mean(LD3)))

ld12 <- ggplot(co) + 
  geom_point(aes(LD1, LD2, color=pop, cex = 3)) + 
  # geom_text_repel(data=cents, x=cents$mean1, y=cents$mean2, label = cents$pop, max.overlaps = 30) +
  geom_text(data=cents, x=cents$mean1, y=cents$mean2, label = cents$pop) +
  theme(legend.position = "left")

# ld23 <- ggplot(co) + 
#   geom_point(aes(LD2, LD3, color=pop, cex = 3)) + 
#   # geom_text_repel(data=cents, x=cents$mean1, y=cents$mean2, label = cents$pop, max.overlaps = 30) +
#   geom_text(data=cents, x=cents$mean2, y=cents$mean3, label = cents$pop) +
#   theme(legend.position = "left")

# Name for figure file
pdfFile = "dapcL80Ind20hapsWildPL.pdf"
pdf(pdfFile, 9,6)
ld12
dev.off()

# frs assignments  --------------------------------------------------------------------

## get coancestry   --------------------------------------------------------------------

frsFileAll = "/Volumes/CatDisk/SDZooRAD/2024/frs/revision/allSamples/run3/hapsInd50_frs_reordered_chunks.out"
coancAll <- read.table(frsFileAll, header = T, row.names = "Recipient", skip = 1)

frsFileWild = "/Volumes/CatDisk/SDZooRAD/2024/frs/revision/wildSamples/run1/hapsWild50Ind_frs_reordered_chunks.out"
coancWild <- read.table(frsFileWild, header = T, row.names = "Recipient", skip = 1)

## subset coanc rows to garden and columns to wild only  --------------------------------------------------------------------
# assign status
wildInd <- which(haps@strata$Status == "Wild")
gardInd <- which(haps@strata$Status != "Wild")

wildNames <- indNames(haps)[wildInd]
gardNames <- indNames(haps)[gardInd]

coancGxW <- coancAll[gardNames,wildNames]

## get closest relatives and populations-------------

# wild populations
fiveClosestWild <- getClosestKin(coancWild)

# garden plants
fiveClosestGard <- getClosestKin(coancGxW)

## compare with dapc-------------

cbind(rownames(fiveClosestGard), rownames(assign.gard$posterior))

## are the two assignments equal?-------------
check <- fiveClosestGard$Pop1 == assign.gard$assign
# rate of equal assignment
assRate <- sum(check)/length(check)
# 0.75
# 0.698 revision

## make one data frame-------------
gardAssBoth <- data.frame(fiveClosestGard, 
                          "dapc" = assign.gard$assign,
                          assign.gard$posterior,
                          "check" = fiveClosestGard$Pop1 == assign.gard$assign)

## save-------------
# write.csv(fiveClosestGard, "gardenFRSAssignmentsWM.csv", quote = F)
write.csv(gardAssBoth, "gardenAssignmentsWM.csv", quote = F)
write.csv(fiveClosestWild, "wildAssignment.csv", quote = F)

# Within Pop Diversity ---------------------------------------------------------------

## set data partition for private alleles -----
# use locality as population not frs clusters b/c these are the units we care about
giObj <- hapsWild
# giObj <- hapsWildPL

# only private alleles within known E. latifrons populations.
# only GH and TR
# giObj <- giObj[pop = c("GH", "TR"), drop = T]

# with WM
# giObj <- giObj[pop = c("GH", "TR", "WM"), drop = T]
# 
# gpObj <- genind2genpop(giObj)
# 
# hierObj <- genind2hierfstat(giObj)

## private alleles----

# output files
paOutfile <- "hapsWildPLLat.pa.csv"
paOutfile <- "hapsWildLat.pa.csv"

# Private alleles with poppr
pa <- writePrivAlleles(giObj,gpObj,hierObj,outfile = paOutfile)

## adegenet summary-----
adSum <- summary(giObj)

## set data partition for diversity stats -----
# use locality as population not frs clusters b/c these are the units we care about
giObj <- hapsWild
giObj <- hapsWildPL

# remove B for within div stats
levels(pop(giObj))
giObj <- giObj[pop = -8, drop = T]

gpObj <- genind2genpop(giObj)

hierObj <- genind2hierfstat(giObj)

## hier.fstat for diversity: F, H, AR-----
popNames(giObj)

inDiv <- getInDiv(giObj, hierObj)

hier.basics <- inDiv$FullBasics

basics.means <- inDiv$basics.means
basics.sd <- inDiv$basics.sd
basics.meds <- inDiv$basics.meds
basics.mads <- inDiv$basics.mads
ar.means <- inDiv$arMeans
ar.meds <- inDiv$arMeds
ar.mads <- inDiv$arMads


write.csv(basics.meds, file = "allWild_Ind50_basics.popmeds.csv")
write.csv(basics.mads, file = "allWild_Ind50_basics.pop_mads.csv")

###bootstraps-----------------------
FisBootHier <- boot(hier.basics$Fis, medianFunc, R=10000, parallel = "multicore", ncpus = 4)
HoBootHier <- boot(hier.basics$Ho, medianFunc, R=10000, parallel = "multicore", ncpus = 4)
HsBootHier <- boot(hier.basics$Hs, medianFunc, R=10000, parallel = "multicore", ncpus = 4)
# ArBootHier <- boot(arHier$Ar, medianFunc, R=10000, parallel = "multicore")

FhierCIlist <- list()
for (i in 1:length(popNames(giObj))) {
  FhierCIlist[[i]] <- boot.ci(FisBootHier, index = i, type = c("basic", "norm"))
}

allFisCI <- do.call(rbind, (lapply(1:length(FhierCIlist), function (x) FhierCIlist[[x]]$normal)))
write.csv(allFisCI, file="allWild_Ind50_FisCI.csv")

HohierCIlist <- list()
for (i in 1:length(popNames(giObj))) {
  HohierCIlist[[i]] <- boot.ci(HoBootHier, index = i, type = c("basic", "norm"))
}

allHoCI <- do.call(rbind,(lapply(1:length(HohierCIlist), function (x) HohierCIlist[[x]]$normal)))
write.csv(allHoCI, file="allWild_Ind50_HoCI.csv")

HshierCIlist <- list()
for (i in 1:length(popNames(giObj))) {
  HshierCIlist[[i]] <- boot.ci(HsBootHier, index = i, type = c("basic", "norm"))
}

allHsCI <- do.call(rbind,(lapply(1:length(HshierCIlist), function (x) HshierCIlist[[x]]$normal)))
write.csv(allHsCI, file="allWild_Ind50_HsCI.csv")


## plot loci distributions----

# Ho
# pdf("Ho_freq.pdf",8,6)
# pdf("Elong_Ho_freq.pdf",8,6)
pdf("Elat_Ho_freq.pdf",8,6)

par(mfrow=c(1,2), oma = c(2,2,0,0))
par(mfrow=c(3,2), oma = c(2,2,0,0))

lapply(1:length(popNames(giObj)),function(x) {
  hist(inDiv$FullBasics$Ho[,x], main = popNames(giObj)[x], cex.main = 2, xlab = "", ylab = "", breaks = (0:20)/20)
  #abline(v=mean(basics$Ho[,x], na.rm = T))
})
mtext("Ho", side = 1, outer = T, cex = 1.3)
mtext("Frequency", side = 2, outer = T, cex = 1.3)

dev.off()

#Hs
# pdf("Elong_Hs_freq.pdf",8,6)
# par(mfrow=c(1,3), oma = c(2,2,0,0))

pdf("Elat_Hs_freq.pdf",8,6)
par(mfrow=c(1,2), oma = c(2,2,0,0))

lapply(1:length(popNames(giObj)),function(x) {
  hist(inDiv$FullBasics$Hs[,x], main = popNames(giObj)[x], cex.main = 2, xlab = "", ylab = "", breaks = (0:20)/20)
})
mtext("Hs", side = 1, outer = T, cex = 1.3)
mtext("Frequency", side = 2, outer = T, cex = 1.3)

dev.off()

#F
# pdf("Elong_F_freq.pdf",8,6)
# par(mfrow=c(1,3), oma = c(2,2,0,0))

pdf("Elat_F_freq.pdf",8,6)
par(mfrow=c(1,2), oma = c(2,2,0,0))

par(mfrow=c(2,3), oma = c(2,2,0,0))

lapply(1:length(popNames(giObj)),function(x) {
  hist(inDiv$FullBasics$Fis[,x], breaks = seq(-1.1,1.1, by=0.2), main = popNames(giObj)[x], cex.main = 2, xlab = "", ylab = "")
  # abline(v=mean(inDiv$FullBasics$Fis[,x], na.rm = T))
})
mtext("F", side = 1, outer = T, cex = 1.3)
mtext("Frequency", side = 2, outer = T, cex = 1.3)

dev.off()


#AR

# pdf("Elong_AR.pdf",8,6)
# par(mfrow=c(1,3), oma = c(2,2,0,0))

pdf("Elat_AR.pdf",8,6)
par(mfrow=c(1,1), oma = c(2,2,0,0))


lapply(1:length(popNames(giObj)),function(x) {
  hist(inDiv$ar[,x], main = popNames(giObj)[x], cex.main = 2, xlab = "", ylab = "")
})
mtext("Allelic Richness", side = 1, outer = T, cex = 1.3)
mtext("Frequency", side = 2, outer = T, cex = 1.3)

dev.off()

# Between pop div ---------------------------------------------------------------



## set data partition for between stats ----------------------------------------
# use locality as population not frs clusters b/c these are the units we care about
giObj <- hapsWildPL
# giObj <- hapsWildfrs
# giObj <- hapsWildLatfrs
# giObj <- hapsWildLongfrs
# giObj <- snpsWild

gpObj <- genind2genpop(giObj)

hierObj <- genind2hierfstat(giObj)

fname <- "elat"
fname <- "elong"
fname <- "all"


## mmod for Gst and D -----------------------------------------------------------
mmodDiff <- diff_stats(giObj)
mmodDiff$global
#        Hs        Ht   Gst_est Gprime_st     D_het    D_mean 
# 0.3056698 0.3962546 0.2286025 0.3643767 0.1491013        NA 

pwDall <- pairwise_D(giObj)
pwGstPall <- pairwise_Gst_Hedrick(giObj)
pwGstNall <- pairwise_Gst_Nei(giObj)

write.csv(mmodDiff$global, paste("mmodDiff_",fname,".csv", sep=""), quote = F)
write.csv(as.matrix(pwDall), paste("pwDall_",fname,".csv", sep=""), quote = F)
write.csv(as.matrix(pwGstPall), paste("pwGstPalll_",fname,".csv", sep=""), quote = F)
write.csv(as.matrix(pwGstNall), paste("pwGstNall_",fname,".csv", sep=""), quote = F)


## poppr for amova---------------------------------------------------------------
# 
# ?dist.dna()
# setPop(loc80ind50) <- ~popAbs
 
amova.res.mean <- poppr.amova(hapsWild, ~popAbs, threads = 16, freq = F, missing = "ignore")
amova.res.mean
# 
# 
# 
# 
## Phi ---------
# nomiss <- missingno(giObj, "mean")
# 
# phiP.res <- Phi_st_Meirmans(nomiss)
# phiP.res$global


# Hoban allele capture---------------------------------------------------------

# species name for table
species_name <- "Elat"

## set data---------------------------------------------------------

# select only E. lat individuals per frs assignment
setPop(haps) <- ~frsGroup
popNames(haps)
Spp_tot_genind <- haps[pop = c("GardLat", "TR", "GH"), drop = T]
popNames(Spp_tot_genind)
# Spp_tot_genind@pop

# set pop slot to gardens and wild pops
setPop(Spp_tot_genind) <- ~popAbs
popNames(Spp_tot_genind)

table(Spp_tot_genind@pop)

# remove single seedling (missing data) and single Strowan individual
Spp_tot_genind <- Spp_tot_genind[pop = -c(1,17), drop = T]

# remove loci not present in all populations/gardens
popList <- seppop(Spp_tot_genind, drop = T)
allLocs <- lapply(popList, locNames)
commonLocs <- Reduce(intersect, allLocs)
length(commonLocs)

# alternative method
# popList <- seppop(Spp_tot_genind, drop = F)
# allele_matrix <- sapply(popList, function(pop) nAll(pop, onlyObserved = T))
# lociKeep <- apply(allele_matrix,1, function(x) all(x>0))
# sum(lociKeep)

Spp_tot_genind <- Spp_tot_genind[loc = commonLocs, drop=T]

wildPops <- c("WM","GH","TR")
gardPops <- popNames(Spp_tot_genind)[which(!popNames(Spp_tot_genind) %in% wildPops)]
wildPopsNums <- which(popNames(Spp_tot_genind) %in% wildPops)

## Part 1: Wild alleles captured in gardens ---------------------------------------------------------------
{ # This block will run the first part, tabulating alleles captured in gardens
  
  # set wd
  setwd("/Volumes/CatDisk/SDZooRAD/2024/Rworking/Revision/Hoban")
  wd <- getwd()
  
  # combine wild pops
  popNames(Spp_tot_genind)[wildPopsNums] <- "Wild"
  popNames(Spp_tot_genind)

  # keep this for the category function
  n_to_drop = 0
  
  # designate pop numbers
  wild_p <- which(popNames(Spp_tot_genind) == "Wild") #2
  garden_p <- which(popNames(Spp_tot_genind) != "Wild") #c(1:13,15) # this includes the 'All' column added below

  n_ind_W<-table(Spp_tot_genind@pop)[wild_p]
  # n_ind_G<-table(Spp_tot_genind@pop)[garden_p];
  Spp_tot_genpop<-genind2genpop(Spp_tot_genind)
  Spp_tot_genind_sep<-seppop(Spp_tot_genind)
  # Spp_tot_genind_sep$"All" <- Spp_tot_genind[pop = 1:13]
  Spp_tot_genind_sep$"All" <- Spp_tot_genind[pop = garden_p]
  
  # add all to garden_p
  garden_p <- which(names(Spp_tot_genind_sep) != "Wild")
  
  # sum the number of copies of each allele in gardens
  # alleles_cap<-colSums(Spp_tot_genind_sep[[garden_p]]@tab,na.rm=T)
  alleles_cap_list <- lapply(Spp_tot_genind_sep[garden_p], function(x) {colSums(x@tab, na.rm = T)})
  
  #Allele categories based only on wild populations (can look at all wild pop'ns or only one if you want)
  # indices of alleles present in each category - a list
  allele_cat_tot<-get.allele.cat(Spp_tot_genpop[wild_p], n_ind_W, n_drop=n_to_drop)
  
  #This goes through each allele category and divides the number captured ex situ (alleles_cap) by the number of alleles existing (allele_cat_tot)
  list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare")
  names(allele_cat_tot) <- list_allele_cat	
  
  alleles_existing_by_sp<-numeric()
  alleles_freq_by_sp<-numeric()
  res <- numeric()  
  
  for (i in 1:length(allele_cat_tot)) {
    # alleles_existing_by_sp[sp,i]<- sum(allele_cat_tot[[i]]>0,na.rm=T)
    alleles_existing_by_sp[i] <- length(allele_cat_tot[[i]])
  }

  # sum the alleles present in gards (>0) / total num wild alleles in cat    
  for (i in 1:length(allele_cat_tot)) { 
    alleles_freq_by_sp[i]<- round(sum(allele_cat_tot[[i]]>0,na.rm=T)/length(allele_cat_tot[[1]]),4)
  }

  # for each garden, for each category get proportion of wild alleles captured
  wild_results_list <- lapply(alleles_cap_list, function (cap) {
      for (l in 1:length(allele_cat_tot)) {
        res[l] <- round(sum(cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
      }
      return(res)
  })
  
  wild_results <- matrix(unlist(wild_results_list), ncol = 5, byrow = T)

  # Compile results
  finalRes <- rbind(alleles_existing_by_sp, alleles_freq_by_sp, wild_results)
  # numPl <- unlist(lapply(Spp_tot_genind_sep, nInd))
  # moveAll <- c(length(numPl),1:(length(numPl)-1))
  # numPl <- numPl[moveAll]
  # finalRes <- cbind(finalRes,numPl)
  rownames(finalRes) <- c(paste(species_name,"wild_num", sep = "_"), paste(species_name,"wild_freq", sep = "_"), names(wild_results_list)) # names(Spp_tot_genind_sep)[garden_p]
  colnames(finalRes) <- c(list_allele_cat) #, "Num Plants")
  finalRes
  write.csv(finalRes,file=paste(wd,"/wildAlleles_in_gardens.csv",sep=""))
  
  # mean proportions across gardens
  colMeans(finalRes[names(wild_results_list),1:5])
  #     global glob_v_com   glob_com glob_lowfr  glob_rare 
  # 0.6229833  0.8523167  0.3456917  0.1999333        NaN 
  
}

## Part 2: compare allele presence between wild pops ---------------------------------------------------------------
{
# set wd
setwd("/Volumes/CatDisk/SDZooRAD/2024/Rworking/Revision/Hoban")
wd <- getwd()

# extract only wild pops
# reset pop names
popNames(Spp_tot_genind)
setPop(Spp_tot_genind) <- ~popAbs

#set wild pops
wildHaps <- Spp_tot_genind[pop = wildPopsNums, drop = T]
popNames(wildHaps)

# list to store allele vectors for each pop
wResList <- vector(mode = "list", length = length(popNames(wildHaps)))
names(wResList) <- popNames(wildHaps)

# count and list alleles in each pop
num_Alleles <- integer()
wRestlist <- list()
for (w in 1:length(popNames(wildHaps))) {
  wp <- wildHaps[pop=w]
  alleleCounts <- colSums(wp@tab, na.rm = T)
  Alleles <- which(alleleCounts>0)
  num_Alleles[w] <- length(Alleles)
  wResList[[w]] <- Alleles
  
}

# find alleles present in one but not another pop

# get all pairwise comparisons, both ways
allComb <- expand.grid(names(wResList),names(wResList))
combList <- apply(allComb,1,as.vector, simplify = F)
lk <- unlist(lapply(combList, function(x) {x[1] != x[2]}))
combList <- combList[lk]
combList

# first in pairs for later
combFirst <- unlist(lapply(combList, function(x) {x[1]}))

# compare each pair
diffList <- lapply(combList, function(x) {setdiff(wResList[[ x[[1]] ]], wResList[[ x[[2]] ]] )})
dnames <- lapply(combList, function(x) paste(x[[1]],"NotIn", x[[2]], sep = "_"))
names(diffList) <- dnames

# number of alleles in 1 and not in 2
numNotShared <- unlist(lapply(diffList, length))

# num of alleles in each pop
numAll <- unlist(lapply(wResList[combFirst], length))

# freq of missing alleles
freqNotShared <- round(numNotShared/numAll,3)

barplot(numNotShared)
barplot(freqNotShared)

wildRes <- rbind(numNotShared,freqNotShared)

write.csv(wildRes, paste(wd,"/wildAllelesNotShared.csv", sep = ""))

}

## Part3 Reverse: find alleles present in Gardens not present in Wild---------------------------------------------------------------
# still using genind from above

{ ## This block will run the third part
  
# set wd
setwd("/Volumes/CatDisk/SDZooRAD/2024/Rworking/Revision/Hoban")
wd <- getwd()

# make sure pop factor is correct
setPop(Spp_tot_genind) <- ~popAbs

# get wild allele list and counts
wAlleleCounts <- colSums(Spp_tot_genind[pop=wildPops]@tab, na.rm = T)
wAlleles <- which(wAlleleCounts>0)
num_wAlleles <- length(wAlleles)

# using seppop to set up loop
Spp_tot_genind_sep<-seppop(Spp_tot_genind[pop = gardPops])
Spp_tot_genind_sep$"All" <- Spp_tot_genind[pop = gardPops]

# calculate freq of gard and wild alleles
# glob_resList <- sapply(Spp_tot_genind_sep, calcAlleleSetFreqs)
global_res <- t(sapply(Spp_tot_genind_sep, calcAlleleSetFreqs))
colnames(global_res) <- c("WildCaptured", "GardenOnly")

# View results
global_res

# save to file
write.csv(global_res, paste(wd, "/gardenAllels_not_in_wild.csv", sep=""))

}


# ---------------------------------------------------------------
# recalculate D. holm G'st

dholm <- read.genepop("/Volumes/SysMaticsLab/2013-2023/HD4/Dioon_RADseq_analysis/Dholm/Rworking/sixPops/haps_7070.gen")
popNames(dholm) <- c("Ixtayutla", "Jamiltepec", "Juchatengo", "Loxicha", "Rancho Limon", "Textitlan")

## mmod for Gst and D -----------------------------------------------------------
dHolm_mmodDiff <- diff_stats(dholm)
dHolm_mmodDiff$global

dHolm_pwDall <- pairwise_D(dholm)
dHolm_pwGstPall <- pairwise_Gst_Hedrick(dholm)
dHolm_pwGstNall <- pairwise_Gst_Nei(dholm)




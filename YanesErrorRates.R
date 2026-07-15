################################################################################
# Nov. 2023
# Optimizing Stacks de novo assembly parameters.
# Script to calculate several error rates from varying values for Stacks parameters m, M, n.
# Compiles these error rates and depth of coverage and plots results
# Requires output from m_test.sh, Mtest.sh, and ntest.sh scripts.
################################################################################

library(stringr)
library(ggplot2)
library(plyr)
library(graphics)


setwd('/Volumes/CatDisk/SDZooRAD/2024/stacks/tests')
homedir <- getwd()

######################################################
# set directory - comment out to specify individual or combined analyses
######################################################

# prefix_m <- paste(homedir,"/little_m/m3-10M2n2hap3/m", sep ="")
prefix_m <- paste(homedir,"/little_m/m3-10M2n2hap3/r80/m", sep ="")

prefix_n <- paste(homedir,"/n/m7M3n4-7hap3/r80/n", sep ="")
# 
prefix_M <- paste(homedir,"/big_m/m7Mn1-5mac3/r80/m", sep ="")
# 
# prefix_opt_nPlus <- paste(homedir, "", sep="")
# 
# prefix_opt <- paste(homedir, "", sep="")
# 
# prefix_OGopt <- paste(homedir, "", sep = "")

######################################################

for ( prefix in c(prefix_m, prefix_M, prefix_n) ) {

    # set of parameter values ### Could get these from the folder names
  if (prefix == prefix_m) {
    ms <- c(3:10)
  }
  if (prefix == prefix_M) {
    ms <- c(1:5)
  }
  if (prefix == prefix_n) {
    ms <- c(4:7)
  }
  # if (prefix == prefix_opt | prefix == prefix_opt_nPlus) {
  #   ms <- c(1:4)
  # }
  # if (prefix == prefix_OGopt) {
  #     ms <- c(1:6)
  # }

# prefix <- "/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/stacks/fullSet/"
# ms <- "opt_m6M2n3"
  samples <- c("Elat12", "Elat34", "Elat55", "Elat73", "Elat95", "Elat111", "Elat154") # "Elat12rep", "Elat34rep", "Elat55rep", "Elat73rep", "Elat95rep", "Elat111rep", "Elat154rep")
  reps <- c("Elat12", "Elat34", "Elat55", "Elat73", "Elat95", "Elat111", "Elat154", "Elat12rep", "Elat34rep", "Elat55rep", "Elat73rep", "Elat95rep", "Elat111rep", "Elat154rep")
  repNums <- c(12, 34, 55, 73, 95, 111, 154)

  # list to hold dfs
  # df.list <- list("47" = data.frame(), "48" = data.frame())
  df.list <- lapply(samples, function(x) {
    x <- data.frame(row.names = ms,
                    parVal = rep(0, length(ms)),
                    num.loci = rep(0, length(ms)),
                    num.missing.loci = rep(0, length(ms)),
                    num.common = rep(0, length(ms)),
                    prop.loci.missing = rep(0, length(ms)),
                    missing.one = rep(0, length(ms)),
                    p.missing.one = rep(0, length(ms)),
                    l.error.rate = rep(0, length(ms)),
                    num.mismatches.alleles = rep(0, length(ms)),
                    a.error.rate = rep(0, length(ms)),
                    num.snps = rep(0, length(ms)),
                    num.snp.matches = rep(0, length(ms)),
                    snp.error.rate = rep(0, length(ms)),
                    meanCov = rep(0, length(ms)),
                    stdevCov = rep(0, length(ms)),
                    minCov = rep(0, length(ms)),
                    maxCov = rep(0, length(ms)),
                    r80Cov = rep(0, length(ms))
    )})
  
  names(df.list) <- samples
  
  ######################################################
  # loop through the parameter values and record info on loci and errors
  
  # index variable
  i=1
  
  for (m in ms){
  
    # HAPLOTYPES
    haps.vcf <- paste(prefix,m,"/populations.haps.vcf.gz", sep="")
    haps <- read.table(file=haps.vcf, sep="\t", header=T, row.names = 1, stringsAsFactors = F, comment.char = "", skip = 15)
  
    # get replicate set only
    haps <- haps[,reps]
  
    # SNPS
    snps.vcf <- paste(prefix,m, "/populations.snps.vcf.gz", sep="")
    snps <- read.table(file=snps.vcf, sep="\t", header=T, stringsAsFactors = F, comment.char = "", skip = 15)
  
    # get replicate set only
    snps <- snps[,reps]

    # Get coverage from log file     
    # covFile <- paste(prefix,m,"/defaultPop/gstacks.log", sep = "")
    covFile <- paste(prefix,m,"/gstacks.log", sep = "")
    
    covFile <- gsub('r80/','', covFile) # for r80 loci
    
    cmd <- paste("grep 'effective per-sample' ", covFile, sep = "")
    cov <- system(cmd, intern=T)
    cov.split <- strsplit(cov, " ")
    cov.vals <- as.numeric(gsub("[a-z,=]+", "", x=cov.split[[1]][6:9], perl=T))
    
    # Get coverage from vcf file
    snpcovfile <- paste(prefix,m,"/snpcov.txt", sep = "")
    snpCovAll <- read.table(snpcovfile, T)
    snpCovMean <- mean(snpCovAll$MEAN_DEPTH)
    
    # 
    # snpDepth <- system2("vcftools", ar, stdout = T)
    # snpDepth <- as.numeric(snpDepth)
    # meanCov <- mean(snpDepth)
    # sdCov <- sd(snpDepth)
    # minCov <- min(snpDepth)
    # maxCov <- max(snpDepth)
    
    
    # Calculate stats	for each rep set
    # for (k in repNums){
    for (k in samples){
  
      # name the samples
      rep1 <- k
      rep2 <- paste(k,"rep", sep = "")
      
      # Remove Info from genotypes
      snps[rep1] <- apply(snps[rep1], 1, function(x) str_split_1(x,":")[1])
      snps[rep2] <- apply(snps[rep2], 1, function(x) str_split_1(x,":")[1])
      
  
      #############################################################################################
      # Calculate error rates of Mastretta-Yanes et al. (2014).
      
      ####### haps ######
      
      # total number loci
      num.loci <- length(haps[[rep1]])
      
      # Number missing loci = count of loci missing in >= 1 replicate
      num.missing.loci <- length(which(haps[[rep1]]=="./." | haps[[rep2]]=="./."))
      
      # Proportion of missing loci = number of missing loci/total number of loci across samples
      prop.loci.missing <- num.missing.loci/num.loci
      
      # Proportion missing in both or only one replicate
      missing.both <- length(which(haps[[rep1]]=="./." & haps[[rep2]]=="./."))
      missing.one <- num.missing.loci-missing.both
      p.missing.diff <- missing.one/num.missing.loci
      
      # Locus error rate = number of shared loci with one locus missing in a replicate / total number of loci across samples
      l.error.rate <- missing.one/num.loci
      
      # loci common to both samples
      common.loci <- haps[which(haps[[rep1]] != "./." & haps[[rep2]] != "./."),] # vcf
      
      # number of common loci
      num.common <- length(common.loci[,1])  #haps$a!="-" & haps$b!="-"))
      
      # identical common loci
      num.matches <- length(which(common.loci[[rep1]]==common.loci[[rep2]]))
      # littlem_df.list[[k]]$num.matches <- num.matches
      
      # mismatches <- which(common.loci[[rep1]]!=common.loci[[rep2]])
      
      num.mismatches <- num.common-num.matches
      
      a.error.rate <- num.mismatches/num.common
  
      ####### snps ######
      
      # snps
      num.snps <- length(snps[[rep1]])
      
      # loci common to both samples
      common.snps <- snps[which(snps[[rep1]] != "./." & snps[[rep2]] != "./."),]
      
      num.common.snps <- length(common.snps[[rep1]])
      
      num.snp.matches <- length(which(common.snps[[rep1]]==common.snps[[rep2]]))
      
      num.snp.mismatches <- num.common.snps-num.snp.matches
      
      snp.error.rate <- num.snp.mismatches/num.common.snps
      
      # Data frame
      df.list[[k]]$parVal[i] <- m
      df.list[[k]]$num.loci[i] <- num.loci
      df.list[[k]]$num.missing.loci[i] <- num.missing.loci 
      df.list[[k]]$num.common[i] <- num.common
      df.list[[k]]$prop.loci.missing[i] <- prop.loci.missing
      df.list[[k]]$missing.one[i] <- missing.one
      df.list[[k]]$p.missing.one[i] <- p.missing.diff
      df.list[[k]]$l.error.rate[i] <- l.error.rate
      df.list[[k]]$num.mismatches.alleles[i] <- num.mismatches
      df.list[[k]]$a.error.rate[i] <- a.error.rate
      df.list[[k]]$num.snps[i] <- num.snps
      df.list[[k]]$num.snp.matches[i] <- num.snp.matches
      df.list[[k]]$snp.error.rate[i] <- snp.error.rate
      df.list[[k]]$meanCov[i] <- cov.vals[1]
      df.list[[k]]$stdevCov[i] <- cov.vals[2]
      df.list[[k]]$minCov[i] <- cov.vals[3]
      df.list[[k]]$maxCov[i] <- cov.vals[4]
      df.list[[k]]$r80Cov[i] <- snpCovMean
      # df.list[[k]]$meanCov[i] <- meanCov
      # df.list[[k]]$stdevCov[i] <- sdCov
      # df.list[[k]]$minCov[i] <- minCov
      # df.list[[k]]$maxCov[i] <- maxCov
  
    }
    
    i <- i+1
    
  }
  
  if (prefix == prefix_m) {
    m.list <- df.list
  }
  if (prefix == prefix_M) {
    M.list <- df.list
  }
  if (prefix == prefix_n) {
    n.list <- df.list
  }
  # if (prefix == prefix_opt) {
  #   opt.list <- df.list
  # }
  # if (prefix == prefix_opt_nPlus) {
  #   opt.list_nPlus <- df.list
  # }
  # if (prefix == prefix_OGopt) {
  #   og.n.list <- df.list
  # }
  

}



m.list
n.list
M.list
opt.list
og.n.list

# get means of each parameter across replicates
                 # convert each df in list to matrix
                 # and return as one array                  # split the array
                                                            # and apply mean to
                                                            # dims 2 and 3      # remove first column 'X'
mean_m <- as.data.frame(aaply(  laply(m.list, as.matrix),         c(2,3), mean))#[,-1]
sd_m <- as.data.frame(aaply(  laply(m.list, as.matrix),         c(2,3), sd)) #[,-1]
mean_M <- as.data.frame(aaply(laply(M.list, as.matrix), c(2,3), mean)) #[,-1]
mean_n <- as.data.frame(aaply(laply(n.list, as.matrix), c(2,3), mean)) #[,-1]
mean_opt <- as.data.frame(aaply(laply(opt.list, as.matrix), c(2,3), mean)) #[,-1]
mean_opt_nPlus <- as.data.frame(aaply(laply(opt.list_nPlus, as.matrix), c(2,3), mean)) #[,-1]
mean_og_n <- as.data.frame(aaply(laply(og.n.list, as.matrix), c(2,3), mean))

mean_m
mean_M



# for single par values
# mean_single <- (aaply(laply(M.list, as.matrix), 2, mean))

write.csv(mean_m, file="mean_m.csv")
write.csv(mean_M, file="meanM.csv")
write.csv(mean_n, file="mean_n.csv")
write.csv(mean_opt_nPlus, file = "mean_opt_nPlus.csv")
write.csv(mean_opt, file = "mean_opt.csv")
write.csv(mean_og_n, file = "mean_opt_ogn.csv")

mean_m <- read.csv(file = "mean_m.csv")
mean_M <- read.csv(file = "meanM.csv")
mean_n <- read.csv(file = "mean_n.csv")
mean_opt_nPlus <- read.csv(file = "mean_opt_nPlus.csv")
mean_opt <- read.csv(file = "mean_opt.csv")

# mean_M <- rbind(mean_single, mean_M[,-1])

#############################################################################################
# Plots
?layout()
dev.off()

plotNloci <- function () {
  plot(mean_m$parVal,mean_m$num.loci,
       xlim = c(0,10),
       main = "Number r80 Loci", xlab = "", ylab = "", type = "n") #ylim=c(5000,20000),yaxt = "n"
  # axis(2, at = seq(5000,20000,2500), labels = seq(5000,20000,2500))
  points(mean_m$parVal,mean_m$num.loci, pch=1, cex=2)
  lines(mean_m$parVal,mean_m$num.loci)
  points(mean_M$parVal,mean_M$num.loci, pch=2, cex=2)
  lines(mean_M$parVal,mean_M$num.loci)
  points(mean_n$parVal,mean_n$num.loci, pch=5, cex=2)
  lines(mean_n$parVal,mean_n$num.loci)
  # points(mean_opt$parVal,mean_opt$num.loci, pch=15, cex=2)
  # lines(mean_opt$parVal,mean_opt$num.loci, lty=4)
  # points(mean_opt_nPlus$parVal,mean_opt_nPlus$num.loci, pch=17, cex=2)
  # lines(mean_opt_nPlus$parVal,mean_opt_nPlus$num.loci, lty=4)
  # points(mean_og_n$parVal,mean_og_n$num.loci, pch=18, cex=2)
  # lines(mean_og_n$parVal,mean_og_n$num.loci, lty=4)
  # lines(x=mean_opt["num.loci"], y=1:10, type = "b", pch=15, lty=4)
  # abline(h=mean_opt["num.loci"], lty = 4)
  # legend(x=0, y=15000, legend = c("m", "M", "n", "Opt", "Opt2"),
  #        pch = c(1,2,5,15,17), 
  #        lty = c(0,0,0,4,4), 
  #        inset = 0.05, 
  #        y.intersp = 0.5)
}
# pt.cex=c(1,1,1,1,1)
plotLerror <- function () {
  plot(mean_m$parVal,mean_m$l.error.rate,xlim = c(0,10), ylim=c(0,0.25), main = "Loci Error Rate", xlab = "", ylab = "", cex=2)
  lines(mean_m$parVal,mean_m$l.error.rate)
  points(mean_M$parVal,mean_M$l.error.rate, pch=2, cex=2)
  lines(mean_M$parVal,mean_M$l.error.rate, pch=2)
  points(mean_n$parVal,mean_n$l.error.rate, pch=5, cex=2)
  lines(mean_n$parVal,mean_n$l.error.rate, pch=5)
  # points(mean_opt$parVal,mean_opt$l.error.rate, pch=15, cex=2)
  # lines(mean_opt$parVal,mean_opt$l.error.rate, lty=4)
  # points(mean_opt_nPlus$parVal,mean_opt_nPlus$l.error.rate, pch=17, cex=2)
  # lines(mean_opt_nPlus$parVal,mean_opt_nPlus$l.error.rate, lty=4)
  # points(mean_og_n$parVal,mean_og_n$l.error.rate, pch=18, cex=2)
  # lines(mean_og_n$parVal,mean_og_n$l.error.rate, lty=4)
  # lines(mean_opt_nPlus$parVal,mean_opt_nPlus$l.error.rate, lty=4)
  
  # abline(h=mean_opt["l.error.rate"], lty=4)
}


plotAerror <- function() {plot(mean_m$parVal,mean_m$a.error.rate,
                               xlim = c(0,10),
                               ylim=c(0,0.25), main = "Allele Error Rate", xlab = "", ylab = "", cex=2)
  lines(mean_m$parVal,mean_m$a.error.rate)
  points(mean_M$parVal,mean_M$a.error.rate, pch=2, cex=2)
  lines(mean_M$parVal,mean_M$a.error.rate, pch=2)
  points(mean_n$parVal,mean_n$a.error.rate, pch=5, cex=2)
  lines(mean_n$parVal,mean_n$a.error.rate, pch=5)
  # points(mean_opt$parVal,mean_opt$a.error.rate, pch=15, cex=2)
  # points(mean_opt_nPlus$parVal,mean_opt_nPlus$a.error.rate, pch=17, cex=2)
  # lines(mean_opt_nPlus$parVal,mean_opt_nPlus$a.error.rate, lty=4)
  # points(mean_og_n$parVal,mean_og_n$a.error.rate, pch=18, cex=2)
  # lines(mean_og_n$parVal,mean_og_n$a.error.rate, lty=4)
  # abline(h=mean_opt["a.error.rate"], lty=4)
}

plotSerror <- function () {plot(mean_m$parVal,mean_m$snp.error.rate,
                                xlim = c(0,10),
                                ylim=c(0,0.05), main = "SNP Error Rate", xlab = "", ylab = "", cex=2)
  lines(mean_m$parVal,mean_m$snp.error.rate)
  points(mean_M$parVal,mean_M$snp.error.rate, pch=2, cex=2)
  lines(mean_M$parVal,mean_M$snp.error.rate, pch=2)
  points(mean_n$parVal,mean_n$snp.error.rate, pch=5, cex=2)
  lines(mean_n$parVal,mean_n$snp.error.rate, pch=5)
  # points(mean_opt$parVal,mean_opt$snp.error.rate, pch=15, cex=2)
  # lines(mean_opt$parVal,mean_opt$snp.error.rate, lty=4)
  # points(mean_opt_nPlus$parVal,mean_opt_nPlus$snp.error.rate, pch=17, cex=2)
  # lines(mean_opt_nPlus$parVal,mean_opt_nPlus$snp.error.rate, lty=4)
  # points(mean_og_n$parVal,mean_og_n$snp.error.rate, pch=18, cex=2)
  # lines(mean_og_n$parVal,mean_og_n$snp.error.rate, lty=4)
  # abline(h=mean_opt["snp.error.rate"], lty=4)
}

plotCov <- function () {plot(mean_m$parVal,mean_m$meanCov,
                             xlim = c(0,10),
                             ylim=c(0,80), main = "Mean Coverage", xlab = "", ylab = "", cex=2)
  lines(mean_m$parVal,mean_m$meanCov)
  points(mean_m$parVal,mean_m$r80Cov, pch=3, cex=2)
  lines(mean_m$parVal,mean_m$r80Cov)
  points(mean_M$parVal,mean_M$meanCov, pch=2, cex=2)
  lines(mean_M$parVal,mean_M$meanCov, pch=2)
  points(mean_n$parVal,mean_n$meanCov, pch=5, cex=2)
  lines(mean_n$parVal,mean_n$meanCov, pch=5)
  # points(mean_opt$parVal,mean_opt$meanCov, pch=15, cex=2)
  # lines(mean_opt$parVal,mean_opt$meanCov, lty=4)
  # points(mean_opt_nPlus$parVal,mean_opt_nPlus$meanCov, pch=17, cex=2)
  # lines(mean_opt_nPlus$parVal,mean_opt_nPlus$meanCov, lty=4)
  # points(mean_og_n$parVal,mean_og_n$meanCov, pch=18, cex=2)
  # lines(mean_og_n$parVal,mean_og_n$meanCov, lty=4)
  # abline(h=mean_opt["meanCov"], lty=4)
  # legend("bottomleft", legend = c("m", "M", "n", "Opt"), 
  #        pch = c(1,2,5,1), pt.cex=c(1,1,1,0), 
  #        lty = c(0,0,0,4), inset = 0.05, y.intersp = 1)
}


pdf(file = "AssemblyParameterTests.pdf",8.5,6)
# pdf(file = "AssemblyLittle_m_Tests.pdf",8.5,6)
mat <- matrix(nrow=2, ncol=3, data=c(1,2,3,4,5,6), byrow = T)
layout(mat)
plotCov()
plotNloci()
plotLerror()
plotAerror()
plotSerror()
plot(1, type="n", axes=F, ylab="", xlab="")
legend(x="center", legend = c("m (M1 N1)", "SNP Coverage", "M=N (m7)", "n (m7 M3)"),
       pch = c(1,3,2,5), 
       lty = c(0,0,0,0), 
       y.intersp = rep(1.5,4), #c(0.7,0.7,0.7,0.7),
       cex = 1.75,
       title = "Parameters",
       title.cex = 2,
       bty = "n")

dev.off()

?legend


x <- pmin(3, pmax(-3, stats::rnorm(50)))
y <- pmin(3, pmax(-3, stats::rnorm(50)))
xhist <- hist(x, breaks = seq(-3,3,0.5), plot = FALSE)
yhist <- hist(y, breaks = seq(-3,3,0.5), plot = FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(-3, 3)
yrange <- c(-3, 3)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(3,1), c(1,3), TRUE)
layout.show(nf)
  
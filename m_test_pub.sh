#! /bin/zsh
setopt extendedglob


#######################################################################################
# Script to run stacks denovo_map.pl iterating m parameter values
# Modified from D. tomasellii study 
# Create folders with parameter names and set directories before running
#######################################################################################



sampdir=/Volumes/CatDisk/SDZooRAD/2024/cleaning/bbduk
infodir=/Volumes/CatDisk/SDZooRAD/2024/stacks/info
outdir=/Volumes/CatDisk/SDZooRAD/2024/stacks/tests/little_m

##############################################################
# Full optimization set
##############################################################

for m in {3..10}; do

# set directory and change to it
	testdir=$outdir/m3-10M2n2hap3/m$m
	cd $testdir

# run pipeline	
	denovo_map.pl \
	--samples $sampdir \
	--popmap $infodir/elat_optPop.txt \
	-o $testdir \
	--paired \
	-T 18 \
	-M 2 \
	-n 2 \
	-m $m \
	-X "populations: --vcf" \
	
	echo
	date
	echo
done &> $outdir/mTest.oe.txt



# run r80 populations and calculate snp depth

for m in {3..10}; do

# set directory and change to it
	dir80=$outdir/m3-10M2n2hap3/r80/m$m
	mkdir -p $dir80
	cd $dir80
	
	/usr/local/bin/populations \
	-P $outdir/m3-10M2n2hap3/m$m \
	-O $dir80 \
	-M $infodir/elat_optPop.txt \
	-r 0.8 \
	-t 8 \
	--vcf-all &> $dir80/pop.oe.txt 
	
	gzip $dir80/*vcf
	
	v=$dir80/populations.snps.vcf.gz
	vcftools --gzvcf $v --site-mean-depth --stdout | cut -f 3 > $dir80/sitecov.txt

done

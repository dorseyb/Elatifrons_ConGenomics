#!/opt/homebrew/bin/bash
# setopt extendedglob


#######################################################################################
# Script to run stacks denovo_map.pl iterating M parameter values
# Modified from D. tomasellii study 
# Create folders with parameter names and set directories before running
#######################################################################################

sampdir=/Volumes/CatDisk/SDZooRAD/2024/cleaning/bbduk
infodir=/Volumes/CatDisk/SDZooRAD/2024/stacks/info
outdir=/Volumes/CatDisk/SDZooRAD/2024/stacks/tests/big_m

# get sample names from popmap
mapfile -t names < <(cat $infodir/elat_optPop.txt | cut -f 1)
# 
# # # for each value of M
for m in {1..5}; do
# # 
# # make a directory to hold output
	testdir=$outdir/m7Mn1-5mac3/m$m
	mkdir -p $testdir
	
	for x in 0 10 20 30 40 50 60; do
		
		echo "batch"$x
		echo
	
		for f in "${names[@]:$x:10}"; do
		
			/usr/local/bin/ustacks \
			--in-type gzfastq \
			--file $sampdir/$f.1.fq.gz \
			--out-path $testdir \
			-m 7 \
			--name $f \
			-t 2 \
			-M $m &> $testdir/$f.ustacks_oe.txt & 
	
		done
		wait
		echo
		echo "batch"$x"done"
		echo
		
	done
	
	wait
 	
# # run the rest of the pipeline
	
	denovo_map.pl \
	--resume \
	--samples $sampdir \
	--popmap $infodir/elat_optPop.txt \
	-o $testdir \
	--paired \
	-T 18 \
	-M $m \
	-n $m \
	-m 3 \
	-X "populations: --vcf" \

done

echo
echo "denovo finished, starting r80 populations"
echo

# run the r80 populations and get snp depth
for m in {1..5}; do

# set directory and change to it
	dir80=$outdir/m7Mn1-5mac3/r80/m$m
# 	mkdir -p $dir80
	cd $dir80
# 	
# 	/usr/local/bin/populations \
# 	-P $outdir/m7Mn1-5mac3/m$m \
# 	-O $dir80 \
# 	-M $infodir/elat_optPop.txt \
# 	-r 0.8 \
# 	-t 18 \
# 	--vcf &> $dir80/pop.oe.txt 
	
# 	gzip "${dir80}"/*vcf
	
	v=$dir80/populations.snps.vcf.gz
	vcftools --gzvcf $v --site-mean-depth --stdout | cut -f 3 > $dir80/snpcov.txt

done






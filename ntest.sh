#!/bin/zsh
setopt extendedglob

sampdir=/Volumes/CatDisk/SDZooRAD/2024/cleaning/bbduk
udir=/Volumes/CatDisk/SDZooRAD/2024/stacks/tests/big_m/m7Mn1-5mac3/m3
tdir=/Volumes/CatDisk/SDZooRAD/2024/stacks/tests/n/m7M3n4-7hap3
pmap=/Volumes/CatDisk/SDZooRAD/2024/stacks/info/elat_optPop.txt


# run pipeline without rerunning ustacks for every value of n
for n in {4..7}; do
	
	# make folder for this test
	testdir=$tdir/n$n
	
	if [[ ! -d $testdir ]]
	then
		mkdir $testdir
	fi
	
	# copy ustacks files to new folder
	
	for file in $udir/[!catalog]*alleles.tsv.gz
		do
			if [[ ! -f $testdir/$file:t ]]
			then
				cp $file $testdir
			fi
		done
	
	for file in $udir/[!catalog]*snps.tsv.gz
		do
			if [[ ! -f $testdir/$file:t ]]
			then
				cp $file $testdir
			fi
		done
	
	for file in $udir/[!catalog]*tags.tsv.gz
		do
			if [[ ! -f $testdir/$file:t ]]
			then
				cp $file $testdir
			fi
		done	

	# run pipeline
	denovo_map.pl \
	--samples $sampdir \
	--popmap $pmap \
	-o $testdir \
	--paired \
	-T 18 \
	-M 3 \
	-n $n \
	-m 7 \
	--resume \
	-X "populations: --vcf"
	
done &> $tdir/nTests.oe.txt


# run r80 populations and get snp depth

for n in {4..7}; do

# set directory and change to it
	dir80=$tdir/r80/n$n
	mkdir -p $dir80
	cd $dir80
	
	/usr/local/bin/populations \
	-P $tdir/n$n \
	-O $dir80 \
	-M $pmap \
	-r 0.8 \
	-t 18 \
	--vcf &> $dir80/pop.oe.txt 
	
	gzip $dir80/*vcf
		
	v=$dir80/populations.snps.vcf.gz
	vcftools --gzvcf $v --site-mean-depth --stdout | cut -f 3 > $dir80/snpcov.txt

done


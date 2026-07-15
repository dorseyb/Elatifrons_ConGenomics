use strict;
use warnings;


my $file = '/Volumes/CatDisk/SDZooRAD/2024/stacks/noReps/revision/popStatsWL/ind50/populations.samples.fa';
my $outdir = '/Volumes/CatDisk/SDZooRAD/2024/phylo/Revision/fasta/diploid/';
my $outdir0 = '/Volumes/CatDisk/SDZooRAD/2024/phylo/Revision/fasta/allele0/';
my $outdir1 = '/Volumes/CatDisk/SDZooRAD/2024/phylo/Revision/fasta/allele1/';

open my $fh, '<', $file or die;

my $fastaref = fasta_to_hash($fh);
my %fasthash = %{$fastaref};

close $fh or die;

my %newhash = ();
my %allele0hash = ();
my %allele1hash = ();

# go through each sequence defline
foreach my $def (keys %fasthash) {


# split the defline to get the locus i.e CLocus_x
	my @d = split /_Sample/, $def;
	my $loc = $d[0];
	my $locnum = $loc;
	$locnum =~ s/CLocus_//;
	
# 	print $locnum."\n";

	my @indArray = split /\[/, $def;
	$indArray[1] =~ s/\]//;
# 	$indArray[1] =~ s/[[]//;
	my $indiv = $indArray[1];
	
#  	print $indiv."\n";

# 	next unless ( (grep /\b$locnum\b/, @markers) && (grep /\b$indiv\b/, @ind) );

# Clean up the defline
	my @defsplit = split / /, $def;
	$defsplit[1] =~ s/[\]]//;
	$defsplit[1] =~ s/[[]//;
	my @splitAllele = split /_/, $defsplit[0];
	my $allele = join("_", @splitAllele[6,7]);
	my $newdef = join("_", $defsplit[1], $allele);
	
# add locus (as ref to anonymous hash) to the hashes for holding hashes of all seqs of each locus

# diploid
	$newhash{$loc} = {} unless exists $newhash{$loc};

# allele 0	
	$allele0hash{$loc} = {} unless exists $allele0hash{$loc};
# allele 1	
	$allele1hash{$loc} = {} unless exists $allele1hash{$loc};

# add defline (key) and seq (value) to new hashes

# diploid
	$newhash{$loc}{$newdef} = $fasthash{$def};

# allele 0	
	if ($splitAllele[7] == 0) {
	print $loc."Allele 0"."\n";
		$allele0hash{$loc}{$newdef} = $fasthash{$def};
	}

# allele 1	
	if ($splitAllele[7] == 1) {
		print $loc."Allele 1"."\n";
		$allele1hash{$loc}{$newdef} = $fasthash{$def};
	}
	
	
}

# go through new locus hash 

# print out diploid seqs
 foreach my $locus (keys %newhash){

# dereference the hash for a locus
	my %lochash = %{ $newhash{$locus} };
# open fasta file	
	open my $fh2, ">>", $outdir.$locus.'.fa';
	
	foreach my $hap (sort(keys %lochash)) {
		print $fh2 ">".$hap."\n";
		print $fh2 $lochash{$hap}."\n";
	}
	close $fh2 or die;
	
}

# print out allele 0		
foreach my $all0 (keys %allele0hash){

	my %a0hash = %{ $allele0hash{$all0} };
	
	open my $fh3, ">>", $outdir0.$all0.'_0.fa';
	
	foreach my $defline (sort(keys %a0hash)) {
		print $fh3 ">".$defline."\n";
		print $fh3 $a0hash{$defline}."\n";
	}

	close $fh3 or die;
}

# print out allele 1	
 foreach my $all1 (keys %allele1hash){
	my %a1hash = %{ $allele1hash{$all1} };
	
	open my $fh4, ">>", $outdir1.$all1.'_1.fa';
	
	foreach my $defline (sort(keys %a1hash)) {
		print $fh4 ">".$defline."\n";
		print $fh4 $a1hash{$defline}."\n";
	}

	close $fh4 or die;
}	




sub fasta_to_hash
	{
	my ($fastafh) = @_; 
	my $header;
	my $seq;
	my %sequence;

	while (<$fastafh>)
		{
		chomp($_);
		
		next if $_ =~ /^\s$/;
		next if $_ =~/#/;
 
		# If the line stars for ">" is a header, and we save the information in $header
		if ($_=~/^>(.+)/)
			{
 			# if we have something in $seq
			if ($seq)
				{
				$sequence{$header}=$seq;
				}
 
 
			$header = $1;
	# 		$counter++;
			# reset sequence!!!
			$seq    = '';
			}
		else
			{
			$seq.=$_;
			}
 
		}
	# Store last sequence
	$sequence{$header}=$seq;

	return \%sequence;
	}
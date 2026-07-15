use strict;
use warnings;


##########################################################################################
# Script to make consensus sequences from fasta files with 2 alleles per individual      #
# Input directory should have individual fasta files for each locus                      #
# These files should come from splitStacksFasta.pl script and have deflines like:        #
# >DhoIxt1-22_1204_Allele_0                        										 #
# Set directories before running. Make sure there is a trailing forward slash.           #
##########################################################################################
	
# all loci
# my $indir = '/Volumes/CatDisk/SDZooRAD/2024/phylo/fasta/diploid/';
# my $outdir = '/Volumes/CatDisk/SDZooRAD/2024/phylo/fasta/cons/';

my $indir = '/Volumes/CatDisk/SDZooRAD/2024/phylo/Revision/fasta/diploid/';
my $outdir = '/Volumes/CatDisk/SDZooRAD/2024/phylo/Revision/fasta/cons/';

# Open directory with fasta files of haplotype sequences
opendir (my $dh, $indir) or die;

# Read files, one at a time
while (readdir $dh) {
	next unless $_ =~ /.*\.fa/;
	
# Set outfile name and open
	my $outname = $_;
	$outname =~ s/\.fa/_con.fa/;
	my $outfile = $outdir.$outname;
	open my $ofh, ">", $outfile or die "Can't open $outfile: $!";

# Open infile - locus fasta
# Convert file to hash
	open my $fh, '<', $indir.$_;
	my $fastaref = fasta_to_hash($fh);
	close $fh or die;
	my %fastahash = %{$fastaref};
 
# Go through deflines, capturing Allele 0 only
# Find Allele 1
# Send both sequences to the consensus subroutine
# Print defline and consensus to outfile
 	foreach my $allele (keys(%fastahash)) {
		my @def = split /_Allele_/, $allele;

		if ($def[1] == 0) {
			my $seq0 = $fastahash{$allele};
			my $def1 = $def[0]."_Allele_1";
			my $seq1 = $fastahash{$def1};
			my $conseq = makeCon($seq0,$seq1);
			print $ofh ">".$def[0]."\n";
			print $ofh $conseq."\n";
		}
		else {
			next;
		}
	}
	close $ofh or die;
}


sub makeCon
	{
	my ($seq0,$seq1) = @_;
	my %iupac = ();
	my @con = ();
	my @hets = qw/AC CA  AT TA  AG GA  CT TC  CG GC  GT TG/;
	my @codes = qw/M M  W W  R R  Y Y  S S  K K/;
	@iupac{@hets} = @codes;
		
	my @seqLetters0 = split //, uc($seq0);
	my @seqLetters1 = split //, uc($seq1);
	foreach my $pos (0..scalar(@seqLetters0)-1)
		{
		my $base0 = $seqLetters0[$pos];
		my $base1 = $seqLetters1[$pos];
		
		if ($base0 =~ /[n,N]/ || $base1 =~ /[n,N]/) {
			$con[$pos] = 'N';
			next;
		}
		elsif ($base0 eq $base1)
			{
			$con[$pos] = $base0;
			}
		else
			{
			my $het = $base0.$base1;
			$con[$pos] = $iupac{$het};
			}
		}
	
	my $conseq = join ('', @con);
	return $conseq;
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
 
		# If the line starts with ">" it is a header, and we save the information in $header
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
			$seq = '';
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
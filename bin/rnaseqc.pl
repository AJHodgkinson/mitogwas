use strict;

my $directory = $ARGV[0]; #Read in directory from command line

opendir my $dir, "$directory/RNAseQC/" or die "Cannot open directory: $!";
my @files = readdir $dir;
closedir $dir;

foreach my $file (@files) {
  if ($file =~ /metrics.tsv$/) {
    open (METRICS, "$directory/RNAseQC/$file") || die "Unable to open metrics file $file to read: $!\n";

    my $reject=0;
    while(<METRICS>) {
      chomp;
      my @array=split;
      if ($_=~/^Intergenic Rate/) {
	$reject=1 if ($array[2]>0.3);
      }
      if ($_=~/Base Mismatch/) {
	$reject=1 if ($array[2]>0.01);
      }
      if ($_=~/rRNA Rate/) {
	$reject=1 if ($array[2]>0.3);
      }
      if ($_=~/Total Reads/) {
	$reject=1 if ($array[2]<5000000);
      }
    }

    close (METRICS);

    my $file2=$file;
    $file2=~s/metrics.tsv/gene_reads.gct/;

    my %mtgenes=("MT-RNR1"=>0,"MT-RNR2"=>0,"MT-ND1"=>0,"MT-ND2"=>0,"MT-CO1"=>0,"MT-CO2"=>0,"MT-ATP8"=>0,"MT-ATP6"=>0,"MT-CO3"=>0,"MT-ND3"=>0,"MT-ND4L"=>0,"MT-ND4"=>0,"MT-ND5"=>0,"MT-ND6"=>0,"MT-CYB"=>0); #genes
    open (COUNTS, "$directory/RNAseQC/$file2") || die "Unable to open counts file $file to read: $!\n";
    
    while(<COUNTS>) {
      chomp;
      my @array=split;
      if (exists $mtgenes{$array[1]}) {
	if ($array[2]==0) {
	  $reject=1;
	}
      }
    }

    close (COUNTS);
      
    my $stub=$file; $stub=~s/.Aligned.sortedByCoord.out.PP.bam.metrics.tsv//;
    print "$stub\n" if ($reject==1);
  }
}
  

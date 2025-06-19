use strict;

#This program generates a count matrix from HTseq files for use in R.  Allows us to look at the distributions for normalisation

my $dir = $ARGV[0]; #Directory to get data from
my $remfile = $ARGV[1]; #IDs that fail QC - RNA codes
my $keepfile = $ARGV[2]; #IDs to use to create matrix creating matrix(from infosheet) - RNA codes
my $outfile = $ARGV[3]; #Outfile for matrix

my %ignore = (); #open file and store name list in hash table
open (REM, "$remfile") || die "Unable to open remove file to read: $!\n";
while (<REM>) {
  chomp;
  my @array=split;
  $ignore{$array[0]}=0;
}
close (REM);

my %keep = (); #open info file with selected individuals and store in hash table if they are not in remove list
open (KEEP, "$keepfile") || die "Unable to open keep file to read: $!\n";
while (<KEEP>) {
  chomp;
  my @array=split;
  if (!($_=~/^DNA/)) {
    $keep{$array[1]}=0 if (!(exists $ignore{$array[1]}));
  }
}
close (KEEP);

opendir my $directory, "$dir" or die "Cannot open directory: $!";
my @files = readdir $directory;
closedir $directory;

my @genes = ();
my @inds = ();
my %data = ();
my $total_inds=0;

###Create gene expression data
foreach my $file (@files) {
  if ($file=~/ReadsPerGene.out.tab/) {
    my $id = $file; $id=~s/.ReadsPerGene.out.tab//;
    if (exists $keep{$id}) {
      $total_inds++;  
      push @inds, $id;
      open (F1, "$dir/$file") || die "Unable to open htseq file to read: $!\n";
      while (<F1>) { #Normal File
	if ($_ =~ /^E/) {
	  my @array=split;
	  my $info = $id."_".$array[0];
	  $data{$info}.=$array[1];
	  push @genes, $array[0] if ($total_inds==1);
	}
      }	
      close (F1);
    }
  }
}


###Print out matrix count file for use in R
open (OUT, ">$outfile") || die "Unable to open $outfile to write to: $!\n";

print OUT "\t$_" foreach @inds;
print OUT "\n";

foreach my $gene (@genes) { #For each gene
  print OUT "$gene";
  foreach my $id (@inds) { #for each ID
    my $info = $id."_".$gene;
    print OUT "\t$data{$info}" if (exists $data{$info}); #Print count data in matrix
    print OUT "\t0" if (!(exists $data{$info})); #Zero if doesn't exist for ID
  }
  print OUT "\n";
}

close (OUT);

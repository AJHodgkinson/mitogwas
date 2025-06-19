use strict;

#This program generates a count matrix from HTseq files for use in R.  Allows us to look at the distributions for normalisation

my $dir = $ARGV[0]; #Directory to get data from
my $remfile = $ARGV[1]; #IDs to remove before creating matrix
my $keepfile = $ARGV[2]; #IDs to use to create matrix creating matrix(from infosheet) - RNA codes
my $outfile = $ARGV[3]; #Outfile for matrix
my $gtf_file = $ARGV[4]; #GTF File

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

###Get gene lengths
my %length = ();
my %exon_starts; my %exon_ends=();

open (FILE, "$gtf_file") || die "Unable to open file to read: $!\n";
while (<FILE>) {
  if ($_=~/exon/) {
    if ($_=~/gene_id \"(ENSG\d{1,20}.\d{1,3})\"/) {
      my $gene=$1;
      my @array=split;
      $exon_starts{$gene}.=",".$array[3];
      $exon_ends{$gene}.=",".$array[4];
    }
  }
}  

close (FILE);

while (my ($a,$b) = each %exon_starts) {
  my $starts=$b; $starts=~s/,//;
  my $ends=$exon_ends{$a}; $ends=~s/,//;
  my @starts=split(/\,/,$starts);
  my @ends=split(/\,/,$ends);
  my %sites=();
  for (my $i=0;$i<@starts;$i++) {
    for (my $j=$starts[$i];$j<=$ends[$i];$j++) {
      $sites{$j}=0;
    }
  }
  my $count=0;
  while (my ($a,$b) = each %sites) {
    $count++;
  }
  $length{$a}=($count/1000);
}


###Get Mito Genes
open (MT, "$gtf_file") || die "Unable to open mito genes file: $!\n";
my %mt_genes = ();

while (<MT>) {
  chomp;
  my $ens; my $gene;
  my @array=split;
  if ($array[0]=~/chrM/) {
    if ($array[2] =~ /gene/) {
      if ($_ =~ /gene_name "(.{1,30})"\;/) {
	$gene = $1;
      }
      if ($_ =~ /gene_id "(.{1,30})"\;/) { #Get ensembl ID
	$ens = $1;
      }
    }
    $mt_genes{$ens} = $gene; # if (!($gene =~ /-T/)); #If it is a gene ID line, store MT gene ID (ensembl)
  }
}

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
      print "$id\n";
      $total_inds++;  
      push @inds, $id;      
      open (F1, "$dir/$file") || die "Unable to open htseq file to read: $!\n";
      my @adjust=(); my @names=();
      while (<F1>) { #Normal File
	if ($_ =~ /^E/) {
	  my @array=split;
	  if (exists $mt_genes{$array[0]}) {
	    my $length=$length{$array[0]}; #get length
	    my $perkb=$array[1]/$length; #normalise by length
	    push @names,$array[0]; #store gene name
	    push @adjust,$perkb; #store normalised count
	    push @genes, $array[0] if ($total_inds==1);
	  }
	}
      }
      close (F1);
      my $sum=0;
      $sum+=$_ foreach @adjust; #Get sum of counts
      $sum=$sum/1000000; #Turn into per million
      for (my $i=0;$i<@adjust;$i++) {
	my $final="NA";
	$final=$adjust[$i]/$sum if ($sum>0); #normalise by library size
	my $info = $id."_".$names[$i];
	$data{$info}.=$final if (!($final=~/NA/));
      }
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

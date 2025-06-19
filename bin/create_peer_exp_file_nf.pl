use strict;

my $infile = $ARGV[0]; #Dataset to convert to peer file: usually batch expression file

my $outfile = $infile; $outfile =~ s/\.txt//; $outfile = $outfile.".peer.csv";
my $headerfile = $infile; $headerfile =~ s/\.txt//; $headerfile = $headerfile.".peer.header.txt";

my $tissuez = $infile; $tissuez=~s/log10_counts_all_tpm_//; $tissuez=~s/\.txt//;

#Get gene expression data
my @new = (); my @inds_rna=();
my $count=0;
open (BATCH, "$infile") || die "Unable to open file $infile to read: $!\n";
while (<BATCH>) {
  chomp; my $line=$_;
  if ($count==0) {
    $line="id".$line;
    $count++;
    @inds_rna=split;
    unshift @inds_rna;
  }
  my @array=split(/\s/,$line);
  for (my $i=0;$i<@array;$i++) {
    $new[$i] .= ",".$array[$i];
  }
}

close (BATCH);

my %batch=();
      
open (OUT, ">$outfile") || die "Unable to open file $outfile to write to: $!\n";
open (HEAD, ">$headerfile") || die "Unable to open file $outfile to write to: $!\n";

my $check=0;

print HEAD "\"\"";

foreach my $line (@new) {
  my $release=$line;
  $release=~s/,//;
  my @a=split(/,/,$release);
  my $i=shift(@a);
  my $string=join(",",@a);
  if (!($release=~/ENSG/)) {
    print OUT "$string\n";
    print HEAD "\t$i";
    print "$i\n";
  }
}

print HEAD "\n";

close (OUT);
close (HEAD);




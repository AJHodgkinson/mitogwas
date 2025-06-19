use strict;

my $dir = $ARGV[0]; #Directory to get data from
my $remfile = $ARGV[1]; #IDs to remove before creating matrix
my $keepfile = $ARGV[2]; #IDs to use to create matrix creating matrix(from infosheet) - RNA codes
my $outfile = $ARGV[3]; #Outfile for matrix

my %p9_meth = (585=>0,1610=>0,4271=>0,5520=>0,7526=>0,8303=>0,9999=>0,10413=>0,12146=>0,12274=>0,14734=>0);
my %meth = (585=>0,1610=>0,3238=>0,4271=>0,4392=>0,5520=>0,5647=>0,5721=>0,5818=>0,5883=>0,7526=>0,8303=>0,9999=>0,10413=>0,12146=>0,12274=>0,14734=>0,15896=>0,15948=>0,2617=>0,13710=>0);
my @meth = (585,1610,3238,4271,4392,5520,5647,5721,5818,5883,7526,8303,9999,10413,12146,12274,14734,15896,15948,2617,13710);
my %data = ();
my @inds = ();

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

foreach my $file (@files) { #For all files in directory
  if ($file=~/MTcountfile\.23\.txt/) { #Get MT countfiles
    my $id = $file; $id=~s/.MTcountfile.23.txt//;
    if (exists $keep{$id}) { #Removes unwanted inds as specified above
      foreach my $site (@meth) { #Create hash for different sites/inds/resampling
	my $tag = $id."_".$site;
	$data{$tag."_full"} = "NA";
	$data{$tag."_s10"} = "NA";
	$data{$tag."_s20"} = "NA";
	$data{$tag."_s30"} = "NA";
	$data{$tag."_s50"} = "NA";
	$data{$tag."_s100"} = "NA";
      }
      print "$id\n";  
      push @inds, $id;
      my $total_count=0; my $total_alt=0; my $total_p9_count10=0; my $total_p9_alt10=0; my $total_p9_count20=0; my $total_p9_alt20=0;
      open (F1, "$dir/$file") || die "Unable to open htseq file to read: $!\n"; #For count file
      while (<F1>) { #Normal File
	chomp;
	my @array=split;
	if (exists $meth{$array[1]}) { #For sites we expect to see methylation (Science paper 2014)
	  my $site = $array[1];
	  my $tag = $id."_".$site;
	  my ($res,$alt,$tot)=&MethFull($_,20); #Subroutine to calc total alternative alleles over all (if >=20 depth)
	  $total_count+=$tot if (exists $p9_meth{$array[1]});
	  $total_alt+=$alt if (exists $p9_meth{$array[1]});
	  my $resample10=&MethResampleMulti($_,10); #Subroutine to randomly sample 10X
	  my $resample20=&MethResampleMulti($_,20); #Subroutine to randomly sample 20X
	  if ($resample10=~/\d/) { #This is just to generate overall meth level at 10X across all P9 sites
	    my $comb_alt=$resample10*10;
	    $total_p9_count10+=10 if (exists $p9_meth{$array[1]});
	    $total_p9_alt10+=$comb_alt if (exists $p9_meth{$array[1]});
	  }
	  if ($resample20=~/\d/) { #This is just to generate overall meth level at 20X across all P9 sites
	    my $comb_alt=$resample20*20;
	    $total_p9_count20+=20 if (exists $p9_meth{$array[1]});
	    $total_p9_alt20+=$comb_alt if (exists $p9_meth{$array[1]});
	  }
	  my $resample30=&MethResampleMulti($_,30);  #Subroutine to randomly sample 30X
	  my $resample50=&MethResampleMulti($_,50);  #Subroutine to randomly sample 50X
	  my $resample100=&MethResampleMulti($_,100);  #Subroutine to randomly sample 100X
	  $data{$tag."_full"} = $res;
	  $data{$tag."_s10"} = $resample10;
	  $data{$tag."_s20"} = $resample20;
	  $data{$tag."_s30"} = $resample30;
	  $data{$tag."_s50"} = $resample50;
	  $data{$tag."_s100"} = $resample100;
	}
      }
      close (F1);
      #Calculate comibined average for p9 sites:
      my $combined="NA";
      $combined=$total_alt/$total_count if ($total_count>0);
      $data{$id."_combined_full"} = $combined;
      my $combined_resample10="NA";
      $combined_resample10=$total_p9_alt10/$total_p9_count10 if ($total_p9_count10==(11*10));
      $data{$id."_combined_s10"} = $combined_resample10;
      my $combined_resample20="NA";
      $combined_resample20=$total_p9_alt20/$total_p9_count20 if ($total_p9_count20==(11*20));
      $data{$id."_combined_s20"} = $combined_resample20;
    }
  }
}
	  
#Create Plink file with data - rows are IDs, Columns are different ways to present methylation level
open (OUT, ">$outfile") || die "Unable to open file to write to: $!\n"; #Now create pheno file for plink

print OUT "FID\tIID";
foreach my $site (@meth) {
  print OUT "\t$site\_full\t$site\_s10\t$site\_s20\t$site\_s30\t$site\_s50\t$site\_s100";
}
print OUT "\tcombined_full\tcombined_s10\tcombined_s20\n";

foreach my $id (@inds) {
  print OUT "$id\t$id";
  my $sum=0; my $count=0; my $sum_alt=0; my $sum_tot=0;
  foreach my $site (@meth) {
    my $tag = $id."_".$site;
    foreach my $extra ("_full","_s10","_s20","_s30","_s50","_s100") {
      my $final=$tag.$extra;
      print OUT "\t$data{$final}"; #Print meth rates
    }
  }
  my $fullcomb=$id."_combined_full"; my $rescomb10=$id."_combined_s10"; my $rescomb20=$id."_combined_s20";
  print OUT "\t$data{$fullcomb}\t$data{$rescomb10}\t$data{$rescomb20}\n"; #report resampled and full combined values
}

close (OUT);

################################

sub MethFull ($$) { #Calculate full mthylation level - alt alleles over total allele
  my ($line,$d) = @_;
  my $res="NA";
  my $tot=0;my $alt=0;
  my @array = split(/\s/,$line);
  my @sides = split (/:/,$array[3]);
  my @nucs = split (/,/,$sides[0]);
  my @counts = split (/,/,$sides[1]);
  my @check; push @check, $counts[$_] for (0..7);
  my $sum = eval join '+', @check;
  if ($sum>=$d) {
    for (my $i=0;$i<8;$i++) {
      $tot+=$counts[$i];
      my $n = uc($nucs[$i]);
      if ($n ne $array[2]) {
	$alt+=$counts[$i];
      }
    }
    $res = $alt/$tot;
  }
  return ($res,$alt,$tot);
}

sub MethResampleMulti ($$) { #Resample to desired depth
  my ($line,$c) = @_;
  my $res_rand;
  my @array = split(/\s/,$line);
  my @sides = split (/:/,$array[3]);
  my @nucs = split (/,/,$sides[0]);
  my @counts = split (/,/,$sides[1]);
  my @check; push @check, $counts[$_] for (0..7);
  my $sum = eval join '+', @check;
  if ($sum<$c) {
    $res_rand = "NA";
  }
  if ($sum>=$c) {
    my @r_counts=(0,0,0,0,0,0,0,0);
    my @all_counts;
    for (my $i=0;$i<8;$i++) {
      push @all_counts, $i for (1..$counts[$i]);
    }
    my %got = ();
    for (my $i=0;$i<$c;$i++) {
    redo_rand:
      my $rand = int(rand(@all_counts));
      if (exists $got{$rand}) {
	goto redo_rand;
      }
      if (!(exists $got{$rand})) {
	$r_counts[$all_counts[$rand]]++;
	$got{$rand}=0;
      }
      my $tot_rand=0;my $alt_rand=0;
      for (my $i=0;$i<8;$i++) {
	$tot_rand+=$r_counts[$i];
	my $n = uc($nucs[$i]);
	if ($n ne $array[2]) {
	  $alt_rand+=$r_counts[$i];
	}
      }
      $res_rand = $alt_rand/$tot_rand;
    }
  }
  return ($res_rand);
}

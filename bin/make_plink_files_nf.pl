use strict;
use Statistics::R;

my $R = Statistics::R->new();

my $dataset = $ARGV[0]; #Dataset to convert
my $meth_file = $ARGV[1]; # "meth_".$dataset.".txt";
my $exp_log_mtnuc_file = $ARGV[2]; #"log10_counts_mtnuc_tpm_".$dataset.".txt";
my $exp_log_mtmt_file = $ARGV[3]; # "log10_counts_mtmt_tpm_".$dataset.".txt";
my $exp_med_mtnuc_file = $ARGV[4]; #"median_log10_counts_mtnuc_tpm_".$dataset.".txt";
my $exp_med_mtmt_file = $ARGV[5]; #"median_log10_counts_mtmt_tpm_".$dataset.".txt";
my $exp_quant_mtnuc_file = $ARGV[6]; #"quantile_log10_counts_mtnuc_tpm_".$dataset.".txt";
my $exp_quant_mtmt_file = $ARGV[7]; #"quantile_log10_counts_mtmt_tpm_".$dataset.".txt";
my $infosheet = $ARGV[8]; #Infosheet with required DNA:RNA match and covariates
my $peerfile =  $ARGV[9]; #location of peer file
my $peerfileheader =  $ARGV[10]; #location of peer file header
my $gpcfile =  $ARGV[11]; #location of gpc file

my $out_exp =  "plink_pheno_".$dataset.".txt";
my $out_exp_unmask =  "plink_pheno_".$dataset."_unmasked.txt";

####Open infosheet and store RNA to DNA code conversion, and then store covar information in a hash table as ID_COVAR=VALUE
##List of covars in @covars
my %convert=(); my %covars=(); my @covars=(); my @header=(); my $line=1;
open (DATA, "$infosheet") || die "Unable to open file to read: $!\n";
while (<DATA>) {
  chomp;
  my @array=split;
  $convert{$array[1]}=$array[0]; #store data as RNAID=DNAID
  if ($line==1) {
    @header=@array;
    for (my $i=2;$i<@array;$i++) {
      push @covars, $array[$i];
    }
  }
  if ($line>1) {
    for (my $i=2;$i<@array;$i++) {
      my $tag=$array[1]."_".$header[$i];
      $covars{$tag}=$array[$i]; #store data as RNAID_COVAR=VALUE
    }
  }
  $line++;
}
close (DATA);

#Use subroutine to open each expression file and store in expression hash table as TYPE_GENE_ID=VALUE
my %expression = (); my %genes = (); my %inds = (); my %expression_unmasked=();
&MakeHash($exp_log_mtnuc_file,"log_mtnuc");
&MakeHash($exp_log_mtmt_file,"log_mtmt");
&MakeHash($exp_med_mtnuc_file,"med_mtnuc");
&MakeHash($exp_med_mtmt_file,"med_mtmt");
&MakeHash($exp_quant_mtnuc_file,"quant_mtnuc");
&MakeHash($exp_quant_mtmt_file,"quant_mtmt");

#Store IDs in @inds and Gene names in @genes
my @genes=(); my @inds=(); my @states = ("log_mtnuc","log_mtmt","med_mtnuc","med_mtmt","quant_mtnuc","quant_mtmt");
while (my ($a,$b) = each %genes) {
  push @genes, $a;
}
while (my ($a,$b) = each %inds) {
  push @inds, $a;
}

#Open methylation file and store in meth hash as ID_SITE=VALUE
my %meth=(); my @meth=(); my $line=1; @header=(); my %meth_unmasked=();
open (METH, "$meth_file") || die "Unable to open file to read: $!\n";
while (<METH>) {
  chomp;
  my @array=split;
  if ($line==1) {
    @header=@array;
    for (my $i=2;$i<@array;$i++) {
      push @meth, $array[$i];
    }
  }
  if ($line>1) {
    for (my $i=2;$i<@array;$i++) {
      my $tag=$array[0]."_".$header[$i];
      $meth{$tag}=$array[$i]; #store data as RNAID_METH=VALUE
      $meth_unmasked{$tag}=$array[$i]; #store data as RNAID_METH=VALUE
    }
  }
  $line++;
}
close (METH);


##Open PEER file and store as ID_PEERX=VALUE
@header = ();
my $count=0;
open (HEAD, "$peerfileheader") || die "Unable to open file $peerfileheader here to read: $!\n";
while (<HEAD>) {
  chomp;
  if ($count==0) {
    @header = split;
    shift @header;
  }
  $count++;
}
close (HEAD);

my %peer_info = (); my $factor=0;
open (PEER, "$peerfile") || die "Unable to open file $peerfile to read: $!\n";
while (<PEER>) {
  chomp;
  $factor++;
  my @array=split(/\,/,$_);
  for (my $i=0;$i<@array;$i++) {
    my $tag=$header[$i]."_PEER".$factor;
    $peer_info{$tag}=$array[$i];
  }
}
close (PEER);

##Open GPC file and store as ID_GPCX=VALUE
my %gpc_info = ();
open (PCA, "$gpcfile") || die "Unable to open file $gpcfile to read: $!\n";
while (<PCA>) {
  chomp;
  if (!($_=~/#/)) {
    my @array=split;
    my @a=split(/\:/,$array[0]);
    for (my $i=1;$i<11;$i++) {
      my $tag=$a[0]."_GPC".$i;
      $gpc_info{$tag}=$array[$i];
    }
  }
}

close (PCA);


##########MASK outlier values with NA

foreach my $g (@genes) {
  foreach my $s (@states) {
    my $string;
    foreach my $i (@inds) {
      my $tag=$s."_".$g."_".$i;
      my $level = "NA";
      $level = $expression{$tag} if (exists $expression{$tag});
      $string .=",".$level;
    }
    $string=~s/,//;
    my $make= "dat<-c($string)\nlow=quantile(dat,na.rm=T)[2]\nmed=quantile(dat,na.rm=T)[3]\nhigh=quantile(dat,na.rm=T)[4]";
    $R->send($make);
    my $low = $R->get('low');
    my $med = $R->get('med');
    my $high = $R->get('high');
    $low=~s/\s//g; $med=~s/\s//g; $high=~s/\s//g;
    $low=~s/0\%//g; $med=~s/50\%//g; $high=~s/100\%//g;
    my $range=($high-$low)*3; #Outlier defined as 3 interquartile ranges above/below 3rd and 1st quartile respectively
    my $top=$high+$range;
    my $bot=$low-$range;
    foreach my $i (@inds) {
      my $tag=$s."_".$g."_".$i;
      my $level = "NA";
      $level = $expression{$tag} if (exists $expression{$tag});
      $expression{$tag}="NA" if (($level<$bot)||($level>$top));
    }
  }
}

foreach my $sites (@meth) {
  my $string;
  foreach my $i (@inds) {
    my $tag=$i."_".$sites;
    my $level = "NA";
    $level = $meth{$tag} if (exists $meth{$tag});
    $string .=",".$level;
  }
  $string=~s/,//;
  my $make= "dat<-c($string)\nlow=quantile(dat,na.rm=T)[2]\nmed=quantile(dat,na.rm=T)[3]\nhigh=quantile(dat,na.rm=T)[4]";
  $R->send($make);
  my $low = $R->get('low');
  my $med = $R->get('med');
  my $high = $R->get('high');
  $low=~s/\s//g; $med=~s/\s//g; $high=~s/\s//g;
  $low=~s/0\%//g; $med=~s/50\%//g; $high=~s/100\%//g;
  my $range=($high-$low)*3; #Outlier defined as 3 interquartile ranges above/below 3rd and 1st quartile respectively
  my $top=$high+$range;
  my $bot=$low-$range;
  foreach my $i (@inds) {
    my $tag=$i."_".$sites;
    my $level = "NA";
    $level = $meth{$tag} if (exists $meth{$tag});
    $meth{$tag}="NA" if (($level<$bot)||($level>$top));
  }
}


#########CREATE boxplots:

my $exp_string; my $gene_string;
foreach my $g (@genes) {
  foreach my $i (@inds) {
    my $tag="log_mtnuc_".$g."_".$i;
    my $level = "NA";
    $level = $expression{$tag} if (exists $expression{$tag});
    $exp_string .=",".$level;
    $gene_string .=",\"".$g."\"";
  }
}
$exp_string=~s/,//; $gene_string=~s/,//;
my $make="data<-c($exp_string)\ngenes<-c($gene_string)";
print "$make\n";
$R->send($make);
my $make="pdf(\"${dataset}.expression.pdf\",height = 8, width = 12)\npar(mar=c(12,5,1,1))\npar(mgp = c(3.5, 1, 0))\nboxplot(data~genes,ylab=\"Gene Expression Log10(TPM)\",xlab=NULL,cex.axis=1,cex.names=1.5,las=2,col=\"light blue\",whisklty = 1,whisklwd = 2,staplelwd = 2,outpch = 16, outcex = 0.5,cex.lab=1.5,frame.plot=FALSE)\nmtext(\"Mitochondrial Genes\", side=1, line=10, cex=1.5)\ndev.off()";
$R->send($make);

my $mod_string; my $site_string;
foreach my $s (@meth) {
  foreach my $i (@inds) {
    my $tag=$i."_".$s;
    if ($s=~/full/) {
      my $level = "NA";
      $level = $meth{$tag} if (exists $meth{$tag});
      $mod_string .=",".$level;
      $site_string .=",\"".$s."\"";
    }
  }
}
$mod_string=~s/,//; $site_string=~s/,//;
my $make="data1<-c($mod_string)\nsites<-c($site_string)";
print "$make\n";
$R->send($make);
my $make="pdf(\"${dataset}.modification.pdf\",height = 8, width = 12)\npar(mar=c(7,5,1,1))\npar(mgp = c(3.5, 1, 0))\nboxplot(data1~sites,ylab=\"RNA Modification Levels\",xlab=NULL,cex.axis=1,cex.names=1.5,las=2,col=\"light blue\",whisklty = 1,whisklwd = 2,staplelwd = 2,outpch = 16, outcex = 0.5,cex.lab=1.5,frame.plot=FALSE)\nmtext(\"Location on Mitochondrial Genome\", side=1, line=6, cex=1.5)\ndev.off()";
$R->send($make);



#######CREATE PLINK FILE with all data and covars

open (OUT1, ">$out_exp") || die "Unable to open file $out_exp to write to: $!\n";

print OUT1 "FID\tIID";

foreach my $g (@genes) {
  foreach my $s (@states) {
    my $tag = $g."_".$s;
    print OUT1 "\t$tag";
  }
}

foreach my $g (@meth) {
  print OUT1 "\t$g";
}

foreach my $g (@covars) {
  print OUT1 "\t$g";
}

print OUT1 "\tPEER1\tPEER2\tPEER3\tPEER4\tPEER5\tPEER6\tPEER7\tPEER8\tPEER9\tPEER10\tPEER11\tPEER12\tPEER13\tPEER14\tPEER15\tPEER16\tPEER17\tPEER18\tPEER19\tPEER20\tGPC1\tGPC2\tGPC3\tGPC4\tGPC5\tGPC6\tGPC7\tGPC8\tGPC9\tGPC10\tRNA\n";

foreach my $i (@inds) {
  my $dna=$convert{$i};
  my $rna=$i;
  print OUT1 "$dna\t$dna";
  foreach my $g (@genes) {
    foreach my $s (@states) {
      my $tag = $s."_".$g."_".$i;
      my $level = "NA";
      $level = $expression{$tag} if (exists $expression{$tag});
      print OUT1 "\t$level";
    }
  }
  foreach my $sites (@meth) {
    my $tag=$i."_".$sites;
    my $level = "NA";
    $level = $meth{$tag} if (exists $meth{$tag});
    print OUT1 "\t$level";
  }
  foreach my $cov (@covars) {
    my $tag=$i."_".$cov;
    my $level = "NA";
    $level = $covars{$tag} if (exists $covars{$tag});
    print OUT1 "\t$level";
  }
  for (my $j=1;$j<=20;$j++) {
    my $tag=$i."_PEER".$j;
    my $level = "NA";
    $level = $peer_info{$tag} if (exists $peer_info{$tag});
    print OUT1 "\t$level";
  }
  for (my $j=1;$j<=10;$j++) {
    my $tag=$dna."_GPC".$j;
    my $level = "NA";
    $level = $gpc_info{$tag} if (exists $gpc_info{$tag});
    print OUT1 "\t$level";
  }    
  print OUT1 "\t$rna\n";
}

close (OUT1);

#######CREATE PLINK FILE with all data and covars (UNMASKED)

open (OUT2, ">$out_exp_unmask") || die "Unable to open file $out_exp to write to: $!\n";

print OUT2 "FID\tIID";

foreach my $g (@genes) {
  foreach my $s (@states) {
    my $tag = $g."_".$s;
    print OUT2 "\t$tag";
  }
}

foreach my $g (@meth) {
  print OUT2 "\t$g";
}

foreach my $g (@covars) {
  print OUT2 "\t$g";
}

print OUT2 "\tPEER1\tPEER2\tPEER3\tPEER4\tPEER5\tPEER6\tPEER7\tPEER8\tPEER9\tPEER10\tPEER11\tPEER12\tPEER13\tPEER14\tPEER15\tPEER16\tPEER17\tPEER18\tPEER19\tPEER20\tGPC1\tGPC2\tGPC3\tGPC4\tGPC5\tGPC6\tGPC7\tGPC8\tGPC9\tGPC10\tRNA\n";

foreach my $i (@inds) {
  my $dna=$convert{$i};
  my $rna=$i;
  print OUT2 "$dna\t$dna";
  foreach my $g (@genes) {
    foreach my $s (@states) {
      my $tag = $s."_".$g."_".$i;
      my $level = "NA";
      $level = $expression_unmasked{$tag} if (exists $expression_unmasked{$tag});
      print OUT2 "\t$level";
    }
  }
  foreach my $sites (@meth) {
    my $tag=$i."_".$sites;
    my $level = "NA";
    $level = $meth_unmasked{$tag} if (exists $meth_unmasked{$tag});
    print OUT2 "\t$level";
  }
  foreach my $cov (@covars) {
    my $tag=$i."_".$cov;
    my $level = "NA";
    $level = $covars{$tag} if (exists $covars{$tag});
    print OUT2 "\t$level";
  }
  for (my $j=1;$j<=20;$j++) {
    my $tag=$i."_PEER".$j;
    my $level = "NA";
    $level = $peer_info{$tag} if (exists $peer_info{$tag});
    print OUT2 "\t$level";
  }
  for (my $j=1;$j<=10;$j++) {
    my $tag=$dna."_GPC".$j;
    my $level = "NA";
    $level = $gpc_info{$tag} if (exists $gpc_info{$tag});
    print OUT2 "\t$level";
  }    
  print OUT2 "\t$rna\n";
}

close (OUT2);

#####

sub MakeHash ($$) {
  my ($f,$t) = @_;
  open (FILE, "$f") || die "Unable to open $exp_log_mtnuc_file to read: $!\n";
  my @head = (); 
  while (<FILE>) {
    chomp;
    if (!($_=~/ENS/)) {
      @head=split;
      my @new=();
      foreach my $item (@head) {
	my $hh=$item; $hh=~s/\"//g;
	push @new, $hh;
      }
      @head=@new;
      unshift(@head, 'EMPTY');
    }
    if ($_=~/ENS/) {
      my @array=split;
      my $gene = $array[0];
      $gene =~ s/\"//g;
      $genes{$gene}=0;
      for (my $i=1;$i<@array;$i++) {
	my $tag = $t."_".$gene."_".$head[$i];
	$inds{$head[$i]}=0;
	$expression{$tag}=$array[$i];
	$expression_unmasked{$tag}=$array[$i];
      }
    } 
  }
}

close (FILE1);

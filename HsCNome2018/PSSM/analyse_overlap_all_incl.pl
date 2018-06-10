#!/usr/bin/perl
use strict;

open(fr1,"PosRefSet_w_expt_proof_Jan2017")||print "cannot open HsCN_ints_HIPPIE $!\n";   # Reading PosRefSet
my @prs=<fr1>;chomp(@prs);

print "Enter the in silico prediction file (path & filename)\n";my $file=<STDIN>;chomp($file); # Reading predictions
open(fr2,$file)||print "cannot open $file $!\n";
my @pssm=<fr2>;chomp(@pssm);

open(fr3,"SWATH_n_DDA_WT_BFDR_le0.05_count_matrix_w_log2")||print "cannot open input $!\n";   # Reading BioID
my @bioid=<fr3>;chomp(@bioid);
print "Enter the Log2(WT:Mut) ratio cut-off if you wish, (Enter na if you don't want filtering)\t";my $cut=<STDIN>;chomp($cut);if($cut=~/na/i){$cut="NA";}

my @splitted_line;my %prs;my %pssm;my %bioid;my %grand;my $i;

foreach my $line(@prs)
{
 if($line!~/^Acc/)
 {
  @splitted_line=split("\t",$line);
  $prs{$splitted_line[0]}="";
  $grand{$splitted_line[0]}="";
 }
}
foreach my $line(@pssm)
{
 if($line!~/^pssm_p/)
 {
  @splitted_line=split("\t",$line);
  $pssm{$splitted_line[5]}="";
  $grand{$splitted_line[5]}="";
 }
}
for($i=3;$i<scalar @bioid;$i++)
{
 @splitted_line=split("\t",$bioid[$i]);
 if($cut eq "NA"){$bioid{$splitted_line[1]}="";$grand{$splitted_line[1]}="";}
 else
 {
  if(($splitted_line[11]>$cut)||($splitted_line[12]>$cut))
  {
   $bioid{$splitted_line[1]}="";
   $grand{$splitted_line[1]}="";
  }
 }
}
my $prs=keys %prs;my $pssm=keys %pssm;my $bioid=keys %bioid;
print "$prs proteins in PosRefSet, $pssm proteins in PSSM and $bioid proteins in BioID\n";

#### Find protein common to all ####
my %all;
foreach my $acc1(sort keys %prs)
{
 foreach my $acc2(sort keys %pssm)
 {
  if($acc1 eq $acc2)
  {
   foreach my $acc3(sort keys %bioid)
   {
    if($acc1 eq $acc3)
    {
     $all{$acc1}="";
    }
   }
  }
 }
}
my $all=keys %all;print "$all proteins common to all 3 - PosRefSet, PSSM and BioID (WTCN w BFDR<=0.05)\n";
print "Do you want their UP acc?\t[y/n]: ";my $choice=<STDIN>;chomp($choice);$choice=uc($choice);
if($choice eq "Y"){&print(\%all);}
######################################

#### proteins common to PosRefSet and BIOID ####
my $nomatch1=0;my $nomatch2=0;my %prs_bioid;
foreach my $acc1(sort keys %prs)
{
 foreach my $acc(sort keys %all){if($acc1 ne $acc){$nomatch1++;}}
 foreach my $acc2(sort keys %pssm){if($acc1 ne $acc2){$nomatch2++;}}
 if(($nomatch1 eq $all)&&($nomatch2 eq $pssm))
 {
  foreach my $acc3(sort keys %bioid)
  {
   if($acc1 eq $acc3)
   {
    $prs_bioid{$acc1}="";
   }
  }
 }
 $nomatch1=0;$nomatch2=0;
}
my $prs_bioid=keys %prs_bioid;print "$prs_bioid proteins common to PosRefSet and BioID (WTCN w BFDR<=0.05)\n";
print "Do you want their UP acc?\t[y/n]: ";my $choice=<STDIN>;chomp($choice);$choice=uc($choice);
if($choice eq "Y"){&print(\%prs_bioid);}
################################################

### proteins common to PosRefSet and PSSM ####
my $nomatch1=0;my $nomatch2=0;my %prs_pssm;
foreach my $acc1(sort keys %prs)
{
 foreach my $acc(sort keys %all){if($acc1 ne $acc){$nomatch1++;}}
 foreach my $acc3(sort keys %bioid){if($acc1 ne $acc3){$nomatch2++;}}
 if(($nomatch1 eq $all)&&($nomatch2 eq $bioid))
 {
  foreach my $acc2(sort keys %pssm)
  {
   if($acc1 eq $acc2)
   {
    $prs_pssm{$acc1}="";
   }
  }
 }
 $nomatch1=0;$nomatch2=0;
}
my $prs_pssm=keys %prs_pssm;print "$prs_pssm proteins common to PosRefSet and PSSM\n";
print "Do you want their UP acc?\t[y/n]: ";my $choice=<STDIN>;chomp($choice);$choice=uc($choice);
if($choice eq "Y"){&print(\%prs_pssm);}
################################################

### proteins common to BioID and PSSM ####
my $nomatch1=0;my $nomatch2=0;my %bioid_pssm;
foreach my $acc3(sort keys %bioid)
{
 foreach my $acc(sort keys %all){if($acc3 ne $acc){$nomatch1++;}}
 foreach my $acc1(sort keys %prs){if($acc1 ne $acc3){$nomatch2++;}}
 if(($nomatch1 eq $all)&&($nomatch2 eq $prs))
 {
  foreach my $acc2(sort keys %pssm)
  {
   if($acc3 eq $acc2)
   {
    $bioid_pssm{$acc3}="";
   }
  }
 }
 $nomatch1=0;$nomatch2=0;
}
my $bioid_pssm=keys %bioid_pssm;print "$bioid_pssm proteins common to BioID and PSSM\n";
print "Do you want their UP acc?\t[y/n]: ";my $choice=<STDIN>;chomp($choice);$choice=uc($choice);
if($choice eq "Y"){&print(\%bioid_pssm);}
################################################

my $prs_only=$prs-($all+$prs_pssm+$prs_bioid);
my $pssm_only=$pssm-($all+$prs_pssm+$bioid_pssm);
my $bioid_only=$bioid-($all+$prs_bioid+$bioid_pssm);
my $grand=keys %grand;
print "$prs_only proteins remain in PosRefSet\n$pssm_only proteins remain in PSSM AND\n$bioid_only proteins remain in BioID\nIf you take Union of these three sets, you have $grand proteins\n";
print "Do you want their UP acc?\t[y/n]: ";my $choice=<STDIN>;chomp($choice);$choice=uc($choice);
if($choice eq "Y"){&print(\%grand);}

sub print()
{
 my %hash=%{$_[0]};
 foreach my $acc(sort keys %hash)
 {
  print "$acc\n";
 }
}

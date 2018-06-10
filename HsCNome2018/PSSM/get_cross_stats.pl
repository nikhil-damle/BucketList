# This prog asks things like how many prots w Px predicted wo any filters also have a Lx predicted w diff filters? OR
# how many Lx predicted using accessibility and IUPRED filter also have a Px predicted with highest confidence?

#!/usr/bin/perl
use strict;

my $dir="/HsCNome_dataRepository";

open(fr1,"$dir/top5K_hits")||print "cannot open top5K_hits $!\n";
my @data1=<fr1>;chomp(@data1);
open(fr2,"$dir/accessible_SLiMs_from_top5K_each")||print "cannot open accessible_SLiMs_from_top5K_each $!\n";
my @data2=<fr2>;chomp(@data2);
open(fr3,"$dir/PxIU_me0.3_LxIU_me0.35")||print "cannot open PxIU_me0.3_LxIU_me0.35 $!\n";
my @data3=<fr3>;chomp(@data3);
open(fr4,"$dir/acc_SLiMs_PxIU_ge0.3_pval_le0.00006_n_LxIU_ge0.35_pval_le0.00008")||print "cannot open acc_SLiMs_PxIU_ge0.3_pval_le0.00006_n_LxIU_ge0.35_pval_le0.00008 $!\n";
my @data4=<fr4>;chomp(@data4);

my %px1;my %px2;my %px3;my %px4;my %lx1;my %lx2;my %lx3;my %lx4;
my %common11Px;my %common12Px;my %common13Px;my %common14Px;my %common21Px;my %common22Px;my %common23Px;my %common24Px;my %common31Px;my %common32Px;my %common33Px;my %common34Px;my %common41Px;my %common42Px;my %common43Px;my %common44Px;
my %common11Lx;my %common12Lx;my %common13Lx;my %common14Lx;my %common21Lx;my %common22Lx;my %common23Lx;my %common24Lx;my %common31Lx;my %common32Lx;my %common33Lx;my %common34Lx;my %common41Lx;my %common42Lx;my %common43Lx;my %common44Lx;

print "\nProcessing proteins w predicted PxIxITs first\n-----------------------------------------\n";
%px1=&getPx(\@data1);my $px1=keys %px1;print "No filters fetched $px1 prots\n";
%px2=&getPx(\@data2);my $px2=keys %px2;print "Accessibillity filter fetched $px2 prots\n";
%px3=&getPx(\@data3);my $px3=keys %px3;print "CN-acc and IUPRED filters fetched $px3 prots\n";
%px4=&getPx(\@data4);my $px4=keys %px4;print "CN-acc, IU and pval filters fetched $px4 prots\n";

print "\nNow processing proteins w predicted LxVPs\n-----------------------------------------------\n";
%lx1=&getLx(\@data1);my $lx1=keys %lx1;print "No filters fetched $lx1 prots\n";
%lx2=&getLx(\@data2);my $lx2=keys %lx2;print "CN-accessibility filter fetched $lx2 prots\n";
%lx3=&getLx(\@data3);my $lx3=keys %lx3;print "CN-acc & IUPRED fetched $lx3 prots\n";
%lx4=&getLx(\@data4);my $lx4=keys %lx4;print "CN-acc, IU and pval filters fetched $lx4 prots\n\n";

%common11Px=&compare(\%px1,\%lx1);out(\%common11Px);
%common12Px=&compare(\%px1,\%lx2);out(\%common12Px);
%common13Px=&compare(\%px1,\%lx3);out(\%common13Px);
%common14Px=&compare(\%px1,\%lx4);out(\%common14Px);
%common21Px=&compare(\%px2,\%lx1);out(\%common21Px);
%common22Px=&compare(\%px2,\%lx2);out(\%common22Px);
%common23Px=&compare(\%px2,\%lx3);out(\%common23Px);
%common24Px=&compare(\%px2,\%lx4);out(\%common24Px);
%common31Px=&compare(\%px3,\%lx1);out(\%common31Px);
%common32Px=&compare(\%px3,\%lx2);out(\%common32Px);
%common33Px=&compare(\%px3,\%lx3);out(\%common33Px);
%common34Px=&compare(\%px3,\%lx4);out(\%common34Px);
%common41Px=&compare(\%px4,\%lx1);out(\%common41Px);
%common42Px=&compare(\%px4,\%lx2);out(\%common42Px);
%common43Px=&compare(\%px4,\%lx3);out(\%common43Px);
%common44Px=&compare(\%px4,\%lx4);out(\%common44Px);

%common11Lx=&compare(\%lx1,\%px1);out(\%common11Lx);
%common12Lx=&compare(\%lx1,\%px2);out(\%common12Lx);
%common13Lx=&compare(\%lx1,\%px3);out(\%common13Lx);
%common14Lx=&compare(\%lx1,\%px4);out(\%common14Lx);
%common21Lx=&compare(\%lx2,\%px1);out(\%common21Lx);
%common22Lx=&compare(\%lx2,\%px2);out(\%common22Lx);
%common23Lx=&compare(\%lx2,\%px3);out(\%common23Lx);
%common24Lx=&compare(\%lx2,\%px4);out(\%common24Lx);
%common31Lx=&compare(\%lx3,\%px1);out(\%common31Lx);
%common32Lx=&compare(\%lx3,\%px2);out(\%common32Lx);
%common33Lx=&compare(\%lx3,\%px3);out(\%common33Lx);
%common34Lx=&compare(\%lx3,\%px4);out(\%common34Lx);
%common41Lx=&compare(\%lx4,\%px1);out(\%common41Lx);
%common42Lx=&compare(\%lx4,\%px2);out(\%common42Lx);
%common43Lx=&compare(\%lx4,\%px3);out(\%common43Lx);
%common44Lx=&compare(\%lx4,\%px4);out(\%common44Lx);

sub getPx()
{
 my @data=@{$_[0]};my @splitted_line;my %px;my $len;
 foreach my $line(@data)
 {
  if($line!~/^pssm_p/)
  {
   @splitted_line=split("\t",$line);
   if($splitted_line[12] eq "PxIxIT")
   {
    $px{$splitted_line[5]}="";
   }
  }
 }
 return %px;
}

sub getLx()
{
 my @data=@{$_[0]};my @splitted_line;my %lx;
 foreach my $line(@data)
 {
  if($line!~/^pssm_p/)
  {
   @splitted_line=split("\t",$line);
   if($splitted_line[12] eq "LxVP")
   {
    $lx{$splitted_line[5]}="";
   }
  }
 }
 return %lx;
}

sub compare()
{
 my %hash1=%{$_[0]};my %hash2=%{$_[1]};my %common;
 my $len1=keys %hash1;my $len2=keys %hash2;print "You gave $len1 prots in Pxset1 and $len2 prots in Lxset\t";

 foreach my $acc1(sort keys %hash1)
 {
  foreach my $acc2(sort keys %hash2)
  {
   if($acc1 eq $acc2)
   {
    $common{$acc1}="";
   }
  }
 }
 return %common;
}

sub out()
{
 my %hash=%{$_[0]};
 my $prots=keys %hash;
 print "$prots proteins in common between the two sets\n";
}

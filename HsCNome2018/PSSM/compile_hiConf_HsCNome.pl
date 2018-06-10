#!usr/bin/perl
use strict;

open(prs,"/PosRefSet_w_expt_proof_Jan2017")||print "cannot open PosRefSet_w_expt_proof_Jan2017 $!\n";
my @prs=<prs>;chomp(@prs);
open(propdPx,"PxIxIT_file.tab")||print "cannot open PxIxIT_file.tab $!\n";
my @propd_Px=<propdPx>;chomp(@propd_Px);
open(propdLx,"LxVP_file.tab")||print "cannot open LxVP_file.tab $!\n";
my @propd_Lx=<propdLx>;chomp(@propd_Lx);
open(pred,"acc_SLiMs_PxIU_ge0.3_pval_le0.00006_n_LxIU_ge0.35_pval_le0.00008")||print "cannot open predictions $!\n";
my @preds=<pred>;chomp(@preds);

open(bioid,"SWATH_n_DDA_WT_BFDR_le0.05_count_matrix_w_log2_ratio_gt0.5")||print "cannot open SWATH_n_DDA_WT_BFDR_le0.05_count_matrix_w_log2_ratio_gt0.5 $!\n";
my @bioid=<bioid>;chomp(@bioid);
open(cnacc,"accessible_SLiMs_from_top5K_each")||print "cannot open accessible_SLiMs_from_top5K_each $!\n";
my @cnacc=<cnacc>;chomp(@cnacc);

my @splitted_line;my $i;my %UPs;my %PRS;my %proPx;my %proLx;my %preds;my %pxpred;my %lxpred;my %cnacc;my %bioid;my %cnacc_in_bioid;my @motifarr;my %propd;
my $noprs=0;my $nopropd=0;my $nopred=0;my $nomatch=0;my $pxflag=0;my $lxflag=0;

for($i=1;$i<scalar @prs;$i++)
{
 @splitted_line=split("\t",$prs[$i]);
 $PRS{$splitted_line[0]}="";
 $UPs{$splitted_line[0]}="";
}
my $prs=keys %PRS;print "PosRefSet had $prs proteins\n";

foreach my $line(@propd_Px)
{
 @splitted_line=split("\t",$line);
 $proPx{$splitted_line[1]}{"PxIxIT"}="";
 $UPs{$splitted_line[1]}="";
 $propd{$splitted_line[1]}="";
}
foreach my $line(@propd_Lx)
{
 @splitted_line=split("\t",$line);
 $proLx{$splitted_line[1]}{"LxVP"}="";
 $UPs{$splitted_line[1]}="";
 $propd{$splitted_line[1]}="";
}
my $propd=keys %propd;
my $proPx=keys %proPx;
my $proLx=keys %proLx;
print "CN-binding SLiMs from $propd proteins went into training\n";
print "Px from $proPx prots went into Px training\n";
print "Lx from $proLx prots went into Lx training\n";

for($i=1;$i<scalar @preds;$i++)
{
 @splitted_line=split("\t",$preds[$i]);
 $preds{$splitted_line[5]}="";
 if($splitted_line[12] eq "PxIxIT")
 {
  $pxpred{$splitted_line[5]}{$splitted_line[12]}="";
  $UPs{$splitted_line[5]}="";
 }
 if($splitted_line[12] eq "LxVP")
 {
  $lxpred{$splitted_line[5]}{$splitted_line[12]}="";
  $UPs{$splitted_line[5]}="";
 }
}
my $pred=keys %preds;
my $pxpred=keys %pxpred;
my $lxpred=keys %lxpred;
print "CN-binding SLiMs predicted in $pred proteins of which $pxpred have Px and $lxpred have Lx\n";

for($i=1;$i<scalar @cnacc;$i++)
{
 @splitted_line=split("\t",$cnacc[$i]);
 $cnacc{$splitted_line[5]}{$splitted_line[12]}="";
}
my $cnacc=keys %cnacc;print "CNacc SLiMs predicted in $cnacc proteins\n";
foreach my $line(@bioid)
{
 @splitted_line=split("\t",$line);
 if($splitted_line[1] ne "")
 {
  $bioid{$splitted_line[1]}="";
 }
}
my $bioid=keys %bioid;print "Log2(WT:Mut)>=0.5 for $bioid proteins in BioID\n";
foreach my $acc1(sort keys %bioid)
{
 foreach my $acc2(sort keys %cnacc)
 {
  if($acc1 eq $acc2)
  {
   $cnacc_in_bioid{$acc1}{$cnacc{$acc2}}="";
   $UPs{$acc1}="";
  }
 }
}
my $cnacc_in_bioid=keys %cnacc_in_bioid;print "CNacc SLiMs predicted in $cnacc_in_bioid proteins suh that Log2(WT:Mut)>=0.5 in BioID\n";
my $hscnome=keys %UPs;print "$hscnome prots in HsCNome\n";
print "UPacc\tPosRefSet\tProPD_train\thighConfPred\tBioID_pred\tPxIxIT\tLxVP\n";
foreach my $acc1(sort keys %UPs)
{
 print "$acc1\t";
 foreach my $acc2(sort keys %PRS)
 {
  if($acc1 eq $acc2){print "Y\t";}
  else
  {
   $noprs++;
   if($noprs eq $prs){print "N\t";}
  }
 }
 foreach my $acc2(sort sort keys %propd)
 {
  if($acc1 eq $acc2)
  {
   print "Y\t";
   foreach my $motif(sort keys %{$proPx{$acc1}})
   {
    if($motif eq "PxIxIT"){$pxflag=1;}
   }
   foreach my $motif(sort keys %{$proLx{$acc1}})
   {
    if($motif eq "LxVP"){$lxflag=1;}
   }
  }
  else
  {
   $nopropd++;
   if($nopropd eq $propd){print "N\t";}
  }
 }
 foreach my $acc2(sort sort keys %preds)
 {
  if($acc1 eq $acc2)
  {
   print "Y\t";
   foreach my $motif(sort keys %{$pxpred{$acc1}})
   {
    if($motif eq "PxIxIT"){$pxflag=1;}
   }
   foreach my $motif(sort keys %{$lxpred{$acc1}})
   {
    if($motif eq "LxVP"){$lxflag=1;}
   }
  }
  else
  {
   $nopred++;
   if($nopred eq $pred){print "N\t";}
  }
 }
 foreach my $acc2(sort sort keys %cnacc_in_bioid)
 {
  if($acc1 eq $acc2)
  {
   print "Y\t";
   foreach my $motif(sort keys %{$cnacc_in_bioid{$acc1}})
   {
    if($motif eq "PxIxIT"){$pxflag=1;}
    if($motif eq "LxVP"){$lxflag=1;}
   }
  }
  else
  {
   $nomatch++;
   if($nomatch eq $cnacc_in_bioid){print "N\t";}
  }
 }
 if($pxflag eq 1){print "Y\t";}else{print "N\t";}
 if($lxflag eq 1){print "Y\n";}else{print "N\n";}
 $noprs=0;$nopropd=0;$nopred=0;$nomatch=0;$pxflag=0;$lxflag=0;
}

#!/usr/bin/perl
use strict;

print "Enter the working directory [SWATH/DDA]\t";my $dir=<STDIN>;chomp($dir);$dir=uc($dir);
open(wt,"$dir/CN_WT_irresp_o_BFDR")||print "cannot opne WT $!\n";
my @wt=<wt>;chomp(@wt);
open(nir,"$dir/CN_NIR_irresp_o_BFDR")||print "cannot opne NIR $!\n";
my @nir=<nir>;chomp(@nir);
open(wf,"$dir/CN_WF_irresp_o_BFDR")||print "cannot opne WF $!\n";
my @wf=<wf>;chomp(@wf);

my %all;my $nomatch=0;my $count=0;my %complete;my $acc;my $gene;
my @wt_counts;my @nir_counts;my @wf_counts;my $wt_count;my $nir_count;my $wf_count;my @splitted_line;

open(fw,">>$dir/count_matrix_$dir\_irresp_o_BFDR");print fw "Prot_acc\tGene\tWT\tWT_BFDR\tNIR\tNIR_BFDR\tWF\tWF_BFDR\n";

my %wt=&get_hash(\@wt,"WT");
my $wt=keys %wt;print "WT has $wt proteins\n";
my %nir=&get_hash(\@nir,"NIR");
my $nir=keys %nir;print "NIR has $nir proteins\n";
my %wf=&get_hash(\@wf,"WF");
my $wf=keys %wf;print "WF has $wf proteins\n";

my $complete=keys %complete;print "$complete prots found int w Cn in different expts\n";

foreach my $str(sort keys %complete)
{
 $str=~/(.*)NIKHIL(.*)/;
 $acc=$1;$gene=$2;print fw "$acc\t$gene\t";
 foreach my $line(@wt)
 {
  @splitted_line=split("\t",$line);
  if($acc eq $splitted_line[1])
  {
   print fw "$splitted_line[5]\t$splitted_line[15]\t";
  }
  else
  {
   $nomatch++;
   if($nomatch eq scalar @wt)
   {
    print fw "0\tNA\t";$nomatch=0;
   }
  }
 }
 $nomatch=0;
 foreach my $line(@nir)
 {
  @splitted_line=split("\t",$line);
  if($acc eq $splitted_line[1])
  {
   print fw "$splitted_line[5]\t$splitted_line[15]\t";
  }
  else
  {
   $nomatch++;
   if($nomatch eq scalar @nir)
   {
    print fw "0\tNA\t";$nomatch=0;
   }
  }
 }
 $nomatch=0;
 foreach my $line(@wf)
 {
  @splitted_line=split("\t",$line);
  if($acc eq $splitted_line[1])
  {
   print fw "$splitted_line[5]\t$splitted_line[15]\n";
  }
  else
  {
   $nomatch++;
   if($nomatch eq scalar @wf)
   {
    print fw "0\tNA\n";$nomatch=0;
   }
  }
 }
 $nomatch=0;#exit();
}

foreach my $acc(sort keys %complete)
{
 if($wt{$acc}){@wt_counts=@{$wt{$acc}};}#print "$acc\t@wt_counts\n";exit();
 if($nir{$acc}){@nir_counts=@{$nir{$acc}};}
 if($wf{$acc}){@wf_counts=@{$wf{$acc}};}
 
 if(scalar @wt_counts >= 1){$wt_count=&get_avgCount(\@wt_counts);}else{$wt_count=0;}
 if(scalar @nir_counts >= 1){$nir_count=&get_avgCount(\@nir_counts);}else{$nir_count=0;}
 if(scalar @wf_counts >= 1){$wf_count=&get_avgCount(\@wf_counts);}else{$wf_count=0;}
 #print fw "$acc\t$complete{$acc}\t$wt_count\t$nir_count\t$wf_count\n";
 undef @wt_counts;undef @nir_counts; undef @wf_counts;
}
sub get_avgCount()
{
 my @counts=@{$_[0]};my $num=@counts;my $total_count=0;my $avg;
 foreach my $val(@counts)
 {
  $total_count=$total_count+$val;
 }
 $avg=$total_count/$num;
 return $avg;#print "$acc\t@wt_counts\t$avg\n";exit();
}

foreach my $acc1(sort keys %wt)
{
 foreach my $acc2(sort keys %nir)
 {
  if($acc1 eq $acc2)
  {
   foreach my $acc3(sort keys %wf)
   {
    if($acc1 eq $acc3)
    {
     $all{$acc1}="";
    }
   }
  }
 }
}
my $all=keys %all;print "$all proteins found to interact w all forms of CN w BFDR <= 0.05\n";

get_uniq(\%wt,\%nir,\%wf,\%all,"WT");
get_uniq(\%nir,\%wt,\%wf,\%all,"NIR");
get_uniq(\%wf,\%wt,\%nir,\%all,"WF");

sub get_uniq()
{
 my %hash1=%{$_[0]};my %hash2=%{$_[1]};my %hash3=%{$_[2]};my %all=%{$_[3]};my $description=$_[4];my $nomatch=0;my $nomatch2=0;my $nomatch3=0;my %uniq;
 my $hash2=keys %hash2;my $hash3=keys %hash3;my $all=keys %all;
 foreach my $acc1(sort keys %hash1)
 {
  foreach my $acc(sort keys %all)
  {
   if($acc1 ne $acc)
   {
    $nomatch++;
    if($nomatch eq $all)
    {
     foreach my $acc2(sort keys %hash2)
     {
      if($acc1 ne $acc2)
      {
       $nomatch2++;
       if($nomatch2 eq $hash2)
       {
        foreach my $acc3(sort keys %hash3)
        {
         if($acc1 ne $acc3)
         {
          $nomatch3++;
          if($nomatch3 eq $hash3)
          {
           $uniq{$acc1}="";
          }
         }
        }
        $nomatch3=0;
       }
      }
     }
     $nomatch2=0;
    }
   }
  }
  $nomatch=0;
 }
 my $uniq=keys %uniq;print "$uniq proteins are uniqly found in $description expt\n";
 #foreach my $acc(sort keys %uniq)
 #{
 # print "$acc\n";
 #}
}
print "**********************************************\n";
compare(\%wt,\%nir,\%all,"WT-NIR");
compare(\%wt,\%wf,\%all,"WT-WF");
compare(\%nir,\%wf,\%all,"NIR-WF");

sub compare()
{
 my %hash1=%{$_[0]};my %hash2=%{$_[1]};my %all=%{$_[2]};my $description=$_[3];my $nomatch=0;
 my $all=keys %all;
 foreach my $acc1(sort keys %hash1)
 {
  foreach my $acc(sort keys %all)
  {
   if($acc1 ne $acc)
   {
    $nomatch++;
    if($nomatch eq $all)
    {
     foreach my $acc2(sort keys %hash2)
     {
      if($acc1 eq $acc2)
      {
       #print "$acc1 is common b/w $description\n";
       $count++;
      }
     }
    }
   }
  }
  $nomatch=0;
 }
 print "$count prots are common b/w only $description\n";
 $count=0;
}

sub get_hash()
{
 my @data=@{$_[0]};my @splitted_line;my %hash;my $str;my $type=$_[1];print "Checking $type\n";
 foreach my $line(@data)
 {
  if($line!~/^Bait/)
  {
   @splitted_line=split("\t",$line);
   $str=$splitted_line[1]."NIKHIL".$splitted_line[2];
   $complete{$str}="";
   if($splitted_line[5] ne 0)
   {
    push(@{$hash{$splitted_line[1]}},$splitted_line[5]);#="";
    #$complete{$splitted_line[1]}++;
   }
  }
 }
 return %hash;
}

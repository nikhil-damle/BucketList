#!/usr/bin/perl
use strict;

open(fr1,"/Hs_disorderome_pepPD_Ylva/CNWT_ipoCa/CNWT_ipoCa")||print "cannot open CNWT_ipoCa $!\n";
my @data1=<fr1>;chomp(@data1);my $len1=@data1;
open(fr2,"/Hs_disorderome_pepPD_Ylva/CNWT_iaoCa/CNWT_iaoCa")||print "cannot open CNWT_ipoCa $!\n";
my @data2=<fr2>;chomp(@data2);my $len2=@data2;
open(fr3,"/Hs_disorderome_pepPD_Ylva/CNNIR_ipoCa/CNNIR_ipoCa")||print "cannot open CNWT_ipoCa $!\n";
my @data3=<fr3>;chomp(@data3);my $len3=@data3;
open(fr4,"/Hs_disorderome_pepPD_Ylva/CNNIR_iaoCa/CNNIR_iaoCa")||print "cannot open CNWT_ipoCa $!\n";
my @data4=<fr4>;chomp(@data4);my $len4=@data4;
open(fr5,"/Hs_disorderome_pepPD_Ylva/CNWF_ipoCa/CNWF_ipoCa")||print "cannot open CNWT_ipoCa $!\n";
my @data5=<fr5>;chomp(@data5);my $len5=@data5;
open(fr6,"/Hs_disorderome_pepPD_Ylva/CNWF_iaoCa/CNWF_iaoCa")||print "cannot open CNWT_ipoCa $!\n";
my @data6=<fr6>;chomp(@data6);my $len6=@data6;

open(fw,">>/Hs_disorderome_pepPD_Ylva/all_counts_12Sept");

my @splitted_line;my %all;my $pep_counter=0;my $prot_counter=0;my $nomatch1=0;my $nomatch2=0;my $nomatch3=0;my $nomatch4=0;my $nomatch5=0;my $nomatch6=0;
my $tot_pep1=0;my $tot_pep2=0;my $tot_pep3=0;my $tot_pep4=0;my $tot_pep5=0;my $tot_pep6=0;my $percent1;my $percent2;my $percent3;my $percent4;my $percent5;my $percent6;

foreach my $line(@data1)
{
 if($line!~/^Peptide.*/)
 {
  @splitted_line=split("\t",$line);
  $all{$splitted_line[2]}{$splitted_line[1]}{$splitted_line[0]}="";
  $tot_pep1=$tot_pep1+$splitted_line[3];
 }
}
foreach my $line(@data2)
{
 if($line!~/^Peptide.*/)
 {
  @splitted_line=split("\t",$line);
  $all{$splitted_line[2]}{$splitted_line[1]}{$splitted_line[0]}="";
  $tot_pep2=$tot_pep2+$splitted_line[3];
 }
}
foreach my $line(@data3)
{
 if($line!~/^Peptide.*/)
 {
  @splitted_line=split("\t",$line);
  $all{$splitted_line[2]}{$splitted_line[1]}{$splitted_line[0]}="";
  $tot_pep3=$tot_pep3+$splitted_line[3];
 }
}
foreach my $line(@data4)
{
 if($line!~/^Peptide.*/)
 {
  @splitted_line=split("\t",$line);
  $all{$splitted_line[2]}{$splitted_line[1]}{$splitted_line[0]}="";
  $tot_pep4=$tot_pep4+$splitted_line[3];
 }
}
foreach my $line(@data5)
{
 if($line!~/^Peptide.*/)
 {
  @splitted_line=split("\t",$line);
  $all{$splitted_line[2]}{$splitted_line[1]}{$splitted_line[0]}="";
  $tot_pep5=$tot_pep5+$splitted_line[3];
 }
}
foreach my $line(@data6)
{
 if($line!~/^Peptide.*/)
 {
  @splitted_line=split("\t",$line);
  $all{$splitted_line[2]}{$splitted_line[1]}{$splitted_line[0]}="";
  $tot_pep6=$tot_pep6+$splitted_line[3];
 }
}
my $genes=keys %all;

foreach my $gene(sort keys %all)
{
 foreach my $prot(sort keys %{$all{$gene}})
 {
  $prot_counter++;
  foreach my $pep(sort keys %{$all{$gene}{$prot}})
  {
   $pep_counter++;
   print fw "$gene\t$prot\t$pep\t";

   foreach my $line(@data1)
   {
    @splitted_line=split("\t",$line);
    if(($splitted_line[2] eq $gene)&&($splitted_line[1] eq $prot)&&($splitted_line[0] eq $pep))
    {
     $percent1=($splitted_line[3]/$tot_pep1)*100;$percent1=sprintf("%.3f",$percent1);
     #print fw "$percent1\t";
     print fw "$splitted_line[6]\t";
    }
    else
    {
     $nomatch1++;
     if($nomatch1 eq $len1)
     {
      print fw "0\t";
     }
    }
   }
   $nomatch1=0;

   foreach my $line(@data2)
   {
    @splitted_line=split("\t",$line);
    if(($splitted_line[2] eq $gene)&&($splitted_line[1] eq $prot)&&($splitted_line[0] eq $pep))
    {
     $percent2=($splitted_line[3]/$tot_pep2)*100;$percent2=sprintf("%.3f",$percent2);
     #print fw "$percent2\t";
     print fw "$splitted_line[6]\t";
    }
    else
    {
     $nomatch2++;
     if($nomatch2 eq $len2)
     {
      print fw "0\t";
     }
    }
   }
   $nomatch2=0;

   foreach my $line(@data3)
   {
    @splitted_line=split("\t",$line);
    if(($splitted_line[2] eq $gene)&&($splitted_line[1] eq $prot)&&($splitted_line[0] eq $pep))
    {
     $percent3=($splitted_line[3]/$tot_pep3)*100;$percent3=sprintf("%.3f",$percent3);
     #print fw "$percent3\t";
     print fw "$splitted_line[3]\t";
    }
    else
    {
     $nomatch3++;
     if($nomatch3 eq $len3)
     {
      print fw "0\t";
     }
    }
   }
   $nomatch3=0;

   foreach my $line(@data4)
   {
    @splitted_line=split("\t",$line);
    if(($splitted_line[2] eq $gene)&&($splitted_line[1] eq $prot)&&($splitted_line[0] eq $pep))
    {
     $percent4=($splitted_line[3]/$tot_pep4)*100;$percent4=sprintf("%.3f",$percent4);
     #print fw "$percent4\t";
     print fw "$splitted_line[3]\t";
    }
    else
    {
     $nomatch4++;
     if($nomatch4 eq $len4)
     {
      print fw "0\t";
     }
    }
   }
   $nomatch4=0;

   foreach my $line(@data5)
   {
    @splitted_line=split("\t",$line);
    if(($splitted_line[2] eq $gene)&&($splitted_line[1] eq $prot)&&($splitted_line[0] eq $pep))
    {
     $percent5=($splitted_line[3]/$tot_pep5)*100;$percent5=sprintf("%.3f",$percent5);
     #print fw "$percent5\t";
     print fw "$splitted_line[3]\t";
    }
    else
    {
     $nomatch5++;
     if($nomatch5 eq $len5)
     {
      print fw "0\t";
     }
    }
   }
   $nomatch5=0;

   foreach my $line(@data6)
   {
    @splitted_line=split("\t",$line);
    if(($splitted_line[2] eq $gene)&&($splitted_line[1] eq $prot)&&($splitted_line[0] eq $pep))
    {
     $percent6=($splitted_line[3]/$tot_pep6)*100;$percent6=sprintf("%.3f",$percent6);
     #print fw "$percent6\n";
     print fw "$splitted_line[3]\n";
    }
    else
    {
     $nomatch6++;
     if($nomatch6 eq $len6)
     {
      print fw "0\n";
     }
    }
   }
   $nomatch6=0;
  }
 }
}
print "You just checked $pep_counter peps from $genes genes encoding $prot_counter proteins\n";

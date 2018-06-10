#!/usr/bin/perl
use strict;

my $dir="/HsCNome_analyses_Feb2018/GO_analyses/BPs";
#my $dir="/HsCNome_analyses_Feb2018/GO_analyses/MFs";
#my $dir="/HsCNome_analyses_Feb2018/GO_analyses/CCs";

open(fr,"$dir/gene_term_associations")||print "cannot open gene_BP_association $!\n";
my @data=<fr>;chomp(@data);

my @splitted_line;my %uniq;my $num;my %genes;my %pairs;my %multi;

print "Which set of genes do you want?\ngenes associated with only one term (print 1)\ngenes associated with exactly 2 terms (print 2)\ngenes associated with >2 terms (print x)\n";
my $inp=<STDIN>;chomp($inp);my $inp=uc($inp);
if(($inp ne 1)&&($inp ne 2)&&($inp ne "X")){print "You must obey my orders\n";exit();}

foreach my $line(@data)
{
 @splitted_line=split("\t",$line);
 $genes{$splitted_line[0]}="";
 if(($splitted_line[2] eq 1)&&($inp eq 1))
 {
  push(@{$uniq{$splitted_line[3]}},$splitted_line[1]);
 }
 if(($splitted_line[2] eq 2)&&($inp eq 2))
 {
  push(@{$pairs{$splitted_line[3]}},$splitted_line[1]);
  #print "$line\n";
 }
 if(($splitted_line[2]>2)&&($inp eq "X"))
 {
  push(@{$multi{$splitted_line[3]}},$splitted_line[1]);
  #print "$line\n";
 }
}
my $genes=keys %genes;
print "\nYou have $genes genes in your daatset\n";

if($inp eq 1){out(\%uniq,$inp,$dir);}
if($inp eq 2){out(\%pairs,$inp,$dir);}
if($inp eq "X"){out(\%multi,$inp,$dir);}

sub out()
{
 my %hash=%{$_[0]};my $p1;my $p2;my $str;my $inp=$_[1];my $dir=$_[2];

 foreach my $BP(sort keys %hash)
 {
  if(($BP=~/(.*)\s(.*)/)&&($inp eq 2))
  {
   $p1=$1;$p2=$2;$str=$p1."_".$p2;
   open(fw,">>$dir/$str\_paired_nodes");
  }
  elsif($inp eq 1)
  {
   open(fw,">>$dir/$BP\_uniq_nodes");
  }
  else
  {
   open(fw,">>$dir/$BP\_nodes");
  }
  $num=@{$hash{$BP}};
  print "$BP\t$num\n";#@{$uniq{$BP}}\n";
  foreach my $gene(sort @{$hash{$BP}})
  {print "$gene\n";
   print fw "$gene\n";
  }
 }
}

#!/usr/bin/perl
use strict;

print "Enter the PSSM predictions\t";my $file=<STDIN>;chomp($file);
open(fr,$file)||print "cannot open $file $!\n";
my @data=<fr>;chomp(@data);
system("rm -f PSSM_len_vs_num_o_SLiMs");
open(fw,">>PSSM_len_vs_num_o_SLiMs");
my $dir="features_combined";   # A dir w $acc.fasta files from HsProteome

my @splitted_line;my $len;my %hash;my @data2;

foreach my $line(@data)
{
 if($line!~/^pssm_p/)
 {
  @splitted_line=split("\t",$line);
  $hash{$splitted_line[5]}++;       # $splitted_line[5]=protein acc
 }
}
my $prots=keys %hash;print "PSSM predicted SLiMs in $prots proteins\n";

foreach my $acc(sort keys %hash)
{
 open(fr2,"$dir/$acc\.fasta")||print "cannot open $acc\.fasta $!\n";
 @data2=<fr2>;chomp(@data2);
 $len=length($data2[1]);
 print fw "$acc\t$len\t$hash{$acc}\n";
}

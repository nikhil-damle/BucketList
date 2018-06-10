#!/usr/bin/perl
use strict;

print "Give me an entire path to the list of acc #\n";my $infile=<STDIN>;chomp($infile);
open(fr,$infile)||print "cannot open $infile $!\n";
my @data=<fr>;chomp(@data);

my %prots;my $seq="";my $totPx=0;my $totLx=0;my $aa_content;my $len;

foreach my $line(@data)
{
 $line=~s/\r//g;$line=~s/\n//g;
 $prots{$line}="";
}
my $prots=keys %prots;print "You gave me $prots prots\n";

foreach my $acc(sort keys %prots)
{
 open(fr,"/features_combined/$acc\.fasta")||print "canot open $acc\.fasta $!\n";
 @data=<fr>;chomp(@data);

 $seq=$data[1];$aa_content=$aa_content.$data[1];
 $len=length($seq);
 $totPx=$totPx+length($seq)-12+1;
 $totLx=$totLx+length($seq)-7+1;
 $seq="";
}
my $aa=length($aa_content);
my $prots=keys %prots;print "The i/p has $prots proteins w $aa aa\n";
print "Possible 12mers: $totPx\n";
print "Possible 7mers: $totLx\n";

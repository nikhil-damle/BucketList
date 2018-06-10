#!/usr/bin/perl
use strict;

print "Enter the acclist of entire proteome\t";my $proteome_acc=<STDIN>;chomp($proteome_acc);
print "Enter the acclist of proteins of interest\t";my $acclist=<STDIN>;chomp($acclist);
open(fr2,$proteome_acc)||print "cannot open proteome acclist $!\n";
my @proteome=<fr2>;chomp(@proteome);

open(fr,$acclist)||print "cannot open PosRefSet_acclist $!\n";
my @prs=<fr>;chomp(@prs);

system("rm -f Ro_entireP_acclist");
open(fw,">>Ro_entireP_acclist");

my $acc;my %proteome;my %prs;my $nomatch=0;my %rest;

foreach my $line(@proteome)
{
 $line=~s/\n//g;$line=~s/\r//g;
 {$proteome{$line}="";}
}
my $proteome=keys %proteome;print "Proteome has $proteome proteins\n";

foreach my $line(@prs)
{
 $line=~s/\n//g;$line=~s/\r//g;
 $prs{$line}="";
}
my $prs=keys %prs;print "$prs proteinsof interest\n";

foreach my $acc1(sort keys %proteome)
{
 foreach my $acc2(sort keys %prs)
 {
  if($acc1 ne $acc2)
  {
   $nomatch++;
   if($nomatch eq $prs)
   {
    print fw "$acc1\n";
    $rest{$acc1}="";
   }
  }
 }
 $nomatch=0;
}
my $rest=keys %rest;print "Rest of the proteome has $rest proteins\n";

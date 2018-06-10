#!/usr/bin/perl
use strict;

open(fw,">>/data/nikhil/Desktop/Stanford_PDF/HsCN_ints/HsCN_abcr1r2_ints_copy");

open(fr,"/data/nikhil/Desktop/Stanford_PDF/HsCN_ints/BIOGRID-CN_ints.tab2.txt")||print "cannot open BioGRID $!\n";
my @data=<fr>;chomp(@data);

my @splitted_line;my %CN;my $i;my %method;my %geneIDs;my $counter=0;my $depod_only;my $partner;my $pubmed;my $acc;
my %BioGRID;my %depod;my %intact;my %hippie;

for($i=1;$i<scalar @data;$i++)
{
 @splitted_line=split("\t",$data[$i]);
 if(($splitted_line[7] eq "PPP3CA")||($splitted_line[7] eq "PPP3CB")||($splitted_line[7] eq "PPP3CC")||($splitted_line[7] eq "PPP3R1")||($splitted_line[7] eq "PPP3R2"))
 {
  $CN{$splitted_line[7]}="";
  #$geneIDs{$splitted_line[8]}{$splitted_line[14]}="";
  push(@{$geneIDs{$splitted_line[8]}},$splitted_line[14]);
  $BioGRID{$splitted_line[8]}{$splitted_line[14]}="";
 }
 if(($splitted_line[8] eq "PPP3CA")||($splitted_line[8] eq "PPP3CB")||($splitted_line[8] eq "PPP3CC")||($splitted_line[8] eq "PPP3R1")||($splitted_line[8] eq "PPP3R2"))
 {
  $CN{$splitted_line[8]}="";
  #$geneIDs{$splitted_line[7]}{$splitted_line[14]}="";
  push(@{$geneIDs{$splitted_line[7]}},$splitted_line[14]);
  $BioGRID{$splitted_line[7]}{$splitted_line[14]}="";
 }
 $method{$splitted_line[11]}="";
}
my $CN=keys %CN;
my $biogrids=keys %BioGRID;
my $methods=keys %method;
print "$biogrids interactors of CN ($CN geneID) using $methods methods exist in BioGRID\n";

foreach my $meth(keys %method)
{
 print "$meth\n";
}
#exit();
#foreach my $CN(sort keys %CN){print "********\n$CN\n********\n";}

open(fr2,"/data/nikhil/Desktop/Stanford_PDF/HsCN_ints/DEPOD_PP2Balpha")||print "cannot open DEPOD $!\n";
my @data2=<fr2>;chomp(@data2);

foreach my $line(@data2)
{
 @splitted_line=split("\t",$line);
 if($splitted_line[3]=~/Homo sapiens/)
 {
  #$geneIDs{$splitted_line[1]}{$splitted_line[6]}="";#print "$splitted_line[1]\n";
  push(@{$geneIDs{$splitted_line[1]}},$splitted_line[6]);
  $depod{$splitted_line[1]}{$splitted_line[6]}="";
  $counter++;
 }
}
my $depods=keys %depod;print "$depods ints exist in DEPOD\n";

open(fr3,"/data/nikhil/Desktop/Stanford_PDF/HsCN_ints/IntAct.txt")||print "cannot open IntAct $!\n";
my @data3=<fr3>;chomp(@data3);

foreach my $line(@data3)
{
 @splitted_line=split("\t",$line);
 if(($splitted_line[0] eq "uniprotkb:Q08209")||($splitted_line[0] eq "uniprotkb:P16298")||($splitted_line[0] eq "uniprotkb:P48454")||($splitted_line[0] eq "uniprotkb:P63098")||($splitted_line[0] eq "uniprotkb:Q96LZ3"))
 {
  $splitted_line[1]=~/uniprotkb:(.*)/;$acc=$1;
  $splitted_line[5]=~/uniprotkb:(.*)\(gene name\)\|/;$partner=$1;
  $splitted_line[8]=~/pubmed:(\d*)/;my $pubmed=$1;
  $intact{$acc}{$pubmed}="";
  #$geneIDs{$partner}{$pubmed}="";#print "Partner = $partner\t$pubmed\n";exit();
  push(@{$geneIDs{$partner}},$pubmed);
 }
 if(($splitted_line[1] eq "uniprotkb:Q08209")||($splitted_line[1] eq "uniprotkb:P16298")||($splitted_line[1] eq "uniprotkb:P48454")||($splitted_line[1] eq "uniprotkb:P63098")||($splitted_line[1] eq "uniprotkb:Q96LZ3"))
 {
  $splitted_line[0]=~/uniprotkb:(.*)/;$acc=$1;
  $splitted_line[4]=~/uniprotkb:(.*)\(gene name\)\|/;$partner=$1;
  $splitted_line[8]=~/pubmed:(\d*)/;my $pubmed=$1;#print "$splitted_line[8]\n";exit();
  $intact{$acc}{$pubmed}="";
  #$geneIDs{$partner}{$pubmed}="";#print "Partner = $partner\t$pubmed\n";exit();
  push(@{$geneIDs{$partner}},$pubmed);
 }
}
my $intacts=keys %intact;print "$intacts ints exist in IntAct\n";

open(fr4,"/data/nikhil/Desktop/Stanford_PDF/HsCN_ints/HIPPIE_Sept2K15.txt")||print "cannot open HIPPIE_Sept2K15.txt $!\n";
my @data4=<fr4>;chomp(@data4);

foreach my $line(@data4)
{
 @splitted_line=split("\t",$line);
 if(($splitted_line[0] eq "PP2BA_HUMAN")||($splitted_line[0] eq "PP2BB_HUMAN")||($splitted_line[0] eq "PP2BC_HUMAN")||($splitted_line[0] eq "CANB1_HUMAN")||($splitted_line[0] eq "CANB2_HUMAN"))
 {
  $splitted_line[2]=~/(.*)\_HUMAN/;$partner=$1;
  $hippie{$partner}="";
  #$geneIDs{$partner}{$splitted_line[5]}="";
  push(@{$geneIDs{$partner}},$splitted_line[5]);
 }
 if(($splitted_line[2] eq "PP2BA_HUMAN")||($splitted_line[2] eq "PP2BB_HUMAN")||($splitted_line[2] eq "PP2BC_HUMAN")||($splitted_line[2] eq "CANB1_HUMAN")||($splitted_line[2] eq "CANB2_HUMAN"))
 {
  $splitted_line[0]=~/(.*)\_HUMAN/;$partner=$1;
  $hippie{$partner}="";
  #$geneIDs{$partner}{$splitted_line[5]}="";
  push(@{$geneIDs{$partner}},$splitted_line[5]);
 }
}
my $hippie=keys %hippie;print "$hippie ints exist in Sept 2015 release of HIPPIE\n";
my $tot=keys %geneIDs;
print "Which partners do you want to print (type biogrid/depod/intact/hippie/all/none)\t";my $input=<STDIN>;chomp($input);

if($input eq "biogrid"){out(\%BioGRID);}
if($input eq "depod"){out(\%depod);}
if($input eq "intact"){out(\%intact);}
if($input eq "hippie"){out(\%hippie);}
if($input eq "all"){out(\%geneIDs);}
if($input eq "none"){exit();}

print "You compiled $tot unique human interactors of CN\n";


sub out()
{
 my %hash=%{$_[0]};my $val=keys %hash;
 foreach my $ints(sort keys %hash)
 {
  print fw "$ints\t@{$hash{$ints}}\n";
  #print fw "$ints\t";
  #foreach my $pubmed(keys %{$hash{$ints}})
  #{
  # print fw "$pubmed, ";
  #}
  #print fw "\n";
 }
 print "You compiled $val unique human interactors of CN in $input\n";
}

=headopen(fr,"/data/nikhil/Desktop/Stanford_PDF/HsCN_ints/HsCN_abcr1r2_ints_corrected.tab")||print "cannot open HsCN_abcr1r2_ints_corrected.tab $!\n";
my @data=<fr>;chomp(@data);

my @splitted_line;my %UPs;my %prots;my $counter=0;
foreach my $line(@data)
{
 if($line!~/^Entry/)
 {
  @splitted_line=split("\t",$line);
  $UPs{$splitted_line[0]}++;
  $prots{$splitted_line[1]}++;
 }
}
my $UPs=keys %UPs;
my $prots=keys %prots;
print "$UPs acc numbers of $prots proteins in your HsCN interactome\n";
=cut
#foreach my $acc(sort {$UPs{$b}<=>$UPs{$a}} keys %UPs)
#{
# print "$acc\t$UPs{$acc}\n";
# $counter++;
# if($counter eq 10)
# {
#  exit();
# }
#}

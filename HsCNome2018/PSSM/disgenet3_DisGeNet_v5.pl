#### This prog computes stat enrichment of dis-associated genes among predicted Cn-partners #######
# Step1: Map all UP/SP IDs onto Entrez geneIDs - BioMart @ Ensembl genome browser v90, dataset - Hs- GRCh-38.p10 - file: refer "fr2"
# Run clean_martExport.pl before running this prog

#!usr/bin/perl
use strict;

open(fr2,"HsCNome_UPacc_NCBIGeneID_mapping")||print "cannot open HsCNome_UPacc_NCBIGeneID_mapping $!\n";
my @data=<fr2>;chomp(@data);

my %UPs;my %GIs;my %unmapped;my %UP_tots;my %GI_tots;
my @splitted_line;my %GIDs;my %prots;my @ids;my @acc;

foreach my $line(@data)
{
 if($line!~/^UniProtKB/)
 {
  @splitted_line=split("\t",$line);
  if($splitted_line[0] ne ""){$UP_tots{$splitted_line[0]}="";}
  if($splitted_line[1] ne ""){$GI_tots{$splitted_line[1]}="";}
  if(($splitted_line[0] ne "")&&($splitted_line[1] ne ""))
  {
   $UPs{$splitted_line[0]}{$splitted_line[1]}="";
   $GIs{$splitted_line[1]}{$splitted_line[0]}="";
  }
  if($splitted_line[1] eq ""){$unmapped{$splitted_line[0]}="";}
 }
}
my $UPs=keys %UPs;
my $GIs=keys %GIs;
my $unmapped=keys %unmapped;
my $UP_tots=keys %UP_tots;
my $GI_tots=keys %GI_tots;
print "Got $UP_tots and $GI_tots\n$UPs UPs were mapped onto $GIs GIs\n$unmapped GIs remained unmapped\n";

print "Genuine multiple mapping events\n---------------------\nAcc\tgeneID\n";
foreach my $acc(sort keys %UPs)
{
 @ids=keys %{$UPs{$acc}};
 if(scalar @ids > 1)
 {
  print "$acc\t@ids\n";
 }
}
print "------------------------\ngeneID\tAcc\n";
foreach my $id(sort keys %GIs)
{
 @acc=keys %{$GIs{$id}};
 if(scalar @acc > 1)
 {
  print "$id\t@acc\n";
 }
}
print "\n************Step 1 of mapping complete***************\n\n";#exit();

# Step2: Compute the total number of disease-associated genes (geneIDs) from DisGeNet satisfying a certain score cut-off

print "Give me cut-off score (Plz refer to http://www.disgenet.org/web/DisGeNET/menu/dbinfo#score to determine appropriate score cut-off)\t";my $score=<STDIN>;chomp($score);
print "Give me min genes for a dis to be considered\t";my $genelimit=<STDIN>;chomp($genelimit);

open(fr,"/DisGeNet_Db_v5/all_gene_disease_associations.tsv")||print "Whr's DisGeNet? $!\n";
my @data=<fr>;chomp(@data);

my @splitted_line;my %dis;my $disID;my %genes;my %scores;my $counter=0;
my %for_enrichment;my %genes_for_enrichment;my $nomatch=0;

foreach my $line(@data)
{
 if(($line!~/\#/)&&($line!~/^geneId/))
 {
  @splitted_line=split("\t",$line);
  $disID=$splitted_line[2];
  if($splitted_line[4]>=$score)
  {
   $dis{$disID}{$splitted_line[0]}="";
   $genes{$splitted_line[0]}="";
  }
 }
}
my $dis=keys %dis;
my $genes=keys %genes;
print "$genes genes(GIs) associated w $dis diseases such that GDA score >=$score\n";

# Step3: Compute how many of the genes from DisGeNet selected above are prot-coding genes

open(fr2,"/DisGeNet_Db_v5/Homo_sapiens.gene_info")||print "cannot open mapping $!\n";   # This file obtained from NCBI FTP site (ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/)
my @data2=<fr2>;chomp(@data2);

my @splitted_line2;my %coders;my %coding;my %common;my $nomatch=0;my %common2;my %dis2;my @arr;my $len;

foreach my $line(@data2)
{
 @splitted_line2=split("\t",$line);
 if($splitted_line2[9] eq "protein-coding")
 {
  $coders{$splitted_line2[1]}="";
 }
}
my $coders=keys %coders;print "Hs proteome has $coders protein-coding genes\n";

foreach my $dis(sort keys %dis)                 # diseases from DisGeNet you selected based on score cut-off alone
{
 foreach my $gene(sort keys %{$dis{$dis}})      # these are the genes associated w those disease ssatisfying the GDA score-cut; but irresp whether they are protein-coding or not
 {
  foreach my $coder(sort keys %coders)          # all protein coding genes in Hs
  {
   if($gene eq $coder)
   {
    $dis2{$dis}{$coder}="";                     # %dis2 stores diseases with only prot-coding genes associated w them
   }
  }
 }
}
foreach my $dis(sort keys %dis2)
{
 @arr=keys %{$dis2{$dis}};
 $len=@arr;
 if($len>=$genelimit)                           # Only those diseases will b considered where the associated prot-coding genes>=$genelimit
 {
  foreach my $gene(sort keys %{$dis2{$dis}})
  {
   $for_enrichment{$dis}{$gene}="";
   $genes_for_enrichment{$gene}=""; 
  }
 }
}
my $diseases=keys %for_enrichment;
my $genes_for_enrichment=keys %genes_for_enrichment;
print "I will consider only $diseases diseases associated w $genes_for_enrichment protein-coding genes for enrichment analyses\n";

=headmy $m;my @m;
foreach my $dis(sort keys %for_enrichment)
{
 @m=keys %{$for_enrichment{$dis}};#print "$dis\t$m\n";
 $m=@m;print "$dis\t@m\n";
}
exit();
=cut

=headforeach my $gene1(sort keys %genes_for_enrichment)
{
 foreach my $gene2(sort keys %coders)
 {
  if($gene1 eq $gene2)
  {
   $coding{$gene1}="";
  }
 }
}
my $genes2=keys %coding;
print "$genes2 of the total $genes_for_enrichment disease-associated genes are part of $coders protein-coding genes\n";
=cut

# Step4: Compute how many of the predicted CN-partners are listed as prot-coding genes

foreach my $id1(sort keys %GIs)           # geneIDs of predicted Cn-partners
{
 foreach my $id2(sort keys %coders)       # total "protein-coding" genes in Hs proteome
 {
  if($id1 eq $id2)
  {
   $common{$id1}="";
  }
  else
  {
   $nomatch++;
   if($nomatch eq $coders)
   {
    print "$id1\n";
   }
  }
 }
 $nomatch=0;
}
my $common=keys %common;print "$common geneIDs of predicted Cn-partners are listed as protein-coding in Hs-genome\ngeneIDs listed above are pseudogenes as per Gene record; but prots are reviewed in UP\n";

# Step5: Determine if prot-coding genes among CN-partners are enriched for prot-coding and dis-associated genes to be considered for EA in gen

foreach my $id1(sort keys %GIs)                         # geneIDs of predicted Cn-partners
{
 foreach my $id2(sort keys %genes_for_enrichment)       # only protein-coding genes to be considered for enrichment
 {
  if($id1 eq $id2)
  {
   $common2{$id1}="";
  }
  else
  {
   $nomatch++;
   if($nomatch eq $genes_for_enrichment)
   {
    #print "predicted Cn-partner $id1 is not a part of prot-coding genes to be considered for EA\n";
   }
  }
 }
 $nomatch=0;
}
my $common2=keys %common2;print "$common2 geneIDs of predicted Cn-partners are part of $genes_for_enrichment protein-coding genes from Hs-genome to be considered for enrichment analyses\n";
print "\n-----------Step5 of determining whether predicted Cn-partners are enriched in diseases complete-----------------\n";

# Step6: Rank diseases in decreasing order of enrichment p-values
# Step6.1: Prepare a table

my @m;my $m;my $x;my $n;my @arr2;my %associated_diseases;

system("rm -f GDA_enrichment_analyses");
open(fw,">>GDA_enrichment_analyses");

foreach my $dis(sort keys %for_enrichment)                  # stores each dis and associated prot-coding genes after applying the filter genes>=$genelimit
{
 @m=keys %{$for_enrichment{$dis}};#print "$dis\t$m\n";      # @m stores associated genes such that all are prot-coding and min # of genes is $genelimit
 @m=sort {$a<=>$b} @m;
 $m=@m;#print "$dis\t@m\n";
 #$n=$coders-$m;
 $n=$genes_for_enrichment-$m;                               # genes associated w all other diseases (NOT associated w disease being considered)
 foreach my $id1(sort keys %{$for_enrichment{$dis}})        # same as for every gene in @m
 {
  foreach my $id2(sort keys %GIs)
  {
   if($id1 eq $id2)
   {
    $x++;                                                   # actual # of CN-partners associated w the disease being considered
    push(@arr2,$id1);
    $associated_diseases{$dis}="";                          # And the disease being considered
   }
  }
 }
 @arr2=sort {$a<=>$b} @arr2;
 print fw "$dis\t$x\t@arr2\t";
 foreach my $gi(sort @arr2)
 {
  foreach my $acc(sort keys %{$GIs{$gi}})
  {
   print fw "$acc, ";
  }
 }
 print fw "\t$m\t@m\t$n\t$common2\n";
 $x=0;undef @arr2;#exit();
}
my $associated_dis=keys %associated_diseases;print "$associated_dis associated diseases\n";

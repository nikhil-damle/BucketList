#!/usr/bin/perl
use strict;

open(fr,"/HsCNome_analyses_Feb2018/GO_analyses/MFs/PK_act.txt")||print "cannot open HsCNome_Feb2018 $!\n";
my @pred=<fr>;chomp(@pred);

my @splitted_line;my %predkins;my %predkinsubs;my %kins;my %KS;my $nosub=0;my $nokin=0;
my %subs;my %pairs;my $pair;my $PKgene;my $subgene;my %submapping;my %kinmapping;my %kingene;my %substrate_genes;my %kinmapping2;my %submapping2;my @acc;my $canonical_sub;my $canonical_kin;
my %subsELM;
my %preds;my $i;my $kinacc;my $subacc;my $kinflag=0;my $subflag=0;my %pairs_in_hscnome;my %kins_in_hscnome;my %subs_in_hscnome;
my $nomatch=0;my $nomatch2=0;
my $int;my %PNet;my $p1;my $p2;
my $len;my $subcount=0;my %unmapped_kins;my %kinonly;my %CNK;my %CNS;
my $match=0;my %KS_others;my %other_targets;my $nomatch2=0;my $nomatch3=0;my %kinrest;my %subrest;my %KS_rest;my %CNS2;my %remaining_pSubs;

foreach my $line(@pred)
{
 @splitted_line=split("\t",$line);
 $predkins{$splitted_line[1]}="";
}
my $predkins=keys %predkins;print "HsCNome has $predkins proteins annotated with PK_act as MF GO term\n";

open(fr,"/PO4_Site_plus/Kinase_Substrate_Dataset")||print "cannot open Kinase_Substrate_Dataset $!\n";
my @data=<fr>;chomp(@data);

foreach my $line(@data)
{
 @splitted_line=split("\t",$line);
 if(($splitted_line[3] eq "human")&&($splitted_line[8] eq "human"))
 {
  if($splitted_line[6]=~/\-/){$canonical_sub=$`;}#print "$canonical_sub\n";exit();}
  else{$canonical_sub=$splitted_line[6];}
  if($splitted_line[2]=~/\-/){$canonical_kin=$`;}
  else{$canonical_kin=$splitted_line[2];}

  $kins{$canonical_kin}{$canonical_sub}="";
  $subs{$canonical_sub}{$canonical_kin}="";
  $pair=$canonical_kin."_".$canonical_sub;
  $pairs{$pair}="";
  $PKgene=uc($splitted_line[0]);
  $subgene=uc($splitted_line[7]);
  $kingene{$PKgene}{$subgene}="";
  $substrate_genes{$subgene}{$PKgene}="";
  $kinmapping{$canonical_kin}=$PKgene;
  $kinmapping2{$PKgene}{$canonical_kin}="";
  $submapping{$canonical_sub}=$subgene;
  $submapping2{$subgene}{$canonical_sub}="";
 }
}
my $kins=keys %kins;my $subs=keys %subs;my $pairs=keys %pairs;my $kingene=keys %kingene;my $substrates=keys %substrate_genes;my $mapped_subs=keys %submapping;my $mapped_kins=keys %kinmapping;
print "PSP has $kins and $subs human PKs and their resp subs (target PKs kwn) forming $pairs PK-subs pairs\nThese correspond to $kingene PK genes and $substrates substrate genes\n";

print "Multiple acc mapping for PK genes\n";
foreach my $gene(sort keys %kinmapping2)
{
 @acc=keys %{$kinmapping2{$gene}};
 if(scalar @acc >1)
 {
  print "$gene\t@acc\n";
 }
}
print "Multiple acc mapping for SUBSTRATE genes\n";
foreach my $gene(sort keys %submapping2)
{
 @acc=keys %{$submapping2{$gene}};
 if(scalar @acc >1)
 {
  print "$gene\t@acc\n";
 }
}

open(fr,"HsCNome_Feb2018")||print "cannot open HsCNome_Feb2018 $!\n";   # Tab-separated file w first coloumn being the prot-acc
my @pred=<fr>;chomp(@pred);

for($i=1;$i<scalar @pred;$i++)
{
 @splitted_line=split("\t",$pred[$i]);
 $preds{$splitted_line[0]}="";
}
my $pred=keys %preds;print "You have $pred entries in HsCNome\n";
my %match;
foreach my $acc1(sort keys %preds)
{
 foreach my $acc2(sort keys %subs)    # %subs{subAcc}{kinAcc}=""; obtained from PSP above
 {
  if($acc1 eq $acc2)
  {
   $match{$acc1}="";
  }
 }
}
my $match=keys %match;print "Predicted CN-partners contain $match PProts such that upstream PKs are kwn in PSP\n";#exit();

foreach my $acc1(sort keys %predkins)  # predicted PKs in HsCNome
{
 foreach my $acc2(sort keys %kins)     # All PKs in PSP %kins{kinAcc}{subAcc}=""
 {
  if($acc1 eq $acc2)                   # if a PK in PSP is also in the HsCNome
  {
   foreach my $sub(sort keys %{$kins{$acc1}})   # Checking every substrate reported in PSP of a PK in HsCNome
   {
    #$predkinsubs{$acc1}{$sub}="";      # Only those PKs whose subs are also present in HsCNome; but this hash stores all subs of that PK in PSP

    foreach my $pred(sort keys %preds)  # %pred => all entries in HsCNome
    {
     if($pred eq $sub)                  #if the sub is also in HsCNome
     {
      $predkinsubs{$acc1}{$sub}="";        # Kin-sub hash: specifically extracts the subs present in the HSCNome of PKs also present in HsCNome
      $subcount++;
      $subs_in_hscnome{$sub}{$acc1}="";    # Sub-kin hash: specifically extracts the subs present in the HSCNome of PKs also present in HsCNome
     }
    }
   }
   if($subcount eq 0) # No subs of that PK is in the HsCNome
   {
    $kinonly{$acc1}="";
    #print "$acc1 is a CN-partner but none of its substrates are\n";
   }
   $subcount=0;
  }
  else
  {
   $nomatch++;
   if($nomatch eq $kins)
   {
    $unmapped_kins{$acc1}="";
    print "$acc1 has no kwn substrates in PSP\n";
   }
  }
 }
 $nomatch=0;
}
my $predkinsubs=keys %predkinsubs;
my $subs_in_hscnome=keys %subs_in_hscnome;
my $kinonly=keys %kinonly;
my $unmapped_kins=keys %unmapped_kins;
print "HsCNome contains $predkinsubs kinases and their $subs_in_hscnome substrates\n$kinonly PKs such that kinases are CN-partners but none of their substrates are\nAND $unmapped_kins PKs for which no human substrates are reported in PSP\n";
#foreach my $acc(sort keys %kinonly){print "$acc\n";}
foreach my $kin(sort keys %predkinsubs){print "$kin\n";}
#### Layer1: Every PK AND/OR its sub is in HsCNome #####
print "\nLayer1 starting...Fetching all the kinases in HsCNome and their $subs_in_hscnome subs that are also present in HsCNome\n";
################## Getting CN-sub pairs ##################
foreach my $acc1(sort keys %subs_in_hscnome)    # %subs_in_hscnome{subAcc}{kinAcc}=""; Every subAcc and every kinAcc present in HsCNome
{
 foreach my $acc2(sort keys %submapping)        # %submapping{subAcc}=subGene
 {
  if($acc1 eq $acc2)
  {
   if(($submapping{$acc2} ne "")&&($submapping{$acc2} ne "	"))
   {
    $int="CN_$submapping{$acc2}";
    $PNet{$int}="";
    $CNS{$int}="";
    $subcount++;
   }
   else{print "Map $acc2 manually\n";}
  }
  else
  {
   $nomatch++;
   if($nomatch eq $mapped_subs)
   {
    print "$acc1 was not mapped onto any GeneName\n";
   }
  }
 }
 $nomatch=0;
################ Getting CN-Kin pairs #################
 foreach my $kinacc(sort keys %{$subs_in_hscnome{$acc1}})
 {
  foreach my $acc2(sort keys %kinmapping)
  {
   if($kinacc eq $acc2)
   {
    if(($kinmapping{$acc2} ne "")&&($kinmapping{$acc2} ne "      "))
    {
     $int="CN_$kinmapping{$acc2}";
     $PNet{$int}="";
     $CNK{$int}="";
    }
    else{print "Map $acc2 manually\n";}
   }
   else
   {
    $nomatch++;
    if($nomatch eq $mapped_kins)
    {
     print "$kinacc was not mapped onto any GeneName\n";
    }
   }
  }
  $nomatch=0;
 }
############### Getting Kin-Subs pairs #################
 foreach my $acc2(sort keys %submapping)
 {
  if($acc1 eq $acc2)
  {
   if(($submapping{$acc2} ne "")&&($submapping{$acc2} ne "      "))
   {
    foreach my $kinacc(sort keys %{$subs_in_hscnome{$acc1}})
    {
     foreach my $kinacc2(sort keys %kinmapping)
     {
      if($kinacc eq $kinacc2)
      {
       if(($kinmapping{$kinacc2} ne "")&&($kinmapping{$kinacc2} ne "      "))
       {
        $int=$submapping{$acc2}."_".$kinmapping{$kinacc2};
        $PNet{$int}="";
        $KS{$int}="";
       }
      }
     }
    }
   }
  }
 }
}
my $PNet=keys %PNet;
my $KS=keys %KS;
my $CNS=keys %CNS;
my $CNK=keys %CNK;
print "Layer1: $subs_in_hscnome substrates of $predkinsubs kinases form $KS KS-pairs, $CNS CN-subs pairs AND $CNK CN-PK pairs\n";
######## %PNet holds PKs w predicted SLiMs and their subs w predicted SLiMs as CN-K, CN-S and KS pairs
#foreach my $pair(sort keys %PNet)
#foreach my $pair(sort keys %CNS)
#foreach my $pair(sort keys %CNK)
#foreach my $pair(sort keys %KS)
#{
# print "$pair\n";
#}
#exit();

########### Adding CN-Konly pairs - PKs have predicted SLiMs but none of their subs ######
=headforeach my $kinacc(sort keys %kinonly)
{
 foreach my $acc2(sort keys %kinmapping)
 {
  if($kinacc eq $acc2)
  {
   if(($kinmapping{$acc2} ne "")&&($kinmapping{$acc2} ne "      "))
   {
    $int="CN_$kinmapping{$acc2}";
    $PNet{$int}="";
    $CNK{$int}="";
   }
   else{print "Map $acc2 manually\n";}
  }
  else
  {
   $nomatch++;
   if($nomatch eq $mapped_kins)
   {
    print "$kinacc was not mapped onto any GeneName\n";
   }
  }
 } 
 $nomatch=0;
}
my $PNet=keys %PNet;
my $KS=keys %KS;
my $CNS=keys %CNS;
my $CNK=keys %CNK;
print "After adding kinases w predicted SLiMs such that none of their subs have predicted SLiMs\n$subs_in_hscnome substrates of $predkinsubs kinases form $KS KS-pairs, $CNS CN-subs pairs AND $CNK CN-PK pairs\n";
=cut
print "************Layer1 is thus complete*************\n";
#foreach my $kin(sort keys %predkinsubs)
#{
 #print "$kin\n";
#}
#exit();
##### Layer2: Other PKs in PSP not in HsCNome but still PO4rylate aforementioned subs #######
print "\nStarting Layer2....Fetching all other PKs from PSP not in the HsCNome but still PO4-rylate these $subs_in_hscnome subs\n";
print "NOTE: $subs_in_hscnome subs are such that some of their upstream kinases, are also present in HsCNome & now I'll fetch their remaining upstream kinases that are not present in HsCNome. I AM NOT LOOKING FOR UPSTREAM PKs OF ALL PO4-PROTS IN HsCNome\n";
foreach my $acc1(sort keys %subs_in_hscnome)    # subs in HsCNome such that the upstream PKs are also in HsCNome
{
 foreach my $acc2(sort keys %subs)              # all subs in PSP
 {
  if($acc1 eq $acc2)
  {
   foreach my $acc3(sort keys %submapping)
   {
    if($acc1 eq $acc3)
    {
     if(($submapping{$acc1} ne "")&&($submapping{$acc1} ne "      "))    # if that sub-acc has a proper geneName
     {
      foreach my $kinacc(sort keys %{$subs{$acc1}})    # all PKs PO4rylating that sub
      {
       foreach my $predkin(sort keys %predkins)        # all PKs in HsCNome PO4rylating that sub 
       {
        if($kinacc ne $predkin)
        {
         $match++;
         if($match eq $predkins)                       # the PK in question is other than the ones already present in the HsCnome
         {
          foreach my $kinacc2(sort keys %kinmapping)
          {
           if($kinacc eq $kinacc2)
           {
            if(($kinmapping{$kinacc2} ne "")&&($kinmapping{$kinacc2} ne "      "))  # if that PK has a proper geneName
            {
             $int=$submapping{$acc1}."_".$kinmapping{$kinacc2};    # KS
             $PNet{$int}="";                                       # Adding those KS pairs to PNet
             $KS_others{$int}="";
             $other_targets{$kinmapping{$kinacc2}}="";
            }
           }
           else{$nomatch++;if($nomatch eq $mapped_kins){print "Kinase $kinacc could not be mapped\n";$nomatch=0;}}
          }
         }
        }
        else{last();}#print "$kinacc is a CN-partner\n";last();}
       }
       $nomatch=0;$match=0;
      }
     }
    }
    else
    {
     $nomatch2++;
     if($nomatch2 eq $mapped_subs){print "substrate $acc1 did not get mapped\n";$nomatch2=0;}
    }
   }
   $nomatch2=0;
  }
  else{$nomatch3++;if($nomatch3 eq $subs){print "substrate $acc1 is not present in PSP\n";$nomatch3=0;last();}}
 }
 $nomatch3=0;
}
my $other_PKs=keys %other_targets;
my $KS_others=keys %KS_others;
print "These substrates are also PO4-rylated by $other_PKs other PKs that are not in the HsCNome\n";#exit();
#foreach my $pair(sort keys %PNet)
#{
# print "$pair\n";
#}
#print "*************Layer2 complete************\n";
#exit();
#### Rest of PSP PO4prots in HsCNome and their upstream PKs irrepsctive of their presence in HsCNome ####
##### NOTE: Remaining PO4prots are in HsCNome, not necessarily their upstream PKs ####

print "Layer3: Fetching all remaining PO4-prots in HsCNome and their upstream PKs\n";
$nomatch=0;$nomatch2=0;
foreach my $acc1(sort keys %preds)   # HsCNome
{
 foreach my $sub(sort keys %subs_in_hscnome) # subs in HsCNome such that the upstream PKs are also in HsCNome
 {
  if($acc1 ne $sub)
  {
   $nomatch++;
   if($nomatch eq $subs_in_hscnome)  # these may not get PO4rylated or they may be PProts w/o kwn upstream PK or the upstream PKs are not in the HsCNome
   {
    $subrest{$acc1}="";
    foreach my $acc2(sort keys %subs)  # PO4prots in PSP w kwn upstream PKs
    {
     if($acc1 eq $acc2)
     {
      $remaining_pSubs{$acc1}="";      # last case: subs in HsCNome w kwn upstream PKs but none of these upstream PKs are in HsCNome (so these subs did not get accomodated in layer1 or 2)
      foreach my $acc3(sort keys %submapping)
      {
       if($acc1 eq $acc3)
       {
        if(($submapping{$acc1} ne "")&&($submapping{$acc1} ne "      "))
        {
         $int="CN_$submapping{$acc1}";
         $PNet{$int}="";
         $CNS2{$int}="";#print "$submapping{$acc1}\n";}#exit();
        }
        else{print "Map $acc1 manually\n";}

        foreach my $kinacc(sort keys %{$subs{$acc1}})
        {
         foreach my $kinacc2(sort keys %kinmapping)
         {
          if($kinacc eq $kinacc2)
          {
           if(($kinmapping{$kinacc2} ne "")&&($kinmapping{$kinacc2} ne "      "))
           {
            $int=$submapping{$acc1}."_".$kinmapping{$kinacc2};
            $PNet{$int}="";
            $KS_rest{$int}="";
            $kinrest{$kinmapping{$kinacc2}}="";
           }
          }
          else
          {
           $nomatch2++;
           if($nomatch2 eq $mapped_kins)
           {
            print "Kinase $kinacc was not mapped\n";
           }
          }
         }
         $nomatch2=0;
        }
       }
       else
       {
        $nomatch++;
        if($nomatch eq $mapped_subs)
        {
         print "$acc1 was not mapped onto any GeneName\n";
        }
       }
      }
      $nomatch=0;
     }
    }
   }
  }
 }
 $nomatch=0;
}
my $PNet=keys %PNet;
my $KS_rest=keys %KS_rest;
my $CNS2=keys %CNS2;
my $subrest=keys %subrest;
my $kinrest=keys %kinrest;
my $remaining_pSubs=keys %remaining_pSubs;
print "$remaining_pSubs of the $subrest remaining subs in HsCNome are PO4-rylated by $kinrest kinases form $CNS2 KS-pairs\n";
#exit();

=headforeach my $acc1(sort keys %subs_in_hscnome)
{
 #if($acc1 eq "P49790"){print "Checkign $acc1 in HsCNome\t";
 foreach my $acc2(sort keys %submapping)
 {
  if($acc1 eq $acc2)
  {print "i.e. $acc2 i.e.\t";
   if(($submapping{$acc2} ne "")&&($submapping{$acc2} ne "      "))
   {
    foreach my $gene(sort keys %substrate_genes)
    {
     if($gene eq $submapping{$acc2})
     {print "$gene\n";
      foreach my $kingene(sort keys %{$substrate_genes{$gene}})
      {
       print "$gene\tpp\t$kingene\n";
       #print "$kingene\n";
       $int="$gene\_$kingene";
       $PNet{$int}="";
      }exit();
     }
     else
     {
      $nomatch2++;
      if($nomatch2 eq $substrates)
      {
       print "$acc2 i.e. $submapping{$acc2} did not have a kwn target PK\n";$nomatch2=0;
      }
     }
    }
   }
  }
  else
  {
   $nomatch++;
   if($nomatch eq $mapped_subs)
   {
    print "$acc1 was not mapped onto any GeneName\n";
   }
  }
  $nomatch2=0;
 }
 $nomatch=0;#}
}
my $PNets=keys %PNet;
=cut
open(fw,">>/HsCNome_analyses_Feb2018/all_PSP_PKs_targeting_subs_in_HsCNome");
foreach my $int(sort keys %PNet)
{
 $int=~/(.*)\_(.*)/;
 $p1=$1;$p2=$2;
 print fw "$p1\tpp\t$p2\n";
}

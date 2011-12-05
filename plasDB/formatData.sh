#!/bin/bash
#######
# 1        Gene ID	
# 2   x 1  Genomic Sequence ID	
# 3   x 2  Organism	
# 4   x 3  Chromosome	
# 5   x 4  Chromosome Order	
# 6   x 5  Genomic Location	
# 7   x 6  Gene Strand	
# 8   x 7  Gene Type	
# 9     8  # Exons	
# 10  - 9  Is Pseudo	
# 11  x 10 Transcript Length	
# 12  x 11 CDS Length	
# 13  x 12 Product Description	
# 14    13 Molecular Weight	
# 15    14 Isoelectric Point	
# 16  x 15 Protein Length	
# 17    16 # TM Domains	
# 18  x 17 SignalP Scores	
# 19  x 18 SignalP Peptide	
# 20  x 19 Pf-iRBC min expr ratio (GS array)	
# 21  x 20 Pf-iRBC min expr time (GS array)	
# 22  x 21 Pf-iRBC max expr ratio (GS array)	
# 23  x 22 Pf-iRBC max fold-induction (GS array)	
# 24  x 23 Pf-iRBC max expr %ile (GS array)	
# 25  x 24 Pf-iRBC max expr time (GS array)	
# 26  x 25 Pf-iRBC+Spz+Gam min expr level (Affy)	
# 27  x 26 Pf-iRBC+Spz+Gam min expr stage (Affy)	
# 28  x 27 Pf-iRBC+Spz+Gam max expr level (Affy)	
# 29  x 28 Pf-iRBC+Spz+Gam max expr %ile (Affy)	
# 30  x 29 Pf-iRBC+Spz+Gam max expr stage (Affy)	
# 31  o 30 EC Numbers	
# 32  x 31 Annotated GO Function	
# 33  x 32 Annotated GO Process	
# 34  x 33 Annotated GO Component	
# 35  x 34 Predicted GO Function	
# 36  x 35 Predicted GO Process	
# 37  x 36 Predicted GO Component	
# 38    37 Ortholog count	
# 39    38 Paralog count	
# 40  x 39 Ortholog Group	
# 41    40 Total SNPs All Strains	
# 42    41 Nonsynonymous SNPs All Strains	
# 43    42 Synonymous SNPs All Strains	
# 44    43 Nonsyn/Syn SNP Ratio All Strains	
# 45  x 44 Predicted Protein Sequence	
# 46  x 45 Weight	#always 10	
#.45  x 46 empty string?
# 47    47 Max Fold Induction
######

#join -a1 -j1 -t"	"  <(sort pf5_summary.txt|cut -f1,9,10,14,15,17,38,39,41-44,46)   <(sort llinas_maxOverAll.txt)| sed 's/null/0/g;s/\t/,/g;'

join -a1 -j1 -t"	"  <(sort pf5_summary.txt)   <(sort llinas_maxOverAll.txt) |
 sed 's/No/0/g;s/Yes/1/g;s/null/0/g;s/\[//g;s/\]//g;'|
 perl -Mv5.10 -F"\t"  -slane '
  BEGIN{
      @fields=qw/0 8 13 14 16 37 38 40 41 42 43/; 
      $labelIDX=47; 

      @head;
  }
  if($F[9] == 1){next};  #skip pseudo
  if($.==1){@head=@F[@fields,$labelIDX];} #set titles for first line
  else{push @data, join(",",@F[@fields],($F[$labelIDX]>9)?1:0 )} #join all lines and 0 if not diff, 1 if it is

  END{ 
    open $csv, ">labeled.csv";
    print $csv join(",",@head);
    print $csv $_ for (@data);

    open $arff, ">labeled.arff";
   #attribute line (numeric unless matches binary
    print $arff "\@relation GENES"; # . localtime;
    for (@head) { 
      $type="numeric";
      given($_) {
	when(m/Max Fold Induction/) { $type="{0,1}" }
	when(m/Gene ID/) { $type="string" }
      }
	print $arff "\@attribute \"$_\" $type" unless $_ eq "Gene ID"; #decided not to use gene ID in weka
     }
    print $arff "";
    print $arff "\@data";
    print $arff join(",",(split /,/)[1..$#fields+1]) for (@data); #fancy moves to leave out Gene ID
  } '



#get positions from summary files
perl -F"[\t]" -slane 's/,//g; m/(\d+) - (\d+) \((\+|-)\)/; $s=$1;$e=$2; ($s,$e)=($2,$1) if $3 eq "-";  print "chr",$F[3],":",$s,"-",$e;' pf5_summary.txt 
#get length and gc from that
while read pos; do  echo -ne "$pos\t"| perl -ne 'if(/Not As/){print $_,"\n"}else{print $_}'; 2bit2seq.pl $pos|tr -d "\n" | perl -lne '$gc=tr/cCGg/*/;chomp; print length($_),"\t", sprintf("%.3f",$gc/length($_))'; done <positions.txt > length_g

#get average gc
awk  '(!/Not A/ && NR!=1){sum+=$3;count+=1;}END{print sum/count}' length_gc 
#0.230274
perl -F"\t" -slane '$,="\t"; if(/Not As/){chomp;m/:(\d+)-(\d+)/; print $_,abs($2-$1),"0.23"}else{print @F}' length_gc > length_gc2
mv length_gc2 length_gc
paste <(cut -f1 pf5_summary.txt) <(cut -f2,3 length_gc) > id_length_gc



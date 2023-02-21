#!/usr/bin/env bash 

name=$1
dir=$2 
noise=$3

  seed="4512"
  noise_val=$(echo "$noise" | bc)
  sed -i '/^>/d' ${dir} 
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' ${dir} > lin_fasta.fasta
  awk -v seed=$seed -v noise_val=$noise_val 'BEGIN {
    srand(seed)
  }
  !/^>/ {
    fast_len=length($0)
    i=int(noise_val*fast_len)
    for(j=1;j<=i;j++) {
      position=int(rand()*length($0))+1
      $0=substr($0,1,(position-1)) substr($0,(position+1))
    }
  }
  1' lin_fasta.fasta >  noise_${name}_${noise}.fasta
  sed -i '1s/^/> contig_1\n/' noise_${name}_${noise}.fasta 

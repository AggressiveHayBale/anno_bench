#!/usr/bin/env bash 

name=$1
dir=$2 
noise=$3

seed="4512"
sed -i '/^>/d' ${dir} 
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' ${dir} > tmp_noise.fasta
awk -v seed=$seed -v noise_val=$noise 'BEGIN {
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
  1' tmp_noise.fasta > noise_${noise}_${name}.fasta
sed -i '1s/^/> contig_1\n/' noise_${noise}_${name}.fasta
rm tmp_noise.fasta


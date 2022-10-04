#!/usr/bin/env bash 

name=$1
dir=$2 
noise=$3

awk -v seed=$RANDOM -v noise=$noise 'BEGIN {} !/^>/ {
    srand(seed)                                        
    i=int(noise*(length($0)))                   
    for(j=1;j<=i;j++) {
        position=int(rand()*length($0))+1             
        $0=substr($0,1,(position-1)) substr($0,(position+1))  
    }
}1' $dir >  ${name}_noise.fasta
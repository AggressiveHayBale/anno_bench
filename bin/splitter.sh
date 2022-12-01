#!/usr/bin/env bash 

name=$1
dir=$2 
contigs=$3


cat ${dir} |  grep -v '^>' | sed -e '1i\\>contig_1' > split_${name}.fasta

line_count=$(wc -l split_${name}.fasta | cut -f1 -d' ')

contig_lines=$(shuf -i 1-${line_count} -n ${contigs})

ct=1
for line in ${contig_lines}
do
    echo "$line"
    ((ct++)) 
    sed -i "${line}i >contig_${ct}" split_${name}.fasta
done
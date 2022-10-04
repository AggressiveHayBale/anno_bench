#!/usr/bin/env bash 

name=$1
dir=$2 
noise=$3

char_count=$(wc -m ${dir} | cut -f1 -d' ')

noise_char_count=$(awk "BEGIN {print($noise/100*$char_count)}"); echo "noise_char_count=$noise_char_count"

char_to_delete=$(shuf -i ${char_count} -n ${noise_char_count})

cut -c $char_to_delete ${dir}
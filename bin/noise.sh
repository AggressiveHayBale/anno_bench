#!/usr/bin/env bash 

name=$1
dir=$2 
noise=$3

char_count=$(wc -m ${dir} | cut -f1 -d' ')


awk '{${char_count}*(${noise}/100)}'
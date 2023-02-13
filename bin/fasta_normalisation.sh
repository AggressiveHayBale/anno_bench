#!/usr/bin/env bash 

dir=$1

sed '/^>/!s/.\{61\}/&\
/g' ${dir} 

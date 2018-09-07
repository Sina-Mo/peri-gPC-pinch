#!/bin/bash

# What should you space the jobs by?
countby=8
endrun=680
howmany=$((endrun/countby))

for i in `seq 1 $howmany`;
do

second=$((countby*i))
first=$((second-(countby-1)))
name=$first-$second
echo Print $first to $second runs. 

awk -v var="$first" 'NR==8 {$0="startrun="'"var"'""} 1' runperipinch.job > temp${i}1_1
awk -v var="$second" 'NR==9 {$0="endrun="'"var"'""} 1' temp${i}1_1 > temp${i}1_2
awk -v var="$name" 'NR==2 {$0="#SBATCH --job-name=r"'"var"'""} 1' temp${i}1_2 > go${first}to${second}.job

rm temp${i}1_1 temp${i}1_2

#sbatch go${first}to${second}.job

done


#!/bin/bash

# What should you space the jobs by?
countby=8
startrun=544
endrun=680
howmany=$(((endrun-startrun+1)/countby))

for i in `seq 1 $howmany`;
do

second=$((countby*i+startrun))
first=$((second-(countby-1)))
echo Submit $first to $second runs. 

sbatch go${first}to${second}.job

done


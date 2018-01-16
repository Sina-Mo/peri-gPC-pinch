#!/bin/bash

# Separate main parameter file into three files
# cut -f 1 allpara.txt > Wo.txt
cut -f 2 allpara.txt > pamp.txt
cut -f 3 allpara.txt > Freq.txt

# Count number of lines in files
numlines=$(grep -c "^" allpara.txt)

# initialize variables
pam=0
Freq=0
# For loop that will write files
for i in `seq 1 $numlines`;
do
# Sets Wo based on i
pamp=$(awk -v var="$i" 'NR==var' pamp.txt)
# Writes file to replace line with Wo
cat template-parameters |  awk -v var="$pamp" 'NR==6 {$0="pamp = "'"var"'"  \\ "} 1' template-parameters > parametersw${i}
# Sets Freq based on i
Freq=$(awk -v var="$i" 'NR==var' Freq.txt)
# Writes file to replace line with Freq
cat template-parameters |  awk -v var="$Freq" 'NR==7 {$0="freq = "'"var"'"   \\"} 1' parametersw${i} > parametersw${i}f${i}
cat template-parameters |  awk -v var="$Freq" 'NR==8 {$0="iposn = "'"var"'"   \\"} 1' parametersw${i}f${i} > parametersw${i}f${i}l${i}
# Cleans up folder
rm parametersw${i} parametersw${i}f${i}
mv parametersw${i}f${i}l${i} parameters${i}
done
# above works

# Some other tries down here. 
# cat input2d |  awk -v var="$variable" 'NR==7 {$0="WO = var   // Womersley number"} 1' input2d > input2d2  # This will replace a single line in input2d and copy it to input2d

# cat input2d |  awk 'BEGIN {} { for(i=1; i <=5; i++)  NR==7 $0="WO = $i   // Womersley number" 1}' input2d > input2d2

# awk 'BEGIN { print "START" } { for(i=1; i <=NR; i++) print i} END { print NR }' input2d


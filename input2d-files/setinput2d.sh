#!/bin/bash

# Separate main parameter file into three files
cut -f 1 allpara.txt > Wo.txt
# cut -f 2 allpara.txt > pamp.txt
cut -f 3 allpara.txt > Freq.txt

# Count number of lines in files
numlines=$(grep -c "^" allpara.txt)

# initialize variables
WO=0
Freq=0
# For loop that will write files
for i in `seq 1 $numlines`;
do
# Sets Wo based on i
WO=$(awk -v var="$i" 'NR==var' Wo.txt)
# Writes file to replace line with Wo
awk -v var="$WO" 'NR==7 {$0="WO = "'"var"'"   // Womersley number"} 1' template-input2d > input2dw${i}
# Sets Freq based on i
Freq=$(awk -v var="$i" 'NR==var' Freq.txt)
# Writes file to replace line with Freq
awk -v var="$Freq" 'NR==5 {$0="FRE = "'"var"'"   // Pumping Frequency (Hz)"} 1' input2dw${i} > input2dw${i}f${i}
# Edits input2d to create different IBlog and visit files and folders
awk -v var="$i" 'NR==201 {$0="   log_file_name = \"runs/IB2d.log"'"var"'"\"                //"} 1' input2dw${i}f${i} > input2dw${i}f${i}l${i}
awk -v var="$i" 'NR==207 {$0="   viz_dump_dirname = \"runs/viz_IB2d"'"var"'"\"                //"} 1' input2dw${i}f${i}l${i} > input2dw${i}f${i}l${i}y${i}
test=$Freq
# echo $test
if [ $(echo " $test <= 0.70" | bc) -eq 1 ]
then 
    cat template-input2d |  awk 'NR==24 {$0="END_TIME =      3.5                // Final Simulation Time (s)"} 1' input2dw${i}f${i}l${i}y${i} > input2dw${i}f${i}l${i}y${i}t${i}
elif [ $(echo " $test < 1.0" | bc) -eq 1 ]
then
    cat template-input2d |  awk 'NR==24 {$0="END_TIME =       3.0               // Final Simluation Time (s)"} 1' input2dw${i}f${i}l${i}y${i} > input2dw${i}f${i}l${i}y${i}t${i}
else
    cat template-input2d |  awk 'NR==24 {$0="END_TIME =       2.5                // Final Simulation Time (s)"} 1' input2dw${i}f${i}l${i}y${i} > input2dw${i}f${i}l${i}y${i}t${i}
fi

# Cleans up folder
rm input2dw${i} input2dw${i}f${i} input2dw${i}f${i}l${i} input2dw${i}f${i}l${i}y${i}
mv input2dw${i}f${i}l${i}y${i}t${i} input2d${i}

echo $i
grep "END_TIME =" input2d${i}  && grep "FRE ="  input2d${i}

done
# above works

# Some other tries down here. 
# cat input2d |  awk -v var="$variable" 'NR==7 {$0="WO = var   // Womersley number"} 1' input2d > input2d2  # This will replace a single line in input2d and copy it to input2d

# cat input2d |  awk 'BEGIN {} { for(i=1; i <=5; i++)  NR==7 $0="WO = $i   // Womersley number" 1}' input2d > input2d2

# awk 'BEGIN { print "START" } { for(i=1; i <=NR; i++) print i} END { print NR }' input2d


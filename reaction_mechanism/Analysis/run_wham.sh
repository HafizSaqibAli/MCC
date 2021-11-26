#!/bin/bash -f

if [ -z $../../../whan ]; then
   echo "Please set ../../../wham and try again. Exiting..."
   exit
fi
# 1. Make a 'meta file' referencing the corrected .tra files - note that force-constant needs to be doubled, so set to 200 (due to different definition)
for i in `seq -1.50 0.10 1.50`; do echo "rc$i/rc${i}_1.tra $i 300"; done > meta.dat
# 2. Run wham with extremes and number of bins such that we'll get bins of 0.05 Ang width and midpoints of -1.5, 1.5.
../wham/wham/wham P -1.50 1.50 41 0.0001 300 0 meta.dat wham.txt 0 1 &> wham.log

cat wham.txt | awk '{print$1,$2}' > pmf.txt

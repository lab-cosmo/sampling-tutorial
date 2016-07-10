#!/bin/sh
#
# This script runs in sequence some simulations in order to find the right mcstep

for mcstep in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15 
do
echo "Running mcstep=$mcstep "
echo "&inp
seed = 1357
temp = 0.22
dataxyz = 'lj.xyz'
nstep = 5000000
stridetrj = 50000
stridelog = 10
mcstep = $mcstep
outputf = 'out$mcstep'
&end" > input$mcstep

../../source/mcnvt input$mcstep > log$mcstep
grep -v "#" log$mcstep | awk '{if(NR!=1) print $2}' | autocorr -maxlag 5000 > acf$mcstep

done

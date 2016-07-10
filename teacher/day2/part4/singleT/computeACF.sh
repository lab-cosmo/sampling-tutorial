#!/bin/bash
#
# This script runs in sequence some simulations in order to find the right mcstep

#set -x
drop=5000
temp=0.19
for friction in 0.1 1 10 100 1000
do
 echo " >>> Friction $friction."
 echo "&inp
seed = 87123 
temp = $temp
dataxyz = '../lj$friction-$temp.xyz'
dt = 0.01
nstep = 20000000
langevinWNtau = $friction
mstep = 10
stridetrj = 50000
stridelog = 10
outputf = '2lj$friction-$temp.xyz'
&end" > input$friction-$temp
# directly compute the acf
../../../../../source/mdcode input$friction-$temp | grep -v "#" - | awk '{if(NR>2000){print $2}}' | autocorr -maxlag 50000 -timestep 10 > acf$friction-$temp &
done

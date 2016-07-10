#!/bin/bash
#
# This script runs in sequence some simulations in order to find the right mcstep

#set -x
drop=5000
for temp in 0.17 0.175 0.18 0.185 0.19 0.195 0.20 0.205 0.21 0.215 0.22
do
# perform a first equilibration
echo "Temp $temp :"
mkdir t$temp
cd t$temp
echo "&inp
seed = 1357
temp = $temp
dataxyz = '../lj.xyz'
dt = 0.01 
nstep = 1000000
langevinWNtau = 10 
mstep = 5
stridetrj = 50000
stridelog = 500
outputf = 'ljin.xyz'
&end" > inputeq

../../../../source/mdcode inputeq > logeq

# now I have a starting equilibrated configuration
# let's look for the right friction at this temp
for friction in 0.1 1 10 100  
do
 echo " >>> Friction $friction."
 echo "&inp
seed = 1357
temp = $temp
dataxyz = 'ljin.xyz'
dt = 0.01
nstep = 20000000
langevinWNtau = $friction
mstep = 10
stridetrj = 50000
stridelog = 10
outputf = 'lj$friction-$temp.xyz'
&end" > input$friction-$temp
 # directly compute the acf
 ../../../../source/mdcode input$friction-$temp | grep -v "#" - | awk '{if(NR>2500){print $4}}' | autocorr -maxlag 50000 -timestep 10 > acf$friction-$temp &
 done
# all the previous are running in background...
 friction=1000
 echo " >>> Friction $friction."
 echo "&inp
seed = 1357
temp = $temp
dataxyz = 'ljin.xyz'
dt = 0.01
nstep = 20000000
langevinWNtau = $friction
mstep = 10
stridetrj = 50000
stridelog = 10
outputf = 'lj$friction-$temp.xyz'
&end" > input$friction-$temp
 # directly compute the acf
 date
 ../../../../source/mdcode input$friction-$temp | grep -v "#" - | awk '{if(NR>2500){print $4}}' | autocorr -maxlag 50000 -timestep 10 > acf$friction-$temp
 date
 # extract the a.c times
 for friction in 0.1 1 10 100 1000
 do
   actime=$(head acf$friction-$temp | grep a.c | awk '{print $4/1}')
   echo "$friction $actime" >> acftot-$temp
 done
cd ..
done

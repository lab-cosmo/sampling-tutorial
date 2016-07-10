# !/bin/bash
# call as mkrwfes.sh T1 T2 < cvout > fes
# reweights energy data from a simulation at T=T1 to T=T2, and computes the cumulative average printing to stdout.


kt1=0.17
if [ $# -gt 0 ]; then kt1=$1; fi
kt2=$kt1
if [ $# -gt 1 ]; then kt2=$2; fi
basev=-160

awk -v kt1=$kt1 -v kt2=$kt2 -v v0=$basev '!/#/{w=exp(-($4-v0)*(1.0/kt2-1.0/kt1) + $5/kt1); t+=$4*w; tw+=w; print $1, t/tw;}' 

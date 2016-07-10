# !/bin/bash
# call as mkrwfes.sh T1 T2 < cvout > fes
# reweights CV data from a simulation at T=T1 to T=T2, and computes the 1D FES which is then printed to stdout.

kt1=0.17
if [ $# -gt 0 ]; then kt1=$1; fi
kt2=$kt1
if [ $# -gt 1 ]; then kt2=$2; fi
basev=-160

awk -v kt1=$kt1 -v kt2=$kt2 -v v0=$basev '!/#/{print $2, exp(-($6-v0)*(1.0/kt2-1.0/kt1)+$7/kt1)}' | histogram -w -xi 0 -xf 25 -whard -n 150 -t 0.1 \
    | awk -v kt=$kt2  'BEGIN{print "# n6   F(n6)" } !/#/{ if (NF==0) print ""; else printf "%15.7e %15.7e\n", $1, -kt*log($2) } ' 

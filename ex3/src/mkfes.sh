# !/bin/bash

kt=0.168
if [ $# -gt 0 ]; then kt=$1; fi

awk '!/#/{print $2}' out.cv | histogram -xi 0 -xf 25 -whard -n 150 -t 0.1 \
    | awk -v kt=$kt 'BEGIN{print "# n6   F(n6)" } !/#/{ if (NF==0) print ""; else printf "%15.7e %15.7e\n", $1, -kt*log($2) } ' >  n6.fes
awk '!/#/{print $4}' out.cv | histogram -xi 0 -xf 15 -whard -n 250 -t 0.1 \
    | awk -v kt=$kt 'BEGIN{print "# n8   F(n8)" } !/#/{ if (NF==0) print ""; else printf "%15.7e %15.7e\n", $1, -kt*log($2) } ' >  n8.fes
awk '!/#/{print $2,$3}' out.cv | histogram -xi 0 -xf 25 -whard -n 250 -t 0.1 -avg > n6.phi
awk '!/#/{print $4,$5}' out.cv | histogram -xi 0 -xf 15 -whard -n 150 -t 0.1 -avg > n8.phi

awk '!/#/{print $2,$4}' out.cv | ndhistogram -g -d 2 -xi 0,0 -xf 25,15 -whard -n 250,150 -t 0.1,0.1 -adaptive 0.1 \
    | awk -v kt=$kt 'BEGIN{print "# n6  n8  F(n6,n8)" } !/#/{ if (NF==0) print ""; else printf "%15.7e %15.7e %15.7e\n", $1, $2, -kt*log($3) } ' > n6n8.fes

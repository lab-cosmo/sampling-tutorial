#!/bin/bash
temp=0.19
rm acftot-$temp
for friction in 0.1 1 10 100 1000
do
  actime=$(head acf$friction-$temp | grep a.c | awk '{print $4/1}')
  echo "$friction $actime" >> acftot-$temp
done

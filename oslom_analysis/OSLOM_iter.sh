#!/bin/bash

#run_1 seed --> 42
#run_2 seed --> 93
#run_3 seed --> 212
#run_4 seed --> 5
#run_5 seed --> 1

net=plasmid.net

for p in $(seq 0.01 0.01 1)
do
  echo $p
  awk -v p="$p" -F " " '$3>=p{print $0}' $net >> test_0.3.net
  /home/macman/bin/OSLOM2/oslom_undir -r 250 -hr 0 -t 0.05 -singlet -seed 42 -cp 0 -f test_0.3.net -w &>> ${p}.oslo.log
  mv ${p}.oslo.log ${p}.net_oslo_files
done


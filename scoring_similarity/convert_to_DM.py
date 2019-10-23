#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 18:02:39 2019

@author: macman

1st argument list of names matching the mash, bindash output
2nd argument mash, bindash output
3rd argument name of the output csv file
"""

import sys
import numpy as np

with open(sys.argv[1]) as infile:
   lines = infile.readlines()
lines = [x.rstrip() for x in lines]

names =dict()
for i, item in enumerate(lines):
    if item not in names:
        names[item] = [i]
    else:
        names[item].append(i)

dm = np.empty((len(names),len(names)))
dm[:] = np.nan

counter=0
old_p1=""
old_p2=""
#x=[]

with open(sys.argv[2]) as infile:
    lines = infile.readlines()
    line_nr=len(lines)
    for line in lines: 
        toks = line.strip().split()
        p1 = toks[0]
        p2 = toks[1]
        a, b = tuple(toks[4].split('/'))
        jaccard = float(a) / float(b)
        if(old_p1 != p1):
            index_p1 = names[p1]
        if(old_p2 != p2):
            index_p2 = names[p2]
        
        # filling in the matrix
        dm[index_p1,index_p2] = jaccard
        if(p1 != p2):
            dm[index_p2,index_p1] = jaccard
        
        # recording old indicies
        old_p1 = p1
        old_p2 = p2
        
        # progress bar
        counter=counter+1
        sys.stdout.write("\r%i/%i lines done" % (counter,line_nr))

sys.stdout.write("\n")
sys.stdout.write("Saving table\n")
np.savetxt(sys.argv[3], dm, delimiter=",", fmt='%1.8f')
sys.stdout.write("Done!\n")

    




#!/bin/env python

import sys, os
import pylab
from pylab import *
import matplotlib.pyplot as plt

infile = open('response_time.txt', 'r')
all_lines = infile.readlines()
big_list = []
condition_A = []
condition_B = []
for lines in all_lines:
     lines = lines.strip()
     big_list.append(lines)


for ii in range(len(big_list)):
    if ii%2 == 0:
       condition_B.append(float(big_list[ii]))
    else:
       condition_A.append(float(big_list[ii]))

mean_A = sum(condition_A)/float(len(condition_A))
mean_B = sum(condition_B)/float(len(condition_B))
print "%f \n" %mean_A 
print "%f \n" %mean_B       

boxplot([condition_A, condition_B])
pylab.show()

infile.close()


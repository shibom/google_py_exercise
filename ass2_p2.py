#!/bin/env python

import os, sys

infile = open('library_data.txt', 'r')
all_lines = infile.readlines()
MAX = -10000
for lines in all_lines:
    lines = lines.split()
    student_id = lines[0]
    names = lines[1:len(lines)-1]
    names = " ".join(names)
    owe = float(lines[len(lines)-1])
  #  print "%s %s %f" %(student_id, names, owe)
    if owe > 0.0:
       print student_id
    if owe > MAX:
       MAX = owe 
print MAX, student_id    
 
infile.close()

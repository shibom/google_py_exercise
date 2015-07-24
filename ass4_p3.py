#!/bin/env python

import pylab as p
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt


infile = open('experiment.txt', 'r')
all_lines = infile.readlines()
data = []
for line in all_lines:
    line = line.split( )
    data.append([float(l) for l in line])
data = np.array(data)

height = data[:,0]; weight = data[:,1]; amt_chocolate = data[:,2]; task1 = data[:,3]; task2 = data[:,4]

p.figure(1)
plt.scatter(height, weight)
plt.xlabel('height')
plt.ylabel('weight')

#fitting the scatter plot with a line 

coeff = np.polyfit(height, weight, 1)
poly = np.poly1d(coeff)
plt.plot(height, poly(height), 'r-')

#pearson r calculation between task1 and task2

p_r =  pearsonr(task1, task2)

p.figure(2)

task = [task1, task2]
p.boxplot(task)
p.title('pearson_coeff = %.5f, %5f' %(p_r))
p.show()

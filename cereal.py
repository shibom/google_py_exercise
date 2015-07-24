#!/bin/env python

import numpy as np
import pylab 
from pylab import *
# load the cereal data set

names = []
data = []
col_labels = ["calories (number)", "protein (g)", "fat (g)", "sodium (mg)", "dietary fiber (g)", "complex carbohydrates (g)", "sugars (g)", "display shelf (1, 2, or 3, counting from the floor)", "potassium (mg)", "vitamins and minerals (0, 25, or 100, respectively)", "weight (in ounces) of one serving (serving size)", "cups per serving"]

f = open("cereal.txt")
for l in f:
  fields = l.split()
  names.append(fields[0])
  data.append([float(x) for x in fields[3:]])

data = np.array(data)

# We now have extracted the following:
#    names      List containing the names of the cereals tested 
#    col_labels Labels of the attributes in the data set
#    data       numpy.ndarray containing the numeric values 
#               of the 12 attributes for the 77 cereals

# example plot: fat content vs. calories
fat = data[:,2]
calories = data[:,0]
complex_carbs = data[:,5]
sugars = data[:,6]

# Below is the answer for Q a)
print mean(calories), mean(complex_carbs), mean(sugars)

# Below is the answer for Q b)
ind = np.argmax(max(sugars))
print "cereal with highest sugar: %s" %names[ind]

# Below is the answer for Q c) and d)
shelf_1 = []; shelf_2 = []; shelf_3 = []
sugar_1 = []; sugar_2 = []; sugar_3 = []

for ii in range(77):
    if (data[ii,7] == 1):
       shelf_1.append(data[ii,0])
       sugar_1.append(data[ii,6])
    elif (data[ii,7] == 2):
       shelf_2.append(data[ii,0])
       sugar_2.append(data[ii,6])
    else:
       shelf_3.append(data[ii,0])
       sugar_3.append(data[ii,6])
print np.mean(shelf_1), np.mean(shelf_2), np.mean(shelf_3)

figure(1)
plot(fat, calories ,'kx')
xlabel(col_labels[2])
ylabel(col_labels[0])
grid()

# fit a line through the points
coefficients = np.polyfit(fat, calories, 1)
poly = np.poly1d(coefficients)
pylab.plot(fat,poly(fat), 'r-')

sug_mat = [sugar_1, sugar_2, sugar_3]
figure(2)
boxplot(sug_mat)

pylab.show()

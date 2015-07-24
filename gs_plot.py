#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

z = np.loadtxt('all.txt')
trans = z.T
print trans.shape
print trans[1:,0].max()
print trans[1:,1].min()

plt.plot(np.log(trans[1:210,0]))
plt.plot(np.log(trans[1:210,90]))
plt.xlim(0,210)
plt.show()

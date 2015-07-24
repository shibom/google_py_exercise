#!/bin/env python

import pylab as P
from pylab import *
import matplotlib.pyplot as plt
import numpy as np

x = arange(-3,3,0.1)
y = 2*x**2 + 3
y1 = 5*x**2 + 3
fig = figure()
ax = fig.add_subplot(111)
ax.plot(x, y, 'k--', lw=1.5, label='f(x) = 2x$^2$ + 3')
ax.plot(x, y1, 'r-', lw=1.5, label='f(x) = 5x$^2$ + 3')
ax.legend()
ax.set_title('Some parabolas')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_ylim(0, 50)
ax.grid(True)

P.show()


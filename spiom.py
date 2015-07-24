#!/bin/env python

from pylab import subplot, scatter, hist, grid
import numpy as np

def splom(data, targets=None):
  (N, D) = data.shape
  if targets == None:
    targets = np.zeros(N)
  for i in range(D):
    for j in range(D):
      subplot(D, D, i*D + j + 1)
      if i == j:
        hist(data[:,i], bins=20)
      else:
        scatter(data[:,j], data[:,i], c=targets)
        grid()



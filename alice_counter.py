#!/bin/env python

from collections import Counter
import pylab
from pylab import *
import matplotlib.pyplot as plt

infile = open("alice.txt", 'r')
all_lines = infile.readlines()
num_words = 0
word_count = Counter()
words = {}
for lines in all_lines:
    lines = lines.split()
    for word in lines:
       
        word_count[word.lower()] += 1
for word, count in word_count.iteritems():
    words.update({(word, count)})
#words = sorted(words, key=words.get)
for_hist = words.values()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(for_hist, 1000)
ax.set_xlim(0,400)
plt.show()
pylab.show()

infile.close()

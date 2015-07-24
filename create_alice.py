#!/bin/env python

import urllib2

URL = "http://www.gutenberg.org/files/11/11.txt"
f = urllib2.urlopen(URL)
of = open("alice.txt", "w")
of.write(f.read())
f.close()
of.close()



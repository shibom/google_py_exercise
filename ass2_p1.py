#!/bin/env python

num = 1
while num < 101:
      if (num%3 == 0) and (num%5 == 0):
           print "FizzBuzz \n"
      elif num%3 == 0:
         print "Fizz \n"
      elif num%5 == 0:
         print "Buzz \n"
      else:
         print "%d \n" %num      
      num += 1
      

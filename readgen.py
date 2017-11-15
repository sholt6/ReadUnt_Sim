#!/usr/bin/python3
# This is a helper script for run.py and should be used to determine an
# appropriate value for that script's -c flag, in order to give the desired
# mean read length. The mean for your given value will be printed to STDOUT.
# Some common values:
#   -c Value    Approximate Mean Value
#   400         1000
#   3000        7500 (Default in run.py)
#   40000       100000
#   400000      1000000

import sys
import numpy as np
import scipy.special as sps
import matplotlib.pylab as plt

shape, scale = 2.5, 3000.

if len(sys.argv) == 1:
    try:
        scale = float(sys.argv[1])
    except:
        print(sys.argv[1] + " was not found or is not a number")
else:
    print("\nusage: python3 read_gen.py [integer]")
    print("Proceeding with default value of 3000")

s = np.random.gamma(shape, scale, 30000)
total = 0
count = 0
for i in s:
    total += i
    count += 1

mean = total/count

print("\nMean is " + str(mean))

count, bins, ignored = plt.hist(s, 50, normed=True)
y = bins**(shape-1)*(np.exp(-bins/scale) / (sps.gamma(shape)*scale**shape))
plt.plot(bins, y, linewidth=2, color='r')
plt.show()

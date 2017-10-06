#!/usr/bin/python3
# Produces a 3D scatter plot from a specified .tsv file

import sys
import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy

try:
    fnames = sys.argv[1:]
except:
    print("\nUsage: 3d_plt.py <input1.tsv> <input2.tsv> ...\n")

# Getting the data in:
headPat = '(\w+)\t(Simple)\t(Read\.Until)'
dataPat = '(\d+)\t(\d+)\t(\d+)'

heads = []
param = []
simps = []
rUnts = []
inter = []

for i in range(0, len(fnames)):
    with open(fnames[i], "r") as tsv:
        content = tsv.readlines()

    param.append([])
    simps.append([])
    rUnts.append([])

    for line in content:
        hmatch = re.match(headPat, line)
        if hmatch:
            heads.append(hmatch.group(1))
            heads.append(hmatch.group(2))
            heads.append(hmatch.group(3))

        dmatch = re.match(dataPat, line)
        if dmatch:
            param[i].append(int(dmatch.group(1)))
            simps[i].append(int(dmatch.group(2)))
            rUnts[i].append(int(dmatch.group(3)))

# Graphing the data
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xpon = rUnts
ypon = simps
zpon = param

for i in range(0, len(fnames)):
    ax.scatter(xpon[i], ypon[i], zpon[i],
               color=numpy.random.rand(3,),
               label=fnames[i])
    ax.legend()

ax.set_xlabel(heads[2])
ax.set_ylabel(heads[1])
ax.set_zlabel(heads[0])

plt.show()

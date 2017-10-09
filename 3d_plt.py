#!/usr/bin/python3
# Produces a 3D scatter plot from a specified .tsv file

import sys
import re
import matplotlib.pyplot as plt
import numpy
from mpl_toolkits.mplot3d import Axes3D
from optparse import OptionParser

usage = ("Usage: %prog [flag] [flag] [flag] <input1> <input2> <input3> ...\n"
         "Specify three parameters to graph and .tsv files to provide data")

parser = OptionParser(usage=usage)

parser.add_option("-s", dest="speed", action="store_true",
                  help="Use speed of sequencing as an axis")

parser.add_option("-i", dest="interval", action="store_true",
                  help="Use interval period as an axis")

parser.add_option("-r", dest="rejPen", action="store_true",
                  help="Use rejection penalty as an axis")

parser.add_option("-l", dest="idLag", action="store_true",
                  help="Use no. bases needed to map strand as an axis")

parser.add_option("-n", dest="simple", action="store_true",
                  help="Use duration of non-read until run as an axis")

parser.add_option("-u", dest="rUntil", action="store_true",
                  help="Use duration of read until run as an axis")

parser.add_option("-f", dest="fc", action="store_true",
                  help="Use fold change as an axis")

(options, args) = parser.parse_args()

try:
    fnames = args
    assert len(fnames) > 0
except:
    print("\nNo files specified\n")
    parser.print_help()
    quit()

# Getting the data in:
headPat = ('(Speed)\t(Interval)\t(RejPen)\t(IdLag)\t'
           '(Simple)\t(Read.Until)\t(Fold.Change)')
dataPat = '(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+\.\d+)'

heads = []
addHead = 1

speed = []
interval = []
rejPen = []
idLag = []
simple = []
readUnt = []
foldChange = []

for i in range(0, len(fnames)):
    with open(fnames[i], "r") as tsv:
        content = tsv.readlines()

    speed.append([])
    interval.append([])
    rejPen.append([])
    idLag.append([])
    simple.append([])
    readUnt.append([])
    foldChange.append([])

    for line in content:

        hmatch = re.match(headPat, line)
        if hmatch and addHead:
            for j in range(1, len(hmatch.groups())+1):
                heads.append(hmatch.group(j))
            addHead = 0

        dmatch = re.match(dataPat, line)
        if dmatch:
            speed[i].append(int(dmatch.group(1)))
            interval[i].append(int(dmatch.group(2)))
            rejPen[i].append(int(dmatch.group(3)))
            idLag[i].append(int(dmatch.group(4)))
            simple[i].append(int(dmatch.group(5)))
            readUnt[i].append(int(dmatch.group(6)))
            foldChange[i].append(float(dmatch.group(7)))


# Graphing the data
axisLabs = {"speed": "Speed (b/s)",
            "interval": "Interval between strands (s)",
            "rejPen": "Rejection Penalty (s)",
            "idLag": "Bases needed to map strand (b)",
            "simple": "Duration of non-read until (h)",
            "rUntil": "Duration of read until (h)",
            "fc": "Fold change (RU/NRU)"}
#############THIS DOESN'T WORK FIX IT NERD
for opt, value in options.__dict__.items():
    if options[opt] == value:
        print("Yes")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xpon = speed
ypon = foldChange
zpon = rejPen

for i in range(0, len(fnames)):
    ax.scatter(xpon[i], ypon[i], zpon[i],
               color=numpy.random.rand(3,),
               label=fnames[i])
    ax.legend()

ax.set_xlabel(heads[2])
ax.set_ylabel(heads[1])
ax.set_zlabel(heads[0])

plt.show()

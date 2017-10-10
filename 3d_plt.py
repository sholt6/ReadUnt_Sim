#!/usr/bin/python3
# Produces a 3D scatter plot from a specified .tsv file

import re
import matplotlib.pyplot as plt
import numpy
from mpl_toolkits.mplot3d import Axes3D
from optparse import OptionParser

# Usage and optparser
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

parser.add_option("-u", dest="readUnt", action="store_true",
                  help="Use duration of read until run as an axis")

parser.add_option("-f", dest="foldChange", action="store_true",
                  help="Use fold change as an axis")

(options, args) = parser.parse_args()

# Graphing Functions
def TwoDPlot():
    fig = plt.figure

    xpon = axes[0][0]
    ypon = axes[1][0]

    for i in range(0, len(fnames)):
        plt.scatter(xpon[i], ypon[i],
                    label=fnames[i])
        plt.legend()

    plt.xlabel(axes[0][1])
#    plt.yticks([0, 1, 2])
    plt.ylabel(axes[1][1])

    plt.show()

def ThreeDPlot():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    xpon = axes[0][0]
    ypon = axes[1][0]
    zpon = axes[2][0]

    for i in range(0, len(fnames)):
        ax.scatter(xpon[i], ypon[i], zpon[i],
                   color=numpy.random.rand(3,),
                   label=fnames[i])
        ax.legend()

    ax.set_xlabel(axes[0][1])
    plt.yticks([0, 1, 2])
    ax.set_ylabel(axes[1][1])
    ax.set_zlabel(axes[2][1])

    plt.show()


# Checking files have been specified
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
dataPat = '(\d+)\t(\d+)\t(\d+\.*\d*)\t(\d+)\t(\d+)\t(\d+)\t(\d+\.\d+)'

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
            rejPen[i].append(float(dmatch.group(3)))
            idLag[i].append(int(dmatch.group(4)))
            simple[i].append(int(dmatch.group(5)))
            readUnt[i].append(int(dmatch.group(6)))
            foldChange[i].append(float(dmatch.group(7)))

# Determining which axes to graph and ensuring 3 are specified
axes = []

if options.speed:
    speedAx = [speed, "Speed (b/s)"]
    axes.append(speedAx)
if options.interval:
    interAx = [interval, "Interval between strands (s)"]
    axes.append(interAx)
if options.rejPen:
    rejPeAx = [rejPen, "Rejection Penalty (s)"]
    axes.append(rejPeAx)
if options.idLag:
    idLagAx = [idLag, "Bases needed to map strand (b)"]
    axes.append(idLagAx)
if options.simple:
    simplAx = [simple, "Duration of non-read until (h)"]
    axes.append(simplAx)
if options.readUnt:
    rUntAx = [readUnt, "Duration of read until (h)"]
    axes.append(rUntAx)
if options.foldChange:
    fcAx = [foldChange, "Fold change (NRU/RU)"]
    axes.append(fcAx)

# Plotting the data
if len(axes) == 2:
    TwoDPlot()
elif len(axes) == 3:
    ThreeDPlot()
else:
    print("\nPlease specify two or three variables to graph\n")
    parser.print_help()
    quit()


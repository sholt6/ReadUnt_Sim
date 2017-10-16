#!/usr/bin/python3
# Produces a 3D scatter plot from a specified .tsv file

import re
import matplotlib.pyplot as plt
import numpy as np
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

parser.add_option("-b", dest="bar", action="store_true",
                  help="Plot a bar graph")

(options, args) = parser.parse_args()


# Graphing Functions
def TwoDPlot():
    xpon = axes[0][0]
    ypon = axes[1][0]

    for i in range(0, len(fnames)):
        plt.plot(xpon[i], ypon[i], marker='o',
                    label=fnames[i])
        plt.legend()
#        plt.legend(['Rejection Penalty 0', 'Rejection Penalty 1', 'Rejection Penalty 2'])
#        plt.legend(['Rejection Penalty 0', 'Rejection Penalty 1', 'Rejection Penalty 2', 'Rejection Penalty 0, Interval 0'])
#        plt.legend(['Rejection Penalty 0', 'Rejection Penalty 1', 'Rejection Penalty 2', 'Rejection Penalty 0, 128 Events'])
#        plt.legend(['Rejection Penalty 0', 'Rejection Penalty 1', 'Rejection Penalty 2', 'Rejection Penalty 0, Interval 0', 'Rejection Penalty 0, 128 Events'])
#        plt.legend(['Rejection Penalty 0, 500 Events', 'Rejection Penalty 0, 128 Events'])

#    plt.title("")             # Add a title

    plt.xlabel(axes[0][1])
#    plt.yticks([0, 1, 2])     # Control Y ticks
    plt.ylabel(axes[1][1])

#    plt.ylim((0,44))           # Control Y axis

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
    # plt.yticks([0, 1, 2])
    ax.set_ylabel(axes[1][1])
    ax.set_zlabel(axes[2][1])

    plt.show()


def BarPlot():
    fig = plt.figure()
    ax = fig.add_subplot(111)

    xpon = axes[0][0]
    ypon = axes[1][0]

    for i in range(0, len(fnames)):
        ax.bar(xpon[i], ypon[i])

    plt.legend(['100kb read average, 128b to ID', '7.5kb read average, 128b to ID'])
    ax.set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24])
    ax.set_xticklabels([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y'])
#    ax.set_yticks(np.arange(2, 6, 1))

    ax.set_xlabel(axes[0][1])
    ax.set_ylabel(axes[1][1])

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
dataPat = '(\d+)\t(\d+\.*\d*)\t(\d+\.*\d*)\t(\d+)\t(\d+)\t(\d+)\t(\d+\.\d+)'

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
            interval[i].append(float(dmatch.group(2)))
            rejPen[i].append(float(dmatch.group(3)))
            idLag[i].append(int(dmatch.group(4)))
            simple[i].append(int(dmatch.group(5)))
            readUnt[i].append(int(dmatch.group(6)))
            foldChange[i].append(float(dmatch.group(7)))

# Determining which axes to graph and ensuring 2 or 3 are specified
axes = []

if options.speed:
    speedAx = [speed, "Chromosome Number"]#"Speed (b/s)"]
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
    fcAx = [foldChange, "Read Until Y Times Faster"]
    axes.append(fcAx)

# Plotting the data
if options.bar is True:
    BarPlot()
elif len(axes) == 2:
    TwoDPlot()
elif len(axes) == 3:
    ThreeDPlot()
else:
    print("\nPlease specify two or three variables to graph\n")
    parser.print_help()
    quit()

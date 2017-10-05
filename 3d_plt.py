#!/usr/bin/python3
import re 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fname = 'speed_rpen.tsv'
entry = '(Speed)'

#Getting the data in:
with open(fname, "r") as tsv:
    lines = tsv.readlines()


    
headPat = entry + '\t(Simple)\t(Read\.Until)\t(RejPen)'
dataPat = '(\d+)\t(\d+)\t(\d+)\t(\d+)'

heads = []
param = []
simps = []
rUnts = []
inter = []

for i in lines:
        hmatch = re.match(headPat,i)
        if hmatch:
            heads.append(hmatch.group(1))
            heads.append(hmatch.group(2))
            heads.append(hmatch.group(3))
            heads.append(hmatch.group(4))
        dmatch = re.match(dataPat,i)
        if dmatch:
            param.append(int(dmatch.group(1)))
            simps.append(int(dmatch.group(2)))
            rUnts.append(int(dmatch.group(3)))
            inter.append(int(dmatch.group(4)))


#Graphing the data

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xpon = rUnts
ypon = simps
zpon = param

ax.scatter(xpon, ypon, zpon)
ax.set_xlabel('Read until (h)')
ax.set_ylabel('Simple (h)')
ax.set_zlabel('Speed')

plt.show()


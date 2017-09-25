#!/usr/bin/python
#Version 2.5 started 2017/09/13, Sam Holt
#This script simulates MinION experiments in which multiple different 
#libraries are loaded, in different quantities, and where the genomes
#are of different sizes. The overall aim is to try to determine whether read
#until is useful for optimising this experiment. This may be a means of 
#simulating sequencing of a specific chromosome from a mixture of chromosomes

#This version modifies the way coverage is mapped to the genome, and
#consequently runs much faster, since two operations are performed per read
#where previously this scaled with the read length

import sys
import library as lib
import random
import datetime
import numpy as np
import matplotlib.pyplot as plt


###FUNCTIONS

def Select(simLibs2):
	bag = []

	for i in range(0, len(simLibs2)):
		for j in range(0, simLibs2[i].ratio):
			bag.append(i)

	chosen = random.choice(bag)

	return chosen


def SimpleRun(simLibs1):
#Simulates an experiment without read until

	sel = Select(simLibs1)

	simRunT = Pore(simLibs1[sel])

	return simRunT
	

def ReadUntil(rUnLibs1):
#Simulates an experiment with read until
	#Selection is a reference to a selected Library

	sel = Select(rUnLibs1)

	if rUnLibs1[sel].get_coverage() >= rUnLibs1[sel].get_needed():
		read = rUnLibs1[sel].get_read()
		rUnLibs1[sel].add_duration(rejTime)
		rUnLibs1[sel].add_coverage(500, read[0])
		untRunT = rejTime
	else:
		untRunT = Pore(rUnLibs1[sel])
	
	return untRunT


def Pore(selection):
#Simulates the passage of sequence through pore
	read = selection.get_read()
	readLen = read[1] - read [0]
	seqTime = (readLen / speed) + interval

	selection.add_coverage(readLen, read[0])
	selection.add_duration(seqTime)
	return seqTime


def Incomplete(libraries):
#This method checks a list of library objects to see if all are complete
	for i in range(0, len(libraries)):
		if libraries[i].get_complete() == 0:
			return 1
		else:
			continue
	
	return 0	#Only reached if no library is incomplete


def Hours(seconds):
#Converts seconds to hours, rounded to nearest second. Returns a string.
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	
	tTot = "{0:.0f}:{1:02d}:{2:02d}".format(h, int(m), int(s))

	return(tTot)


def Results(library):
#Provides a string for reporting on the results of a given library
	timeOut = Hours(library.get_duration())
	cover = library.get_cov_ratio()
	

	output = ("{0}: {1} to complete, {2:.1f}x coverage "\
	.format(library.get_name(), timeOut, cover))

	print(output)

	return(output)


def Graphs(lib, suffix):
#Produces and saves a graph for a given library
	
	data = np.zeros( (1, lib.gsize) )
	
	for i in range(1, len(lib.map[0])-1):
		data[0][i] = data[0][i-1] + lib.map[0][i]
	
	plt.figure(figsize = (8, 6))	
	plt.plot(data[0])
	plt.title(lib.get_name() + suffix)
	plt.ylabel('Coverage')
	plt.xlabel('Position (bp)')
	
	name = (lib.get_name() + suffix)
	plt.savefig(name)


###FUNCTIONS END	

#Define speed of sequencing, interval between sequences, time lost to each 
#rejection and coverage desired
speed = 450	#rate of sequencing
interval = 1	#time taken for a pore to acquire new strand
rejPen = 1	#time taken to reject a strand
#covDes = 30	#########probably outdate#####
idLag = 500	#no. bases needed to map a strand
rejTime = (interval + rejPen + (idLag/speed))

#####
#Create requisite libraries:
simLibs = []	#Libraries for simple experiment
rUnLibs = []	#Libraries for read until experiment
inLibs = []	#Library specifications from input file


try:
	inp = open(sys.argv[1], 'r')
except:
	print("No input file given or file not found")

for line in inp:
	inline = line.split()

	for i in range(0, len(inline)):
		try:
			inline[i] = int(inline[i])
		except:
			pass

	inLibs.append(inline)

inp.close()


for i in range(0, len(inLibs)):
	simLibs.append(lib.Library(inLibs[i][0], inLibs[i][1], \
		inLibs[i][2], inLibs[i][3]))
	
	rUnLibs.append(lib.Library(inLibs[i][0], inLibs[i][1], \
		inLibs[i][2], inLibs[i][3]))

outfile = open("logfile", "w")

#####
#Run the simple experiment
simTotT = 0				#Variable for recording total duration

print("Performing Simple Run:")
while Incomplete(simLibs):
	runTime = SimpleRun(simLibs)
	simTotT += runTime

outfile.write("Non-Read Until Results:\n")
for obj in simLibs:
	output = Results(obj)
	outfile.write(output + "\n")

print("Total run time = {0}".format(Hours(simTotT)))

for i in range(0, len(simLibs)):
	Graphs(simLibs[i], "_no_read_until")


#####
#Run the read until experiment
rUnTotT = 0				#Variable for recording total duration

print("\nPerforming Read Until Run:")
while Incomplete(rUnLibs):
	runTime = ReadUntil(rUnLibs)
	rUnTotT += runTime

outfile.write("\nRead Until Results:\n")
for obj in rUnLibs:
	output = Results(obj)
	outfile.write(output + "\n")

print("Total run time = {0}".format(Hours(rUnTotT)))

for i in range(0, len(rUnLibs)):
	Graphs(rUnLibs[i], "_read_until")

outfile.close()

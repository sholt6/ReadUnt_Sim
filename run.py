#!/usr/bin/python3
#This script simulates MinION experiments in which multiple different 
#libraries are loaded, in different quantities, and where the genomes
#are of different sizes. The overall aim is to try to determine whether read
#until is useful for optimising this experiment. This may be a means of 
#simulating sequencing of a specific chromosome from a mixture of chromosomes


import sys
import library as lib
import random
#import datetime
import time
import numpy as np
import matplotlib.pyplot as plt

start = time.time()

###FUNCTIONS


def Incomplete(libraries):
#This method checks a list of library objects to see if all are complete
	for i in range(0, len(libraries)):
		if libraries[i].get_complete() == 0:
			return 1
		else:
			continue
	
	return 0	#Only reached if no library is incomplete


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


def Select(simLibs2):
	bag = 0 
	for lib in simLibs2:
		entries = lib.gsize * lib.ratio
		bag += entries

	choice = random.randint(0, bag)

	cumulative = 0
	i = 0
	for lib in simLibs2:		
		end = (lib.gsize * lib.ratio) + cumulative
		if end >= choice:
			return i
		else:
			cumulative += end
			i += 1


def Pore(selection):
#Simulates the passage of sequence through pore
	read = selection.get_read()
	readLen = read[1] - read [0]
	seqTime = (readLen / speed) + interval

	selection.add_coverage(readLen, read[0])
	selection.add_duration(seqTime)
	return seqTime


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

def Help():
	print("Usage: python3 run.py <input_file>\n")

	print("Input files should contain parameters for required libraries")
	print("and should be in a tab-separated format as follows:")
	print("<number of bases>   <ratio>   <coverage desired>   <name>\n")

	print("For example:")
	print("5000000	3	30	\"Genome_A\"")
	print("7000000	1	30	\"Genome_B\"\n")



###FUNCTIONS END	

#Define speed of sequencing, interval between sequences, time lost to each 
#rejection and coverage desired
speed = 450	#rate of sequencing
interval = 1	#time taken for a pore to acquire new strand
rejPen = 1	#time taken to reject a strand
idLag = 500	#no. bases needed to map a strand
rejTime = (interval + rejPen + (idLag/speed))

#####
#Create requisite libraries:
simLibs = []	#Libraries for simple experiment
rUnLibs = []	#Libraries for read until experiment
inLibs = []	#Library specifications from input file

#Open input file
if (len(sys.argv) > 2):
	print("Too many arguments given\n")
	Help()
	quit()

try:
	inp = open(sys.argv[1], 'r')
except:
	print("No input file given or file not found")
	Help()

#Open output file
fname = sys.argv[1] + "_results"
outfile = open(fname, "w")
outfile.write("Parameters for this run were:\n")
	
#Read input, add to oupfile for posterity
for line in inp:
	outfile.write(line)
	inline = line.split()

	for i in range(0, len(inline)):
		try:
			inline[i] = int(inline[i])
		except:
			pass

	inLibs.append(inline)

inp.close()
outfile.write("\n")

#Initialise library objects from input
for i in range(0, len(inLibs)):
	simLibs.append(lib.Library(inLibs[i][0], inLibs[i][1], \
		inLibs[i][2], inLibs[i][3]))
	
	rUnLibs.append(lib.Library(inLibs[i][0], inLibs[i][1], \
		inLibs[i][2], inLibs[i][3]))




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
	print(str(obj.get_coverage()) + " bases sequenced")
	outfile.write(output + "\n")

print("\nTotal run time = {0}".format(Hours(simTotT)))

#for i in range(0, len(simLibs)):
#	Graphs(simLibs[i], "_no_read_until")


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
	print(str(obj.get_coverage()) + " bases sequenced")
	outfile.write(output + "\n")

print("Total run time = {0}".format(Hours(rUnTotT)))

#for i in range(0, len(rUnLibs)):
#	Graphs(rUnLibs[i], "_read_until")

end = time.time()

outfile.write("\nThis script took {0:.3f} seconds to complete".format(end-start))
print("\nThis script took {0:.3f} seconds to complete".format(end-start))

outfile.close()

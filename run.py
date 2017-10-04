#!/usr/bin/python3
#This script simulates MinION experiments in which multiple different 
#libraries are loaded, in different quantities, and where the genomes
#are of different sizes. The overall aim is to try to determine whether read
#until is useful for optimising this experiment. This may be a means of 
#simulating sequencing of a specific chromosome from a mixture of chromosomes


import sys
import library as lib
import random
import time
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

start = time.time()

# Defining options:
usage = ("usage: %prog [input_file] args\n\n"
	"Input files should contain parameters for required libraries\n"
	"and should be in a tab-separated format as follows:\n"
	"<number of bases>   <ratio>   <coverage desired>   <name>\n"
	"For example:\n"
	"5000000	3	30	\"Genome_A\"\n"
	"7000000	1	30	\"Genome_B\"\n")

parser = OptionParser(usage = usage)

parser.add_option("-f", "--file", dest="filename",\
	help = "name of tsv for graphing vals", type="string")

parser.add_option("-v", "--var", dest = "var",\
	help = "variable to be recorded, speed by default",\
	type = "string", default = "speed")

parser.add_option("-g", "--graph", dest = "graph",\
	help = "graphs are output only if this flag is specified",\
	action = "store_true", default = False)

parser.add_option("-s", "--speed", dest = "speed",\
	help = "speed of sequencing (b/s), def = 450",\
	type = "int", default = 450)

parser.add_option("-i", "--interval", dest = "interval",\
	help = "interval between strands (s), def = 1",\
	type = "float", default = 1)

parser.add_option("-r", "--rejPen", dest = "rejPen",\
	help = "penalty for rejecting read (s), def = 1",\
	type = "float", default = 1)

parser.add_option("-l", "--idLag", dest = "idLag",\
	help = "bases needed to map strands (b), def = 500",\
	type = "int", default = 500)

(options, args) = parser.parse_args()


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

	#If coverage of the library has been achieved:
	if rUnLibs1[sel].get_coverage() >= rUnLibs1[sel].get_needed():
		read = rUnLibs1[sel].get_read()
		readLen = read[1] - read[0]
		seqTime = (readLen / speed) + interval
		
		#If the read is shorter than the amount needed to id it:
		if readLen < idLag:
			rUnLibs1[sel].add_duration(seqTime)
			rUnLibs1[sel].add_coverage(readLen, read[0])
			untRunT = seqTime

		#Reject the read:
		else:
			rUnLibs1[sel].add_duration(rejTime)
			rUnLibs1[sel].add_coverage(idLag, read[0])
			untRunT = rejTime

	#Sequence the read as normal:
	else:
		untRunT = Pore(rUnLibs1[sel])
	
	#Return time taken to sequence:
	return untRunT


def Select(simLibs2):
#Selects a library to produce a read from
	bag = 0
	for lib in simLibs2:
		entries = lib.gsize * lib.ratio
		bag += entries

	choice = random.randint(0, bag)

	cumulative = 0
	i = 0
	for lib in simLibs2:
		size = lib.gsize * lib.ratio
		cumulative += size

		if cumulative >= choice:
			return i
		else:
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
	if options.graph == False:
		return
	
	data = np.zeros( (1, lib.gsize) )

	for i in range(0, len(data[0])):
		data[0][i] = data[0][i-1]
		if i in lib.map:
			data[0][i] += lib.map[i]
	
	plt.figure(figsize = (8, 6))	
	plt.plot(data[0])
	plt.title(lib.get_name() + suffix)
	plt.ylabel('Coverage')
	plt.xlabel('Position (bp)')
	
	name = (lib.get_name() + suffix)
	plt.savefig(name)
	plt.close()


###FUNCTIONS END	



#Define speed of sequencing, interval between sequences, time lost to each 
#rejection and coverage desired
speed = options.speed		#rate of sequencing - bases/s
interval = options.interval	#time taken for a pore to acquire new strand
rejPen = options.rejPen		#time taken to reject a strand
idLag = options.idLag		#no. bases needed to map a strand
rejTime = (interval + rejPen + (idLag/speed))

varis = {"speed":speed,
	"interval":interval,
	"rejpen":rejPen,
	"idlag":idLag
}


#####
#Create requisite libraries:
simLibs = []	#Libraries for simple experiment
rUnLibs = []	#Libraries for read until experiment
inLibs = []	#Library specifications from input file

#Open input file
if (len(args) > 2):
	print("Too many arguments given\n")
	parser.print_help()

try:
	inp = open(args[0], 'r')
except:
	print("No input file given or file not found")
	parser.print_help()

#Open output file
fname = args[0] + "_results"
outfile = open(fname, "w")
outfile.write("Parameters for this run were:\n")
	
#Read input, add to outfile for posterity
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

outfile.write("\nTotal Simple run time = {0}\n".format(Hours(simTotT)))
print("\nTotal Simple run time = {0}\n".format(Hours(simTotT)))

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
	print(str(obj.get_coverage()) + " bases sequenced")
	outfile.write(output + "\n")

outfile.write("\nTotal Read Until run time = {0}\n".format(Hours(rUnTotT)))
print("\nTotal Read Until run time = {0}\n".format(Hours(rUnTotT)))

for i in range(0, len(rUnLibs)):
	Graphs(rUnLibs[i], "_read_until")

simH = int(simTotT / 3600)
rUnH = int(rUnTotT / 3600)

out = varis[options.var]
if options.filename is not None:
	with open(options.filename, "a") as values:
		values.write("{0}\t{1}\t{2}\n".format(out, simH, rUnH))

#Finishing
end = time.time()

outfile.write("\nThis script took {0:.3f} seconds to complete"\
	.format(end - start))
print("\nThis script took {0:.3f} seconds to complete"\
	.format(end - start))

outfile.close()

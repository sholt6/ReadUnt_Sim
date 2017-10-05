#!/usr/bin/python3
#This Python class will simulate an individual nanopore. Given a list of 
#library objects (see library.py) it will select one to draw a read from. It
#will then simulate sequencing of this read, adding coverage and duration to 
#the originating library as appropriate, and add to its own duration counter.
#All reads are stored in a stack after sequencing

class Npore:
	
	def __init__(self, ):
		self.duration = 0
		self.reads = []
	
	##Experiment Methods
	def Simple(self, libList):
		sel = Select(libList)
		


	def ReadUntil(self, libList):
	
	
	def Select(self, libList):
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

	
	##Read Record Methods
	def ReadGen(self, lib)
		
		self.reads.append([])
		
		newRead = lib.get_read()

		self.reads[-1][0] = lib			#libref
		self.reads[-1][1] = newRead[0]		#read start base
		self.reads[-1][2] = newRead[1]		#read end base
		self.reads[-1][3] = self.duration	#read start time
		self.reads[-1][4] = 0			#read duration

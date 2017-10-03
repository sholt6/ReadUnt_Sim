#!/usr/bin/python3

#This class describes a sequencing library to be put into a specific 
#nanopore experiment. It therefore generates a read on demand for a 
#specific genome, and allows for the libraries to be put into the 
#MinION in ratios other than 1:1


import random
import numpy as np

class Library:

	def __init__(self, gsize, ratio, covDes, name="Unamed Library"):
		self.gsize = gsize
		self.ratio = ratio
		self.covDes = covDes
		self.name = name
		self.coverage = 0
		self.duration = 0
		self.map = {}

	###Read generating method
	def get_read(self):
		start = random.randint(0, self.gsize-1)
		length = np.random.gamma(2.5, 3000., 1)
		end = start + int(length)
		read = [start, end]
		return read

	###Coverage methods
	def add_coverage(self, increment, start):	#Increase coverage
		self.coverage += increment

		try:
			self.map[start] += 1
		except:
			self.map[start] = 1

		end = start + increment + 1
		if end > self.gsize:
			end = self.gsize - 1

		try:
			self.map[end] -= 1
		except:
			self.map[end] = -1


	def get_coverage(self):			#Get current coverage in bases
		return self.coverage

	def get_cov_ratio(self):		#Get coverage as a multiple
		return self.coverage / self.gsize

	def get_needed(self):		#Get required coverage in bases
		return (self.gsize * self.covDes)

	def get_complete(self):			#Check coverage has been hit
		if self.coverage < int(self.get_needed()):
			return 0
		else:
			return 1

	###Time methods
	def add_duration(self, time):
		self.duration += time

	def get_duration(self):
		return self.duration
		#These two are no longer used:
	def get_minutes(self):
		return divmod(self.duration, 60)

	def get_hours(self):
		m = self.duration/60
		return divmod(m, 60)

	###Name methods
	def set_name(self, newName):
		self.name = newName

	def get_name(self):
		return self.name

#!/usr/bin/python3
# This Python class will simulate an individual nanopore. Given a list of
# library objects (see library.py) it will select one to draw a read from. It
# will then simulate sequencing of this read, adding coverage and duration to
# the originating library as appropriate, and add to its own duration counter.
# All reads are stored in a stack after sequencing


class Npore:

    def __init__(self, ):
        self.duration = 0
        self.reads = []

    ## Experiment Methods
    def simple(self, libList):
        sel = Select(libList)
        # Should call Select() on libList; return lib

        # Then call Pore() on selected lib

        # Add duration to read and to self.duration
        # Duration inc. seqTime, interval, any rejPen
        # Need to call add_coverage() Library method

    def readUntil(self, libList):
        # Similar to above, but with read until characteristics

    def select(self, libList):  # Copy-paste from run.py
        bag = 0
        for obj in libList:
            entries = obj.gsize * obj.ratio
            bag += entries

        choice = random.randint(0, bag)

        cumulative = 0
        i = 0
        for obj in libList:
            size = obj.gsize * obj.ratio
            cumulative += size

            if cumulative >= choice:
                return i
            else:
                i += 1

    def pore(self, selection):
        # In run.py as it stands, this is where get_read() is called
        # Perform length calc and time calc; add to Library

    ## Read Record Methods
    def readGen(self, lib)
        
        self.reads.append([])

        newRead = lib.get_read()

        self.reads[-1].append(lib)              # [0] libref
        self.reads[-1].append(newRead[0])       # [1] read start base
        self.reads[-1].append(newRead[1])       # [2] read end base
        self.reads[-1].append(self.duration)    # [3] read start time

        return self.reads[-1]   # Is this a good idea?

    def readRec(self, time):
        # This method should record the duration required by a given read
        # Include its interval and rejPen (where applicable)
        self.reads[-1].append(time)               # [4] read duration

    def prevReadStart(self):
        return self.reads[-1][3]

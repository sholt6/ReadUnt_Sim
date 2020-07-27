#!/usr/bin/python3
# This script simulates MinION experiments in which multiple different
# libraries may be loaded, in different quantities, and where the genomes
# are of different sizes. The overall aim is to try to determine whether read
# until is useful for optimising this experiment. This may be a means of
# simulating sequencing of a specific chromosome from a mixture of chromosomes

import os
import library as lib
import random
import time
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

start = time.time()

def parse_options():
    usage = ("usage: %prog [input_file] args\n\n"
             "Input files should contain parameters for required libraries\n"
             "and should be in a tab-separated format as follows:\n"
             "<number of bases>   <ratio>   <coverage desired>   <name>\n"
             "For example:\n"
             "5000000    3    30    \"Genome_A\"\n"
             "7000000    1    30    \"Genome_B\"\n")

    parser = OptionParser(usage=usage)

    parser.add_option("-f", "--file", dest="filename",
                      help="name of tsv for graphing vals", type="string")

    parser.add_option("-g", "--graph", dest="graph",
                      help="graphs are output only if this flag is specified",
                      action="store_true", default=False)

    parser.add_option("-m", "--map", dest="map",
                      help="store a coverage map of each library",
                      action="store_true", default=False)

    parser.add_option("-s", "--speed", dest="speed",
                      help="speed of sequencing (b/s), def = 450",
                      type="int", default=450)

    parser.add_option("-i", "--interval", dest="interval",
                      help="interval between strands (s), def = 1",
                      type="float", default=1)

    parser.add_option("-r", "--rejPen", dest="rejPen",
                      help="penalty for rejecting read (s), def = 1",
                      type="float", default=1)

    parser.add_option("-l", "--idLag", dest="idLag",
                      help="bases needed to map strands (b), def = 500",
                      type="int", default=500)

    parser.add_option("-c", "--scale", dest="scale",
                      help="scale value for read generator, def = 3000,",
                      type="float", default=3000)

    parser.add_option("-n", "--name", dest="name",
                      help="optional name for this run in .tsv output from -f",
                      type="string", default="NA")

    return parser



# This method checks a list of library objects to see if all are complete
def incomplete(libraries):
    for i in range(0, len(libraries)):
        if libraries[i].get_complete() == 0:
            return 1
        else:
            continue

    return 0    # Only reached if no library is incomplete


# Simulates an experiment without read until
def simple_run(simLibs1, speed, interval):

    sel = select(simLibs1)

    simRunT, readLen = pore(simLibs1[sel], speed, interval)
    # For SimpleRun(), sequenced is identical to readLen but is included for
    # symmetry with ReadUntil()
    sequenced = readLen

    return simRunT, readLen, sequenced


# Simulates an experiment with read until
def read_until_run(rUnLibs1, speed, interval, idLag, rejTime):

    sel = select(rUnLibs1)

    # If coverage of the library has been achieved:
    if rUnLibs1[sel].get_coverage() >= rUnLibs1[sel].get_needed():
        read = rUnLibs1[sel].get_read()
        readLen = read[1] - read[0]
        seqTime = (readLen / speed) + interval

        # If the read is shorter than the amount needed to id it:
        if readLen < idLag:
            rUnLibs1[sel].add_duration(seqTime)
            rUnLibs1[sel].add_coverage(readLen, read[0])
            sequenced = readLen

            try:
                assert(readLen < idLag)
            except:
                print("Read Until method mishandled a >500b read")
                quit()

            return seqTime, readLen, sequenced

        # Reject the read:
        else:
            rUnLibs1[sel].add_duration(rejTime)
            rUnLibs1[sel].add_coverage(idLag, read[0])
            sequenced = idLag

            return rejTime, readLen, sequenced

    # Sequence the read as normal:
    else:
        untRunT, readLen = pore(rUnLibs1[sel], speed, interval)
        sequenced = readLen

    # Return time taken to sequence:
    return untRunT, readLen, sequenced


# Selects a library to produce a read from
def select(simLibs2):
    bag = 0
    for obj in simLibs2:
        entries = obj.gsize * obj.ratio
        bag += entries

    choice = random.randint(0, int(bag))

    cumulative = 0
    i = 0
    for obj in simLibs2:
        size = obj.gsize * obj.ratio
        cumulative += size

        if cumulative >= choice:
            return i
        else:
            i += 1


# Simulates the passage of sequence through pore
def pore(selection, speed, interval):
    read = selection.get_read()
    readLen = read[1] - read[0]
    seqTime = (readLen / speed) + interval

    selection.add_coverage(readLen, read[0])
    selection.add_duration(seqTime)
    return seqTime, readLen


# Converts seconds to hours, rounded to nearest second. Returns a string
def hours(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)

    tTot = "{0:.0f}:{1:02d}:{2:02d}".format(h, int(m), int(s))

    return(tTot)


# Provides a string for reporting on the results of a given library
def results(library):
    timeOut = hours(library.get_duration())
    cover = library.get_cov_ratio()

    output = ("{0}: {1} to complete, {2:.1f}x coverage "
              .format(library.get_name(), timeOut, cover))

    print(output)

    return(output)


# Produces and saves a graph for a given library
def graphs(lib, suffix, graph):
    if graph is False:
        return

    data = np.zeros((1, int(lib.gsize)))

    for i in range(0, len(data[0])):
        data[0][i] = data[0][i-1]
        if i in lib.map:
            data[0][i] += lib.map[i]

    plt.figure(figsize=(8, 6))
    plt.plot(data[0])
    plt.title(lib.get_name() + suffix)
    plt.ylabel('Coverage')
    plt.xlabel('Position (bp)')

    name = (lib.get_name() + suffix)
    plt.savefig(name)
    plt.close()


def main():

    parser = parse_options()
    (options, args) = parser.parse_args()
    # Define speed of sequencing, interval between sequences, time lost to each
    # rejection and coverage desired
    experimentName = options.name
    speed = options.speed          # rate of sequencing - bases/s
    interval = options.interval    # time taken for a pore to acquire new strand
    rejPen = options.rejPen        # time taken to reject a strand
    idLag = options.idLag          # no. bases needed to map a strand
    scale = options.scale          # scale value for read generator
    graph = options.graph
    rejTime = (interval + rejPen + (idLag/speed))

    if options.graph is True:
        options.map = True

    mapReads = options.map

    #####
    # Create requisite libraries:
    simLibs = []    # Libraries for simple experiment
    rUnLibs = []    # Libraries for read until experiment
    inLibs = []     # Library specifications from input file

    # Open input file
    if (len(args) > 2):
        print("Too many arguments given\n")
        parser.print_help()

    try:
        inp = open(args[0], 'r')
    except:
        print("ERROR: No input file given or file not found\n")
        parser.print_help()
        quit()

    # Open output file
    fname = args[0] + "_results"
    outfile = open(fname, "w")
    outfile.write("Parameters for this run were:\n")

    # Read input, add to results file for posterity
    for line in inp:
        outfile.write(line)
        inline = line.split()

        for i in range(0, len(inline)):
            try:
                inline[i] = float(inline[i])
            except:
                pass

        inLibs.append(inline)


    inp.close()
    outfile.write("\n")

    # Initialise library objects from input
    for i in range(0, len(inLibs)):
        simLibs.append(lib.Library(inLibs[i][0], inLibs[i][1],
                       inLibs[i][2], scale, mapReads, inLibs[i][3]))

        rUnLibs.append(lib.Library(inLibs[i][0], inLibs[i][1],
                       inLibs[i][2], scale, mapReads, inLibs[i][3]))


    #####
    # Run the simple experiment
    simTotT = 0                # Variable for recording total duration
    simReads = []              # List of read lengths
    simBases = 0               # Total bases sequenced

    print("Performing Simple Run:")
    while incomplete(simLibs):
        (runTime, read, sequenced) = simple_run(simLibs, speed, interval)
        simTotT += runTime
        simReads.append(read)
        simBases += sequenced

    simAvgRead = int(np.mean(simReads))  # Average read length for simple

    outfile.write("Non-Read Until Results:\n")
    for obj in simLibs:
        output = results(obj)
        print(str(obj.get_coverage()) + " bases sequenced")
        outfile.write(output + "\n")

    outfile.write("\nTotal Simple run time = {0}\n".format(hours(simTotT)))
    print("\nTotal Simple run time = {0}\n".format(hours(simTotT)))

    for i in range(0, len(simLibs)):
        graphs(simLibs[i], "_no_read_until", graph)


    #####
    # Run the read until experiment
    rUnTotT = 0                # Variable for recording total duration
    rUnReads = []              # List of read lengths
    rUnBases = 0               # Total bases sequenced

    print("\nPerforming Read Until Run:")
    while incomplete(rUnLibs):
        (runTime, read, sequenced) = read_until_run(rUnLibs, speed, interval, idLag, rejTime)
        rUnTotT += runTime
        rUnReads.append(read)
        rUnBases += sequenced

    rUnAvgRead = int(np.mean(rUnReads))   # Average read length for read until

    outfile.write("\nRead Until Results:\n")
    for obj in rUnLibs:
        output = results(obj)
        print(str(obj.get_coverage()) + " bases sequenced")
        outfile.write(output + "\n")

    outfile.write("\nTotal Read Until run time = {0}\n".format(hours(rUnTotT)))
    print("\nTotal Read Until run time = {0}\n".format(hours(rUnTotT)))

    for i in range(0, len(rUnLibs)):
        graphs(rUnLibs[i], "_read_until", graph)

    readAvg = np.mean((simAvgRead, rUnAvgRead))
    readAvg = round(readAvg,)

    # Output .tsv if needed
    header = ("Name\tSpeed\tInterval\tRejPen\tIdLag\tSimple.Hours"
              "\tRead.Until.Hours\tSimple.Bases\tRead.Until.Bases\tAvg.Read"
              "\tFold.Change.Hours\tFold.Change.Bases\n")
    simH = int(simTotT / 3600)
    rUnH = int(rUnTotT / 3600)
    fcHours = round((simH / rUnH), 3)
    fcBases = round((simBases / rUnBases), 3)

    if options.filename is not None:
        os.system('touch ' + options.filename)
        if os.stat(options.filename).st_size == 0:
            with open(options.filename, "w") as values:
                values.write(header)

        with open(options.filename, "a") as values:
            values.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}"
                         "\t{7}\t{8}\t{9}\t{10}\t{11}\n"
                         .format(experimentName, speed, interval, rejPen, idLag,
                                 simH, rUnH, simBases, rUnBases, readAvg, fcHours,
                                 fcBases))

    # Finishing
    end = time.time()

    outfile.write("\nThis script took {0:.3f} seconds to complete"
                  .format(end - start))
    print("\nThis script took {0:.3f} seconds to complete"
          .format(end - start))

    outfile.close()

if __name__ == "__main__":
    print("MAIN")
    main()

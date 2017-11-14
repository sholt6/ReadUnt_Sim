# ReadUnt_Sim

## Quick Start

Download run.py and library.py and run from the command line:

    python3 run.py <input file> [OPTS]

An example can be found in the Example folder.

## Read Until Simulator
This script is intended to simulate a nanopore sequencing experiment with and without read until applied, to offer an estimate of the benefit in a given circumstance.
An input file must be provided, to describe the input DNA and its coverage goals. 

### Speed
-s option allows setting of sequencing speed. This is the rate at which DNA is allowed through pores, in bases per second

### Mean Read Length
-c option is used in a gamma distribution function to generate reads but is not itself the mean read length. The script readgen.py can assist in determining an appropriate value

### Interval
-i is the time taken to acquire a new read after one finishes, in seconds

### Rejection Penalty
-r is the time penalty incurred in order to reject a read, in seconds

### Identification lag
-l is the amount of bases which must be sequenced in order for read until to be able to map it to the reference and make a decision on whether to reject

## Input Files
The input file must be a tab-separated value file in which each line describes something in the input library. These may be individual chromosomes, whole genomes, or anything else required. The columns are as follows:

    Size in base pairs    Relative Copy Number      Desired Coverage        Name

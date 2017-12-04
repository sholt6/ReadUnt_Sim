# ReadUnt_Sim
For assistance, contact Sam Holt at sam_holt30@hotmail.co.uk
## Quick Start

Download run.py and library.py and run from the command line:

    python3 run.py <input file> [OPTS]

Examples can be found in example.md.

## Requirements
The following modules must be installed:
- numpy
- matplotlib

## Read Until Simulator
This script is intended to simulate a nanopore sequencing experiment with and without read until applied, to offer an estimate of the benefit in a given circumstance.
An input file must be provided, to describe the input DNA and its coverage goals. 
The hours values output by the script **should not be taken as accurate estimates of the duration of a given experiment**. This is because the script only simulates a single nanopore.

### Speed
-s option allows setting of sequencing speed. This is the rate at which DNA is allowed through pores, in bases per second.

### Mean Read Length
-c option is used in a gamma distribution function to generate reads but is not itself the mean read length. The script readgen.py can assist in determining an appropriate value for a desired mean read length.

### Interval
-i is the time taken to acquire a new read after one finishes, in seconds.

### Rejection Penalty
-r is the time penalty incurred in order to reject a read, in seconds.

### Identification lag
-l is the amount of bases which must be sequenced in order for read until to be able to map it to the reference and make a decision on whether to reject.

## Input Files
The input file must be a tab-separated value file in which each line describes something in the input library. These may be individual chromosomes, whole genomes, or anything else required. The columns are as follows:

    Size in base pairs    Relative Copy Number      Desired Coverage        Name

## Other Files In This Repository
There are other Python scripts in this repository besides run.py and library.py. These include:
- average.py
  - Accepts a list of .tsv files of same shape and size, and attempts to average corresponding values.
- npore.py
  - This was intended to be used in a more detailed simulation which would have included multiple nanopores. It is incomplete.
- readgen.py 
  - Can be used to help determine what scale value should be used in the run.py -c option to give a desired mean read length
- replace.py
  - This was used in the results of some metagenomics experiments to replace genome percentage values with base percentage values. It is not written for generalised usage
- specimen_plot.ipynb
  - A Jupyter Notebook which was used as a base from which to create specific plots of results

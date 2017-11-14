# Read Until Simulation Example
This is intended as an example demonstration of Read Until Simulator

Please download run.py, library.py, and input\_original from the [ReadUnt_Sim repository](https://github.com/sholt6/ReadUnt_Sim/tree/develop). These should be placed in the same directory and run using Python 3

## Example 1
The file input\_original describes an input of two prokaryotic genomes, for species A and B. Species A has a genome of 5Mb, B of 7Mb. Copies of genome A are present at three times the rate of genome B. The goal is to sequence both genomes to 30x coverage. 
To run this experiment with default settings:

    python3 run.py input\_original 

The result should be a description of the amount of time taken to sequence each genome for the non-read until experiment, along with the coverage achieved and number of bases sequenced. This should then be followed by the same information for the read until experiment.
To view the coverage across each genome, add in the -g flag:

    python3 run.py input\_original -g

Once the script completes, the graphs should be present in the working directory.

## Example 2
It may now be desirable to modify the sequencing parameters. Adding -s <integer> sets the speed to a value other than 450. We will try two alternative values, and output the results to example2\_results.tsv in order to more easily compare them. The -n flag gives each run a name, allowing it to be identified in the results file.

    python3 run.py input\_original -f example2_results.tsv -s 400 -n 400
    python3 run.py input\_original -f example2_results.tsv -s 500 -n 400

The file example\_2results.tsv contains the duration and bases sequenced for both of these runs. They can be identified by the Names column, and the Fold.Change.Hours column gives us a measure of how many times faster (or slower!) the experiment was with read until than without.

## Example 3
In the repository (see above for link) can be found input\_human\_genome\_chr\_21\_30x. This allows simulation of a human genome sequencing experiment. Each human chromosome (1-21, X, Y) is present at equal copy number. The goal is set at 30x coverage of chromosome 21.

    python3 run.py input\_human\_genome\_chr\_21\_30x

This run can be expected to take considerably longer. You may like to vary which chromosome is enriched, and observe how the size of the target influences the performance of read until.

#!/usr/bin/python3

# This script imports two identically formatted .tsv files to Pandas dataframes
# and prints a new .tsv in which each cell contains the average of all the input
# files

import sys
import pandas as pd

# Check args have been provided
args = sys.argv[1:]

if not args:
    print("Usage: python3 average.py <tsv1> <tsv2> ...")
    quit()

# Import data
data = []

for fname in args:
    df = pd.read_csv(fname, sep='\t')
    data.append(df)

# Assert more than one dataframe has been put in
try:
    assert(len(data) > 1)
except:
    print("Please specify at least two files")

# Check dataframes are of same structure
for i in range(1, len(data)):
    try:
        assert(data[i].shape == data[i-1].shape)
    except:
        print("Dataframes are not all of the same shape")
        quit()

    for j in range(0, len(data[i].columns)):
        try:
            assert(data[i].columns[j] == data[i-1].columns[j])
        except:
            print("Dataframe columns do not match")
            quit()



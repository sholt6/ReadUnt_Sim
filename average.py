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

# At this point we have a list of appropriately formatted dataframes to average
# Concatenate dataframes and produce means
df_concat = pd.concat(data)
by_row_index = df_concat.groupby(df_concat.index)
df_means = by_row_index.mean()

# If Names column contains strings, they are absent from df_means.
# This block restores them
try:
    df_means['Name']
except KeyError:
    df_means.insert(0, 'Name',
                    pd.Series(data[0]['Name'].values, index=df_means.index))

# Output to new .tsv file
df_means.to_csv('RENAME_ME.tsv', sep='\t', index=False)

#!/usr/bin/python3
# For a results tsv from experiments using my super_simple_metagenome this
# will convert the 'Names' column from a percentage copy number to a percentage
# bases

import pandas as pd
import sys

def SimpleGenomePC(hpc):
    epc = 100 - hpc
    
    basesTotal = (4630707 * epc) + (3234830000 * hpc)
    basesHuman = 3234830000 * hpc
    propHuman = basesHuman / basesTotal
    
    try:
        assert(propHuman > 0)
        assert(propHuman < 1)
        return propHuman
    except:
        print("Invalid human genome percentage received")

files = sys.argv[1:]

for i in range(0, len(files)):
    newDf = pd.read_csv(files[i], sep="\t")

    for val in range(0, len(newDf['Name'])):
        hpc = newDf.loc[val,'Name']
        newVal = SimpleGenomePC(hpc)
        newDf.loc[val,'Name'] = newVal

    splitFileName = files[i].split('.')
    newFileName = splitFileName[0] + "_base_pcs." + splitFileName[1]
    newDf.to_csv(newFileName, sep='\t', index=False)

#!/usr/bin/env python
import pandas as pd
import glob
import os

## Set directory with counts data
os.chdir("../result/4_count")

## List all counts files
extension = 'tsv'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]

## Combine all files
### Blank df to hold results
#index=[]
#all_tsv = pd.DataFrame({'Geneid' : []})
all_tsv = pd.DataFrame(index=[])
for f in all_filenames:
    ### Skip combined file if already exists
    if not f.startswith('combined_feature_counts'):
        ### Read in tsv
        temp = pd.read_csv(f, delimiter="\t", skiprows=1, index_col='Geneid')
        ### Remove gene info columns
        temp = temp.drop(['Chr','Start','End','Strand','Length'], axis=1)
        ### Combine
        all_tsv = all_tsv.join(temp, how='outer')

## Clean sample names
all_tsv.columns = all_tsv.columns.str.replace('result/3_bam_filter/', '', regex=True)
all_tsv.columns = all_tsv.columns.str.replace('_Aligned.sortedByCoord.filter.bam', '', regex=True)

## Save
all_tsv.to_csv("combined_feature_counts.tsv", index=True, encoding='utf-8-sig', sep="\t")

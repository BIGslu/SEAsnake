#!/usr/bin/env python
import pandas as pd
import glob
import os

##################################
## Counts table
##################################

## List all counts files
extension = 'tsv'
all_filenames = [i for i in glob.glob('result/4_count/*.{}'.format(extension))]

## Combine all files
### Blank df to hold results
all_tsv = pd.DataFrame(index=[])

for f in all_filenames:
    print(f)
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
all_tsv.to_csv("result/5_combined/combined_feature_counts.tsv", index=True, encoding='utf-8-sig', sep="\t")


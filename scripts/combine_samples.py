#!/usr/bin/env python
import pandas as pd
import glob
import os
extension = 'tsv'

##################################
## Counts table
##################################

## List all counts files
all_count_file = [i for i in glob.glob('result/4_count/*.{}'.format(extension))]

## Combine all files
### Blank df to hold results
all_count = pd.DataFrame(index=[])

for f in all_count_file:
    print(f)
    ### Read in tsv
    temp = pd.read_csv(f, delimiter="\t", skiprows=1, index_col='Geneid')
    ### Remove gene info columns
    temp = temp.drop(['Chr','Start','End','Strand','Length'], axis=1)
    ### Combine
    all_count = all_count.join(temp, how='outer')

## Clean sample names
all_count.columns = all_count.columns.str.replace('result/3_bam_filter/', '', regex=True)
all_count.columns = all_count.columns.str.replace('_Aligned.sortedByCoord.filter.bam', '', regex=True)

## Save
all_count.to_csv("result/5_combined/combined_feature_counts.tsv", index=True, encoding='utf-8-sig', sep="\t")

##################################
## Flagstat
##################################
## List all flagstat files
all_stat_file = [i for i in glob.glob('result/qc/3_flagstat/*.{}'.format(extension))]

## Combine all files
### Blank df to hold results
all_stat = pd.DataFrame(index=[])

if len(all_stat_file) > 0:
    for f in all_stat_file:
        print(f)
        ### Read in tsv
        temp = pd.read_csv(f, sep='+', header=None, names=['value',2,3])
        ### Create data ID column
        temp = temp.assign(index = ['QC_pass', 'primary', 'secondary', 'supplementary', 'duplicate', 'primary_duplicate', 'mapped', 'primary_mapped', 'paired', 'read1', 'read2', 'paired_proper', 'paired_mapped', 'singleton', 'paired_diff_chr', 'paired_diff_chr5'])
        temp = temp.set_index('index')
        temp = temp.filter(items=['value'])
        ### Transpose
        temp = temp.transpose()
        ### Get library ID
        libID = f.replace('result/qc/3_flagstat/', '', 1)
        libID = libID.replace('_flagstat.tsv', '', 1)
        temp = temp.assign(libID = libID)
        ### Combine
        all_stat = pd.concat([all_stat,temp])
    ## reorder columns
    cols = list(all_stat.columns)
    cols.pop(cols.index('libID')) #Remove libID from list
    all_stat = all_stat[['libID']+cols]
    ## Save
    all_stat.to_csv("result/5_combined/combined_flagstat.tsv", index=False, encoding='utf-8-sig', sep="\t")
else:
    no_dat = pd.DataFrame(message=["samtools flagstat not completed."])
    no_dat.to_csv("result/5_combined/combined_flagstat.tsv", index=False, encoding='utf-8-sig', sep="\t")

##################################
## Picard
##################################
## List all picard files
all_picard_file = [i for i in glob.glob('result/qc/4_picard/*.{}'.format(extension))]

## Combine all files
### Blank df to hold results
all_picard = pd.DataFrame(index=[])


if len(all_picard_file) > 0:
    for f in all_picard_file:
        print(f)
        ### Read in tsv
        temp = pd.read_csv(f, delimiter="\t", skiprows=6)
        ### Remove histogram data
        temp = temp.iloc[[0]]
        ### Get library ID
        libID = f.replace('result/qc/4_picard/', '', 1)
        libID = libID.replace('_picard.tsv', '', 1)
        temp = temp.assign(libID = libID)
        ### Combine
        all_picard = pd.concat([all_picard,temp])
    ## reorder columns
    cols = list(all_picard.columns)
    cols.pop(cols.index('libID')) #Remove libID from list
    all_picard = all_picard[['libID']+cols]
    ## Save
    all_picard.to_csv("result/5_combined/combined_picard.tsv", index=True, encoding='utf-8-sig', sep="\t")
else:
    no_dat = pd.DataFrame(message=["picard not completed."])
    no_dat.to_csv("result/5_combined/combined_picard.tsv", index=False, encoding='utf-8-sig', sep="\t")

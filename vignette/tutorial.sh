#!/bin/sh

#  tutorial.sh
#  
#
#  Created by Kim Dill-McFarland on 12/14/21.
#  


nohup snakemake --snakefile Snakefile_step1 --cores 15 > log/SEAsnake_step1.log 2>&1

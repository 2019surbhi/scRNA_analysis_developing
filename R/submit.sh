#!/bin/bash 
# INPUTS: 
# 1) Name of sbatch script in "scripts" folder to submit
# 2) Name of output folder

## e.g.
# ```
# ./submit.sh submit.pipe.sb ~/beegfs/mouse-bladder/scrna/seurat-results/mm10/subsets/10-8-uro-cleaned-matt-dev/

# Standard name for run should be [date]-[data]-[pipeline]
# e.g. "10-7-immune-harmony"
cd scripts/
export jid=$(sbatch $1 $2 | grep -o "[0-9]*")
sed "s/FOO/$jid/g" document.sb > document_${jid}.sb
sbatch document_${jid}.sb $2
rm document_${jid}.sb 

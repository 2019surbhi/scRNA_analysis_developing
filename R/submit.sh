#!/bin/bash 
# Standard name for run should be [date]-[data]-[pipeline]
# e.g. "10-7-immune-harmony"
cd scripts/
export out_dir="/home/bradlem4/beegfs/mouse-bladder/scrna/seurat-results/mm10/subsets/doc-test/"
export jid=$(sbatch $1 $out_dir | grep -o "[0-9]*")
sed "s/FOO/$jid/g" document.sb > document_${jid}.sb
sbatch document_${jid}.sb $out_dir
rm document_${jid}.sb 

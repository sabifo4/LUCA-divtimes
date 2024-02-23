#!/bin/bash

# Get args
aln=$1     # 1, 2, etc.
pipedir=$2 # Path to pipeline dir
name_wd=$3 # Name of the working directory, e.g., `euk110`

# Replace vars in template bash script for job array
cp pipeline_Hessian_CODEML_template.sh $pipedir/pipeline_Hessian.sh
if [[ $aln -eq 1 ]]
then 
sed -i 's/\#\$\ \-t\ ..*/\#\$\ \-t\ 1/' $pipedir/pipeline_Hessian.sh
else 
sed -i 's/NUM/'${aln}'/' $pipedir/pipeline_Hessian.sh
fi
# Replace name of working directory
upd_wd=$( echo $name_wd | sed 's/\//\\\//g' | sed 's/_/\\_/g' )
sed -i 's/WDNAME/'${upd_wd}'/g' $pipedir/pipeline_Hessian.sh


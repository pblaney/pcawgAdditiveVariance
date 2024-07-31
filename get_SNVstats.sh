#!/bin/bash

echo "###########################################################"
echo "#                        getSNVstats                      #"
echo "###########################################################"
echo 

# Load modules
module add matlab/R2023a

# User defines if null or obs dataset to run
nullOrObs=$1
matScript=$(echo "get_SNVstats_${nullOrObs}.m")
matLog=$(echo "get_SNVstats_${nullOrObs}.log")

echo 
echo "Null or Obs ... ${nullOrObs}"
echo 

# Execute
matCmd="
matlab -nodisplay -nosplash -nojvm < ${matScript} > ${matLog}
"
echo "CMD: ${matCmd}"
eval "${matCmd}"

echo 
echo "D O N E ..."
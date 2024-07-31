#!/bin/bash

echo "###########################################################"
echo "#                     additiveVariance                    #"
echo "###########################################################"
echo 

# Load modules
module add matlab/R2023a

# Execute
matCmd="
matlab -nodisplay -nosplash -nojvm < additive_variance.m > additive_variance.log
"
echo 
echo "CMD: ${matCmd}"
eval "${matCmd}"

echo
echo "D O N E ..."
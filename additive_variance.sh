#!/bin/bash

echo "###########################################################"
echo "#                     additiveVariance                    #"
echo "###########################################################"
echo 

# Load modules
module add matlab/R2023a

# make necessary directories if not made
matlab -nodisplay -nosplash -nojvm < additive_variance.m > additive_variance.log


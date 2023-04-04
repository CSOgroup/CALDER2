#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

TEST_ENV_NAME="R-Test"

mamba create --yes --name ${TEST_ENV_NAME} --channel bioconda --channel conda-forge \
		r-essentials r-devtools cooler r-nloptr
conda activate ${TEST_ENV_NAME}

R -e "devtools::install('.', dependencies=TRUE, Ncpus = 20)"
R -e "devtools::test()"

# conda deactivate 
mamba env remove --yes --name ${TEST_ENV_NAME}

#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

TEST_ENV_NAME="R-Test"

mamba create --yes --name ${TEST_ENV_NAME} --channel bioconda \
		r-essentials r-devtools bioconductor-biocinstaller bioconductor-genomicranges \
		r-nloptr
conda activate ${TEST_ENV_NAME}

R -e "devtools::install('.', dependencies=TRUE)"

conda deactivate 
mamba env remove --yes --name ${TEST_ENV_NAME}

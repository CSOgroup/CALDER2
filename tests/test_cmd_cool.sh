#!/bin/bash

mkdir -p "tests/output"

scripts/calder --input tests/data/test.cool \
			   --type cool \
			   --bin_size 50000 \
			   --genome hg38 \
			   --nproc 10 \
			   --outpath "tests/output/test_cmd_cool_out"
#!/bin/bash


create_test_cool(){
	url=$1
	outpath=$2
	binsize=$3

	if [[ -f ${outpath}/test.cool ]]; then return; fi
	if [[ ! -f ${outpath}/source.pairs.gz ]] 
	then
		echo "[Error] Please manually download the file ${url} to ${outpath}"
		echo "        and give it source.pairs.gz as name"
		exit -1
	fi

	zcat ${outpath}/source.pairs.gz | head -26 | tail -2 | cut -d' ' --output-delimiter=$'\t' -f 2,3 > ${outpath}/test.chrom.sizes
	cooler cload pairs --chrom1 2 --pos1 3 --chrom2 4 --pos2 5 \
						${outpath}/test.chrom.sizes:${binsize} \
						${outpath}/source.pairs.gz \
						${outpath}/test.cool
	rm ${outpath}/source.pairs.gz

	cooler balance --force --max-iters 1000 ${outpath}/test.cool
}


create_test_mcool(){
	cool_path=$1
	outpath=$2
	binsize=$3

	if [[ ! -f ${outpath} ]]
	then
		cooler zoomify --resolutions ${binsize}N \
					   --balance \
					   --balance-args '--force --max-iters 1000' \
					   --out ${outpath} ${cool_path}
	fi
}


test_data_path="tests/testthat/data"

mkdir -p ${test_data_path}

# Test cool file
source_cool_file="https://data.4dnucleome.org/files-processed/4DNFI2EK1IOQ/@@download/4DNFI2EK1IOQ.pairs.gz"
test_cool_binsize=50000

create_test_cool ${source_cool_file} ${test_data_path} ${test_cool_binsize}

# Creating MCOOL from COOL
create_test_mcool ${test_data_path}/test.cool ${test_data_path}/test.mcool ${test_cool_binsize}
#!/bin/sh
vcffile=$1
hlafile=$2
outdir=$3
export PATH=/lustre1/shiyang/software/anaconda2/bin:${PATH}
NETMHC=/lustre1/shiyang/software/netMHCpan-4.0/netMHCpan
CACHE_DIR=/gpfs1/ruibinxi_pkuhpc/
neosv=/lustre1/shiyang/software/neosv/neosv.py
python ${neosv} \
	   -vf ${vcffile} \
	   -hf ${hlafile} \
	   -np ${NETMHC} \
	   -o ${outdir} \
	   -r 75 -pd ${CACHE_DIR}


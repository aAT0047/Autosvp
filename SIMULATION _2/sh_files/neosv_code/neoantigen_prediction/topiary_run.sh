#!/bin/sh
maffile=$1
hlafile=$2
outfile=$3
export PATH=/lustre1/shiyang/software/anaconda2/bin:${PATH}
export PATH=/lustre1/shiyang/software/netMHCpan-4.0:${PATH}
export PYENSEMBL_CACHE_DIR=/gpfs1/ruibinxi_pkuhpc/
topiary --maf ${maffile} \
		--mhc-alleles-file ${hlafile} \
		--genome GRCh37 \
		--mhc-predictor netmhcpan \
		--ic50-cutoff 500 \
		--percentile-cutoff 2.0 \
		--mhc-epitope-lengths 8-11 \
		--output-csv ${outfile}

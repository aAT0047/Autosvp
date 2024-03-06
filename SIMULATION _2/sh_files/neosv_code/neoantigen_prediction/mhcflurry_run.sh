#!/bin/bash

INFILE=$1
OFILE=$2

source /lustre1/shiyang/python2_to_python3_new.sh
source activate mhcflurry_new
mhcflurry-predict ${INFILE} --out ${OFILE}
source deactivate
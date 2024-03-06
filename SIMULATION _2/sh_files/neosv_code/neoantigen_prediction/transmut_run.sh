#!/bin/bash

PEPFILE=$1
FAFILE=$2
OUTDIR=$3

source /lustre1/shiyang/python2_to_python3_new.sh
source activate transmut
cd /home/ruibinxi_pkuhpc/lustre1/shiyang/software/TransPHLA-AOMP-master/TransPHLA-AOMP
python pHLAIformer.py --peptide_file ${PEPFILE} --HLA_file ${FAFILE} --output_dir ${OUTDIR}
source deactivate
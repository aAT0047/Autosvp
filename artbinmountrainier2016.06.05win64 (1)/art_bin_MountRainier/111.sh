#!/bin/bash
#illumina test examples
art=art_illumina
#art=../../bin/MacOS64/art_illumina
#art=../../bin/Linux64/art_illumina

# 1) simulation of single-end reads of 35bp with 10X using the built-in combined quality profile, and without Indels
$art  -i frandom_sequence_1.fasta -p -l 150 -f 50 -m 400 -s 50 -o C:\Users\1\Desktop\code_ppg\artbinmountrainier2016.06.05win64 (1)\art_bin_MountRainier\test\simulated_data


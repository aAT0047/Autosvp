#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
# import argparse

# # 创建命令行解析器
# parser = argparse.ArgumentParser(description="将输入文件处理为VCF格式并将结果输出到指定文件")

# # 添加输入文件路径参数
# parser.add_argument("input_file", type=str, help="输入文件的路径")

# # 添加输出文件路径参数
# parser.add_argument("output_file", type=str, help="输出文件的路径")

# # 解析命令行参数
# args = parser.parse_args()

# Define VCF file's metadata and title line
vcf_header = """##fileformat=VCFv4.2
##source=YourProgramName
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=CNV,Description="Copy number variable region">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""

# Map list types to VCF format's ALT field
alt_mapping = {
    'onesnv':'<SNV>',
    'onedelete': '<DEL>',
    'oneinsert': '<INS>',
    'onetandem_repeat': '<DUP>',
    'oneinversion': '<INV>'
}

# 从输入文件读取数据
with open(r'C:\Users\1\Desktop\code_ppg\sv_sh\sh_files\A_stableCallerPaperSimFlowShell_1.sh') as answerObject:
    answerFile = answerObject.read()

answerFile_rows = answerFile.split('\n')
answerFile_arr = answerFile_rows[0].split(' -')

delAnswerArr = answerFile_arr[5].split(' ')[1].split(':')
# insAnswerArr = answerFile_arr[6].split(' ')[1].split(':')
trAnswerArr = answerFile_arr[6].split(' ')[1].split(':')
invAnswerArr = answerFile_arr[7].split(' ')[1].split(':')
# snvAnswerArr = answerFile_arr[9].split(' ')[1].split(':')
if len(answerFile_arr) > 8:
    parts = answerFile_arr[8].split(' ')
    
    if len(parts) > 1:
        snvAnswerArr = parts[1].split(':')
    else:
        snvAnswerArr = []  # 或其他默认值
else:
    snvAnswerArr = []  # 或其他默认值
allAnswerArr = snvAnswerArr+delAnswerArr  + trAnswerArr + invAnswerArr
your_list = allAnswerArr
# print(your_list)
# 打开输出文件并写入VCF数据
:with open(args.output_file, 'w') as f:
    f.write(vcf_header + '\n')
    for item in your_list:
        fields = item.split(',')
        svtype = fields[0]
        chrom = fields[1]
        pos = fields[2]
        
        if svtype == 'onesnv':
            replaced_nucleotide = fields[4]
            info_field = "REPLACED_BASE={}".format(replaced_nucleotide)
        elif svtype in ['onedelete', 'oneinversion', 'onetandem_repeat']:
            end = fields[3]
            info_field = "END={}".format(end)
        else:
            info_field = "."

        id_field = '.'
        ref = '.'
        qual = '.'
        filter_field = '.'
        alt = alt_mapping.get(svtype, '.')

        vcf_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrom, pos, id_field, ref, alt, qual, filter_field, info_field)
        f.write(vcf_line + '\n')





print("VCF文件已写入 {}".format(args.output_file))

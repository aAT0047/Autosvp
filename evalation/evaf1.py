#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os

def extract_variant_info(vcf_path):
    variant_data = []
    
    with open(vcf_path, 'r') as f:
        for line in f:
            # Skip metadata and header lines
            if line.startswith("#"):
                continue
            
            fields = line.strip().split('\t')
            
            # Ensure the line has enough fields
            if len(fields) < 8:
                # print("Unexpected line format: {}".format(line))
                continue
            
            # Chromosome and position
            chrom = fields[0]
            pos = fields[1]
            
            # Variant type (from ALT column)
            variant_type = fields[4]
            
            # Extract END from INFO field
            info_fields = fields[7].split(';')
            end_position = None
            for info in info_fields:
                if info.startswith("END="):
                    end_position = info.split('=')[1]
                    break
            
            variant_data.append((chrom, pos, end_position, variant_type))
    
    return variant_data


def sort_vcf(input_path, output_path):
    # Extract header, metadata, and variants
    headers, variants = [], []
    with open(input_path, 'r') as f:
        for line in f:
            if line.startswith("##"):
                headers.append(line.strip())
            elif line.startswith("#"):
                column_headers = line.strip()
            else:
                variants.append(line.strip())

    # Sort variants by chromosome and then by position
    sorted_variants = sorted(variants, key=lambda x: (x.split('\t')[0], int(x.split('\t')[1])))

    # Write to output
    with open(output_path, 'w') as f:
        for header in headers:
            f.write(header + "\n")
        f.write(column_headers + "\n")
        for variant in sorted_variants:
            f.write(variant + "\n")

def is_matching_variant(var1, var2):
    _, start1, end1, type1 = var1
    # print(start1, end1, type1)
    _, start2, end2, type2 = var2

    # Convert to integer if not None, otherwise keep as None
    start1 = int(start1) if start1 is not None else None
    end1 = int(end1) if end1 is not None else None
    start2 = int(start2) if start2 is not None else None
    end2 = int(end2) if end2 is not None else None

    # Check conditions with added None checks
    start_condition = abs(start1 - start2) <= 500 if start1 is not None and start2 is not None else False
    end_condition = abs(end1 - end2) <= 500 if end1 is not None and end2 is not None else False
    # print(start_condition and end_condition and type1 == type2)
    return start_condition or end_condition 
        # and type1 == type2


def compute_f1(standard_vcf, called_vcf):
    TP = 0
    FP = 0
    FN = 0
    
    matched_indices = set()  # To keep track of matched variants in called_vcf

    for std_var in standard_vcf:
        matched = False
        for idx, call_var in enumerate(called_vcf):
            if is_matching_variant(std_var, call_var):
                matched = True
                matched_indices.add(idx)
                break
        if matched:
            TP += 1
        else:
            FN += 1

    FP = len(called_vcf) - len(matched_indices)

    # 使用浮点数除法
    precision = float(TP) / (TP + FP) if (TP + FP) != 0 else 0
    # 使用浮点数除法
    recall = float(TP) / (TP + FN) if (TP + FN) != 0 else 0

    print(precision)
    print( recall)
    if precision + recall == 0:  # To handle the case when both precision and recall are zero
        return 0.0
    else:
        f1 = 2.0 * (precision * recall) / (precision + recall)
        return f1


# 导入你的 sort_vcf、extract_variant_info 和 compute_f1 函数


def main():
    parser = argparse.ArgumentParser(description="计算两个VCF文件的F1分数")
    parser.add_argument("input_vcf_1", help="第一个输入VCF文件的路径")
    parser.add_argument("input_vcf_2", help="第二个输入VCF文件的路径")
    parser.add_argument("output_f1_score", help="输出F1分数的文件路径")

    args = parser.parse_args()

    # 排序输入VCF文件
    sort_vcf(args.input_vcf_1, args.input_vcf_1)
    sort_vcf(args.input_vcf_2, args.input_vcf_2)
    # 提取变异信息
    variant_info_1 = extract_variant_info(args.input_vcf_1)
    variant_info_2 = extract_variant_info(args.input_vcf_2)
    # print(variant_info_1)
    # print(variant_info_2)
    # 计算F1分数
    f1_score = compute_f1(variant_info_1, variant_info_2)
    print(f1_score )
    # 将F1分数写入输出文件
# 将F1分数写入输出文件，显示四位小数
    with open(args.output_f1_score, "w") as output_file:
        output_file.write("%.8f" % f1_score)


if __name__ == "__main__":
    main()

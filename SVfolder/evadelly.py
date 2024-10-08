import os
import csv
from evaf1 import *

def get_vcf_files(input_dir):
    """从目录中获取所有 .vcf 文件"""
    vcf_files = [f for f in os.listdir(input_dir) if f.endswith('.vcf')]
    return vcf_files

def parse_vcf(vcf_file):
    """解析 VCF 文件，提取 TRA 变异并找到证据最多的一个"""
    tra_variants = []
    
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue  # 跳过头部信息
            columns = line.strip().split("\t")
            info = columns[7]
            if "SVTYPE=TRA" in info:  # 查找TRA变异
                pe = int([x.split('=')[1] for x in info.split(';') if x.startswith('PE')][0])
                tra_variants.append((columns, pe))
    
    # 如果有TRA变异，选择PE最大的一个
    if tra_variants:
        best_tra = max(tra_variants, key=lambda x: x[1])
        return best_tra[0]
    
    return None

def format_to_output(vcf_data, event_id_base):
    """将提取的VCF数据格式化成指定的输出格式"""
    chrom1, pos1, id1, ref, alt, qual, filt, info, fmt, sample = vcf_data
    info_dict = dict(item.split("=") for item in info.split(";") if "=" in item)
    
    chr2 = info_dict.get("CHR2", "")
    end = info_dict.get("END", "")
    
    # 动态生成事件 ID
    event_id_1 = f"{event_id_base}"
    event_id_2 = f"{event_id_base}"
    
    output = []
    output.append(f"{chrom1}\t{pos1}\t{event_id_1}\tN\t[{chr2}:{end}[N\t.\tPASS\tSVTYPE=SV;STRANDS=--:503;EVENT={event_id_base}_BCL2_1;MATEID={event_id_2};CIPOS=0,0;CIEND=0,0;CIPOS95=0,0;CIEND95=0,0;SU=14.19;PE=0;SR=5238.35;END={end}\tGT:SU:PE:SR\t./.:14.19:0:5238.35")
    output.append(f"{chr2}\t{end}\t{event_id_2}\tN\tN[{chrom1}:{pos1}]\t.\tPASS\tSVTYPE=SV;STRANDS=--:503;EVENT={event_id_base}_BCL2_1;MATEID={event_id_1};CIPOS=0,0;CIEND=0,0;CIPOS95=0,0;CIEND95=0,0;SU=14.19;PE=0;SR=5238.35;END={pos1}\tGT:SU:PE:SR\t./.:14.19:0:5238.35")
    
    return "\n".join(output)

def process_vcf_files(input_dir, output_dir):
    """处理所有 VCF 文件并返回修改后的文件路径和对应的事件 ID"""
    vcf_files = get_vcf_files(input_dir)
    modified_files = []
    
    for vcf_file in vcf_files:
        vcf_path = os.path.join(input_dir, vcf_file)
        
        # 提取 VCF 文件名作为事件 ID 的基础部分
        event_id_base = os.path.splitext(os.path.basename(vcf_file))[0]
        
        best_tra = parse_vcf(vcf_path)
        
        if best_tra:
            formatted_output = format_to_output(best_tra, event_id_base)
            
            # 输出文件名
            output_vcf = os.path.join(output_dir, f"{event_id_base}_formatted.vcf")
            
            # 将格式化后的内容写入输出文件
            with open(output_vcf, 'w') as out_file:
                out_file.write(formatted_output)
            
            # 添加到返回列表
            modified_files.append((output_vcf, event_id_base))
            print(f"Processed {vcf_file} -> {output_vcf}")
    
    return modified_files

# 新增函数用于将评估结果写入 CSV 文件
def write_evaluation_to_csv(csv_file, evaluation_data):
    """将评估数据写入 CSV 文件"""
    with open(csv_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['File Name', 'f1_score', 'precision', 'recall'])  # CSV的表头
        for row in evaluation_data:
            writer.writerow(row)

# 主函数执行流程
input_directory = "/data/home/std_12/ShiHe/results/sarcomabam/delly"  # 输入文件夹路径
output_directory = "/data/home/std_12/ShiHe/results/sarcomabam/delly_output_vcfs"  # 输出文件夹路径
output_truth = "/data/home/std_12/ShiHe/sarcomavcf"  # 用于 evaf1 的真值文件

# 确保输出目录存在
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# 处理所有 VCF 文件并获取修改后的文件路径和事件ID列表
modified_vcf_files = process_vcf_files(input_directory, output_directory)

# 准备CSV文件
csv_file_path = os.path.join(output_directory, "evaluation_results.csv")

# 初始化评估结果列表
evaluation_data = []

# 遍历所有生成的 VCF 文件，并对每个文件进行评估
for output_vcf, event_id_base in modified_vcf_files:
    # 真值文件路径
    output_vcf_path = os.path.join(output_truth, f"{event_id_base}.vcf")
    
    # 使用 evamain 函数进行评估
    TP, FP, FN, f1_score, precision, recall = evamain(output_vcf_path, output_vcf)
    
    # 将文件名和评估结果保存到列表中
    evaluation_data.append([event_id_base, f1_score, precision, recall])

# 将评估结果写入 CSV 文件
write_evaluation_to_csv(csv_file_path, evaluation_data)

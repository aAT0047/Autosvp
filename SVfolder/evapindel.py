import os
import csv
import pysam
from evaf1 import *

import subprocess  # 用于运行外部命令
def get_vcf_files(input_dir):
    """从目录中获取所有 .vcf 文件"""
    vcf_files = [f for f in os.listdir(input_dir) if f.endswith('.vcf.gz')]
    return vcf_files

# 新增函数用于将评估结果写入 CSV 文件
def write_evaluation_to_csv(csv_file, evaluation_data):
    """将评估数据写入 CSV 文件"""
    with open(csv_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['File Name', 'f1_score', 'precision', 'recall'])  # CSV的表头
        for row in evaluation_data:
            writer.writerow(row)

import pysam
import subprocess

def process_vcf(input_vcf_path, output_vcf_path):
    # 打开压缩的 VCF 文件
    vcf_in = pysam.VariantFile(input_vcf_path)

    # 生成索引文件 (使用 tabix)
    try:
        subprocess.run(['tabix', '-p', 'vcf', input_vcf_path], check=True)
        print(f"索引文件已生成：{input_vcf_path}.tbi")
    except subprocess.CalledProcessError as e:
        print(f"Error creating index for {input_vcf_path}: {e}")

    # 输出到新的 VCF 文件
    with open(output_vcf_path, 'w') as vcf_out:
        # 写入原 VCF 文件的头部信息
        vcf_out.write(str(vcf_in.header))

        # 字典用于存储每个染色体上的重排变异对，键是染色体名，值是一个列表，存储成对变异
        variant_pairs = {}

        # 遍历每条记录，按染色体和变异位置进行分组，并过滤掉非重排变异
        for record in vcf_in:
            svtype = record.info.get('SVTYPE', None)

            # 只保留易位、倒位等类型的结构变异（根据 SVTYPE 过滤）
            if svtype in ['TRA', 'INV']:  # TRA (易位), INV (倒位)
                chrom = record.chrom
                start = record.pos
                end = record.info.get('END', start)  # 使用 END 字段确定变异的结束位置
                svlen = abs(record.info.get('SVLEN', 0))  # 计算变异长度
                
                # 存储成对变异信息
                if chrom not in variant_pairs:
                    variant_pairs[chrom] = []
                variant_pairs[chrom].append((record, start, end, svlen))

        # 保留变异长度最大的成对重排变异，仅保留一对
        largest_pairs = []

        for chrom, variants in variant_pairs.items():
            # 按变异的长度 (SVLEN) 进行排序
            variants.sort(key=lambda x: x[3], reverse=True)
            
            # 只保留最大的变异对（最多一对）
            if len(variants) >= 2:
                largest_pairs.append(variants[0])  # 第一条
                largest_pairs.append(variants[1])  # 第二条
                break  # 只保留一对变异，完成后退出循环

        # 将筛选出的最大一对重排变异写入新的 VCF 文件
        for variant, _, _, _ in largest_pairs:
            chrom_start = variant.chrom   # 染色体1
            pos_start = variant.pos     # 起始位置1
            chrom_end = variant.chrom  # 假设同一染色体的变异
            pos_end = variant.info.get('END', pos_start)  # 结束位置2
            sv_type = variant.info.get('SVTYPE', 'UNK')  # 变异类型
            sv_size = abs(variant.info.get('SVLEN', 0))  # 变异大小
            support_reads = 200  # 支持的读段数
            unique_id = "CTX_Event"  # 生成唯一事件ID

            # 生成 VCF 格式的 ALT 和 INFO 字段
            ref = 'N'
            alt_1 = f'[{chrom_end}:{pos_end}[N'
            alt_2 = f'N[{chrom_start}:{pos_start}]'
            qual = '.'
            filter_ = 'PASS'
            
            # 根据不同的 ALT，构建两个不同的 INFO
            info_1 = (f'SVTYPE={sv_type};STRANDS=--:503;EVENT={unique_id};MATEID=1;'
                      f'CIPOS=0,0;CIEND=0,0;CIPOS95=0,0;CIEND95=0,0;SU={support_reads};'
                      f'PE=0;SR=20;END={pos_end}')
            info_2 = (f'SVTYPE={sv_type};STRANDS=--:503;EVENT={unique_id};MATEID=2;'
                      f'CIPOS=0,0;CIEND=0,0;CIPOS95=0,0;CIEND95=0,0;SU={support_reads};'
                      f'PE=0;SR=20;END={pos_start};SECONDARY')
            format_ = 'GT:SU:PE:SR'

            # 写入两条记录到 VCF 文件
            vcf_out.write(f"{chrom_start}\t{pos_start}\t.\t{ref}\t{alt_1}\t{qual}\t{filter_}\t{info_1}\t{format_}\t.\n")
            vcf_out.write(f"{chrom_end}\t{pos_end}\t.\t{ref}\t{alt_2}\t{qual}\t{filter_}\t{info_2}\t{format_}\t.\n")

    # 关闭 VCF 文件
    vcf_in.close()


# 主函数执行流程
input_directory = "/data/home/std_12/ShiHe/results/ALK-RET-ROS1bam/pindel"  # 输入文件夹路径
output_directory = "/data/home/std_12/ShiHe/results//ALK-RET-ROS1bam/pindelvcf"  # 输出文件夹路径
output_truth = "/data/home/std_12/ShiHe//ALK-RET-ROS1vcf"  # 用于 evaf1 的真值文件

# 确保输出目录存在
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

vcf_files = get_vcf_files(input_directory)

# 准备CSV文件
csv_file_path = os.path.join(output_directory, "Pindel.csv")

# 初始化评估结果列表
evaluation_data = []

# 遍历所有生成的 VCF 文件，并对每个文件进行评估
for vcf_file in vcf_files:
    # 提取 VCF 文件名作为事件 ID 的基础部分
    event_id_base = os.path.splitext(os.path.basename(vcf_file))[0]
    event_id_base = os.path.splitext(os.path.basename(event_id_base))[0]
    # 定义输入和输出路径
    input_vcf_path = os.path.join(input_directory, vcf_file)
    modified_vcf_path = os.path.join(output_directory, f"{event_id_base}_modified.vcf")

    # 处理VCF文件：保留最大的一对变异并生成新的 VCF 文件
    process_vcf(input_vcf_path, modified_vcf_path)
    output_vcf_path = os.path.join(output_truth, f"{event_id_base}.vcf")
    # 使用 evamain 函数进行评估 (output_vcf_path 为真值, input_vcf_path 为测试文件)
    TP, FP, FN, f1_score, precision, recall = evamain(output_vcf_path ,modified_vcf_path)
    
    # 将文件名和评估结果保存到列表中
    evaluation_data.append([event_id_base, f1_score, precision, recall])

# 将评估结果写入 CSV 文件
write_evaluation_to_csv(csv_file_path, evaluation_data)

print(f"评估已完成，结果已保存至 {csv_file_path}")

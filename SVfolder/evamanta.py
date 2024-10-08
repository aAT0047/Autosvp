import os
import csv
from evaf1 import *

def get_vcf_files(input_dir):
    """从目录中获取所有 .vcf 文件"""
    vcf_files = [f for f in os.listdir(input_dir) if f.endswith('.vcf')]
    return vcf_files

# 新增函数用于将评估结果写入 CSV 文件
def write_evaluation_to_csv(csv_file, evaluation_data):
    """将评估数据写入 CSV 文件"""
    with open(csv_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['File Name', 'f1_score', 'precision', 'recall'])  # CSV的表头
        for row in evaluation_data:
            writer.writerow(row)
# 主函数执行流程
input_directory = "/data/home/std_12/Manta_Results/All_Samples_VCF/ALK-RET-ROS1bam"  # 输入文件夹路径
output_directory = "/data/home/std_12/Manta_Results"  # 输出文件夹路径
output_truth = "/data/home/std_12/ShiHe/ALK-RET-ROS1vcf"  # 用于 evaf1 的真值文件

# 确保输出目录存在
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

vcf_files = get_vcf_files(input_directory)

    


# 准备CSV文件
csv_file_path = os.path.join(output_directory, "ALK-RET-ROS1evaluation_results.csv")

# 初始化评估结果列表
evaluation_data = []

# 遍历所有生成的 VCF 文件，并对每个文件进行评估
for vcf_file in vcf_files:
    # 提取 VCF 文件名作为事件 ID 的基础部分
    event_id_base = os.path.splitext(os.path.basename(vcf_file))[0]
    # 解决breakdancer后缀问题
    # event_id_base = os.path.splitext(event_id_base)[0]
    output_vcf_path = os.path.join(output_truth, f"{event_id_base}.vcf")
    output_vcf = os.path.join(input_directory, vcf_file)
    # 使用 evamain 函数进行评估
    TP, FP, FN, f1_score, precision, recall = evamain(output_vcf_path, output_vcf)
    
    # 将文件名和评估结果保存到列表中
    evaluation_data.append([event_id_base, f1_score, precision, recall])

# 将评估结果写入 CSV 文件
write_evaluation_to_csv(csv_file_path, evaluation_data)

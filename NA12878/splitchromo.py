import os
import subprocess
import argparse
from tqdm import tqdm

def split_chromosome(chr_name):
    folder_name = chr_name
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    output_file = os.path.join(folder_name, f"{chr_name}.bam")
    cmd = f"samtools view -@ 32 -b {args.bam_file} {chr_name} > {output_file}"
    subprocess.run(cmd, shell=True)

# 使用argparse接受命令行参数
parser = argparse.ArgumentParser(description="Split BAM file by chromosomes.")
parser.add_argument("bam_file", help="Path to the BAM file.")
args = parser.parse_args()

# 首先为BAM文件创建索引
cmd_index = f"samtools index {args.bam_file}"
subprocess.run(cmd_index, shell=True)

# 获取所有染色体的列表
cmd = f"samtools idxstats {args.bam_file}"
result = subprocess.check_output(cmd, shell=True).decode('utf-8')
chromosomes = [line.split('\t')[0] for line in result.splitlines() if int(line.split('\t')[2]) > 0]

# 顺序处理每个染色体并显示进度条
for chr_name in tqdm(chromosomes):
    split_chromosome(chr_name)

print("BAM file has been split by chromosomes!")

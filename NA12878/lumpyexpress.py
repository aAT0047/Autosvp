#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
import subprocess
import re
import pandas as pd
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import csv
import sys
import time

# Define the global variables
output_dir = "/home/cloudam/simulat_2/normal/bam"
vcf_dir = "/home/cloudam/simulat_2/normal/vcf"
def execute_command(cmd):
    return subprocess.check_output(cmd, shell=True).decode('utf-8').strip()

def init_tqdm(pbar):
    def tqdm_callback(*a, **k):
        pbar.update()
    return tqdm_callback

def process_bam_files(sample_name):
    # output_dir = "/home/cloudam/simulat_2/sample_bam"
    
    # 提取discordant配对末端比对数据
    bam_file = os.path.join(output_dir, "{}.bam".format(sample_name))
    
    # 索引初始的BAM文件
    subprocess.call(["samtools", "index", bam_file])
    discordants_unsorted_file = os.path.join(output_dir, "{}.discordants.unsorted.bam".format(sample_name))
    cmd = "samtools view -b -F 1294 {} > {}".format(bam_file, discordants_unsorted_file)
    subprocess.call(cmd, shell=True)

    # Extract the split-read alignments
    splitters_unsorted_file = os.path.join(output_dir, "{}.splitters.unsorted.bam".format(sample_name))
    cmd = 'samtools view -h {} | /home/cloudam/.conda/envs/py2env/share/lumpy-sv-0.2.14a-2/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > {}'.format(bam_file, splitters_unsorted_file)
    subprocess.call(cmd, shell=True)

    # 排序discordant和split-read数据
    discordants_file = os.path.join(output_dir, "{}.discordants.bam".format(sample_name))
    subprocess.call(["samtools", "sort", "-o", discordants_file, discordants_unsorted_file])
    
    splitters_file = os.path.join(output_dir, "{}.splitters.bam".format(sample_name))
    subprocess.call(["samtools", "sort", "-o", splitters_file, splitters_unsorted_file])
    
    # 索引sorted BAM文件
    subprocess.call(["samtools", "index", discordants_file])
    subprocess.call(["samtools", "index", splitters_file])

    # 清理未排序的文件
    os.remove(discordants_unsorted_file)
    os.remove(splitters_unsorted_file)

def process_sample(data):
    output_vcf_path = os.path.join(vcf_dir, "{}.vcf".format(sample_name))
    unique_output = "/home/cloudam/simulat_2/normal/unique_output"
    # 创建日志文件名
    log_file = os.path.join(log_folder, "{}_log.txt".format(sample_name))

    # 打开日志文件并写入内容
    with open(log_file, 'w') as f:
        f.write("开始处理样本 {}".format(sample_name))
    # 初始化CSV文件
    csv_path = "/home/cloudam/simulat_2/normal/csv"
    if not os.path.exists(csv_path):
        os.mkdir(csv_path)
    csv_file = "{}/{}_params.csv".format(csv_path, sample_name)
    with open(csv_file, 'wb') as f:
        f.write("f1_score,precision,recall\n")
        # 创建新文件夹
    new_folder = os.path.join(output_dir, sample_name)
    if not os.path.exists(new_folder):
        os.makedirs(new_folder)
    # 切换到新文件夹
    os.chdir(new_folder)
    # 定义当前运行的结果文件名
    output_vcf = "{}/{}.vcf".format(output_dir, sample_name)
    # lumpyexpress -B 2300.bam -S 2300.splitters.bam -D 2300.discordants.bam -o output.vcf
    cmd = [
    "lumpyexpress",
    "-B", "{}/{}.bam".format(output_dir, sample_name),
    "-S", "{}/{}.splitters.bam".format(output_dir, sample_name),
    "-D", "{}/{}.discordants.bam".format(output_dir, sample_name),
    "-o", output_vcf
    ]
    subprocess.call(cmd)
    unique_output_filename = "output_f1_score_{}_iter.txt".format(sample_name)
    unique_output = os.path.join(log_folder, unique_output_filename)
    subprocess.call("evaf1.py {} {} {}".format(output_vcf_path, output_vcf, unique_output), shell=True)
                                        
    with open(unique_output, 'r') as f:
        lines = f.readlines()
        f1_score = lines[0].strip()
        precision = lines[1].strip()
        recall = lines[2].strip()
    # 指定生成的VCF文件和目标文件夹
    generated_vcf = output_vcf  # 这里应该是生成的VCF文件的完整路径
    target_folder = '/home/cloudam/simulat_2/normal/vcf/'

    # 使用mv命令移动VCF文件到目标文件夹
    mv_cmd = 'mv {} {}'.format(generated_vcf, target_folder)
    subprocess.call(mv_cmd, shell=True)
    with open(csv_file, 'a') as f:
        f.write("{},{},{}\n".format( f1_score, precision, recall))
    print("Sample {} processed.".format(sample_name))



if __name__ == "__main__":
    log_folder = "/home/cloudam/simulat_2/normal/log"
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)

    # sample_list = [1, 3, 4, 7, 8, 9, 12, 13, 14, 15, 16, 18]
    sample_list = [1, 3, 4]
    sample_names =[str(i) for i in sample_list]
    sample_datas =   sample_names 
   
    print('准备discordant和split-read数据')
    with tqdm(total=len(sample_datas), desc="Processing BAM files") as pbar:
        pool = Pool(processes=cpu_count(), initializer=init_tqdm, initargs=(pbar,))
        for _ in pool.imap_unordered(process_bam_files, sample_datas):
            pbar.update()  # 更新进度条
        pool.close()
        pool.join()
    print('discordant和split-read数据完成')
    print('开始运行LUMPY')
    # 使用所有可用的CPU核心，并添加进度条
    with tqdm(total=len(sample_names), desc="Processing samples") as pbar:
        pool = Pool(processes=cpu_count(), initializer=init_tqdm, initargs=(pbar,))
        for _ in pool.imap(process_sample, sample_names):
            pbar.update()  # 更新进度条
        pool.close()  # 关闭池，防止更多的任务被添加到池中
        pool.join()   # 等待所有任务完成
    print("All samples processed.")
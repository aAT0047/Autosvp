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
output_dir = "/home/cloudam/NA12878/sample_bam"
vcf_dir = "/home/cloudam/NA12878/sample_vcf"
# 从CSV文件读取前305个样本的数据


def read_and_clean_csv(filename):
    print('开始清理分布小于10的数据')
    
    # 读取csv文件到DataFrame
    df = pd.read_csv(filename)
    
    # 删除含有空值的行
    df_cleaned = df.dropna()

    # 将处理后的数据保存回csv文件
    df_cleaned.to_csv(filename, index=False)
    print('清理数据结束')
    
    # 将mean和stdev列转换为浮点数
    df_cleaned.loc[:, 'Mean'] = df_cleaned['Mean'].astype(float)
    df_cleaned.loc[:, 'Stdev'] = df_cleaned['Stdev'].astype(float)

   
    #df['Sample'] = df['Sample'].astype(str)

    samples = list(df_cleaned.itertuples(index=False, name=None))
    return samples

def execute_command(cmd):
    return subprocess.check_output(cmd, shell=True).decode('utf-8').strip()

def init_tqdm(pbar):
    def tqdm_callback(*a, **k):
        pbar.update()
    return tqdm_callback

def parse_args():
    parser = argparse.ArgumentParser(description="从CSV文件读取前305个样本的数据")
    parser.add_argument("csv_path", type=str, help="CSV文件的路径")
    return parser.parse_args()

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
    sample_name, mean, stdev = data
    sample_name =str(sample_name)
    # histo_path = os.path.join(output_dir, "{}.lib1.histo".format(sample_name))
    output_vcf_path = os.path.join(vcf_dir, "{}.vcf".format(sample_name))
    # 定义参数范围
    w_values = [10000]  # Modified this line to be a list
    msw_values = [2,3]
    tt_values = [0,10]
    back_distance_values = [10,20]
    min_mapping_threshold_values =[0,1] 
    min_clip_values = [10,20]
    read_length_values = [150,152]
    min_non_overlap_values = [148,185]
    discordant_z_values = [4,5,6]

    # 定义日志文件夹路径
    # log_folder = "/home/cloudam/NA12878/log"
    unique_output = "/home/cloudam/NA12878/unique_output"

    # 创建日志文件名
    log_file = os.path.join(log_folder, "{}_log.txt".format(sample_name))

    # 打开日志文件并写入内容
    with open(log_file, 'w') as f:
        f.write("开始处理样本 {}".format(sample_name))

    
    # 初始化CSV文件
    csv_path = "/home/cloudam/NA12878/csv"
    if not os.path.exists(csv_path):
        os.mkdir(csv_path)
    csv_file = "{}/{}_params.csv".format(csv_path, sample_name)
    # with open(csv_file, 'wb') as f:
    #     f.write("msw,tt,back_distance,min_mapping_threshold,min_clip,read_length,min_non_overlap,discordant_z,f1_score\n")
    with open(csv_file, 'wb') as f:
        f.write("msw,tt,back_distance,min_mapping_threshold,min_clip,read_length,min_non_overlap,discordant_z,f1_score,precision,recall\n")

    total_iterations = len(w_values) * len(msw_values) * len(tt_values) * len(back_distance_values) * len(min_mapping_threshold_values) * len(min_clip_values) * len(read_length_values) * len(min_non_overlap_values) * len(discordant_z_values)
    current_iteration = 0
    # 创建新文件夹
    new_folder = os.path.join(output_dir, sample_name)
    if not os.path.exists(new_folder):
        os.makedirs(new_folder)

    # 切换到新文件夹
    os.chdir(new_folder)
    # 循环运行LUMPY，并记录参数和VCF文件
    for w in w_values:
        for msw in msw_values:
            for tt in tt_values:
                for back_distance in back_distance_values:
                    for min_mapping_threshold in min_mapping_threshold_values:
                        for min_clip in min_clip_values:
                            for read_length in read_length_values:
                                for min_non_overlap in min_non_overlap_values:
                                    for discordant_z in discordant_z_values:

                                        # 定义当前运行的结果文件名
                                        output_vcf = "{}/{}_variants_w{}_msw{}_tt{}_bd{}_mmt{}_mc{}_mean{}_stdev{}_rl{}_mno{}_dz{}.vcf".format(output_dir, sample_name, w, msw, tt, back_distance, min_mapping_threshold, min_clip, mean, stdev, read_length, min_non_overlap, discordant_z)

                                        # 打印当前参数设置到日志文件
                                        with open(log_file, 'a') as f:
                                            f.write("正在运行参数设置：w={}, msw={}, tt={}, back_distance={}, min_mapping_threshold={}, min_clip={}, mean={}, stdev={}, read_length={}, min_non_overlap={}, discordant_z={}\n".format(w, msw, tt, back_distance, min_mapping_threshold, min_clip, mean, stdev, read_length, min_non_overlap, discordant_z))
                                            
                                        # 运行LUMPY并记录到VCF文件
                                        # 构建lumpy命令
                                        cmd = [
                                            "lumpy",
                                            "-w", str(w),
                                            "-msw", str(msw),
                                            "-tt", str(tt),
                                            "-pe", 'id:{},read_group:readgroup1,bam_file:{}/{},histo_file:{}/{},mean:{},stdev:{},read_length:{},min_non_overlap:{},discordant_z:{},back_distance:{},weight:1,min_mapping_threshold:{}'.format(sample_name, output_dir, sample_name + ".discordants.bam", output_dir, sample_name + ".lib1.histo", mean, stdev, read_length, min_non_overlap, discordant_z, back_distance, min_mapping_threshold),
                                            "-sr", 'id:{},bam_file:{}/{},back_distance:{},weight:1,min_mapping_threshold:{}'.format(sample_name, output_dir, sample_name + ".splitters.bam", back_distance, min_mapping_threshold),
                                           
                                        ]
                                        # if proc.returncode != 0:
                                        #     print("Error executing lumpy command.")
                                        # else:
                                        #     print("lumpy command executed successfully.")
                                        # 使用subprocess执行命令
                                        # Execute the command with Python-based redirection
                                        with open(output_vcf, 'w') as out_vcf, open(log_file_path, 'w') as logf:
                                            process = subprocess.Popen(cmd, stdout=out_vcf, stderr=logf)
                                            process.communicate()

                                        # unique_output_filename = "output_f1_score_{}.txt".format(sample_name)  # 使用样本名称作为文件名的一部分
                                        unique_output_filename = "output_f1_score_{}_iter{}.txt".format(sample_name, current_iteration)
                                        unique_output = os.path.join(log_folder, unique_output_filename)
                                        subprocess.call("evaf1.py {} {} {}".format(output_vcf_path, output_vcf, unique_output), shell=True)
                                        
                                        with open(unique_output, 'r') as f:
                                            lines = f.readlines()
                                            f1_score = lines[0].strip()
                                            precision = lines[1].strip()
                                            recall = lines[2].strip()

                                        # 指定生成的VCF文件和目标文件夹
                                        generated_vcf = output_vcf  # 这里应该是生成的VCF文件的完整路径
                                        target_folder = '/home/cloudam/NA12878/vcf/'

                                        # 使用mv命令移动VCF文件到目标文件夹
                                        mv_cmd = 'mv {} {}'.format(generated_vcf, target_folder)
                                        subprocess.call(mv_cmd, shell=True)
                                        # 创建或打开CSV文件以追加内容
                                        # with open(csv_file, 'a') as f:
                                        #     f.write("{},{},{},{},{},{},{},{},{}\n".format( msw, tt, back_distance, min_mapping_threshold, min_clip, read_length, min_non_overlap, discordant_z,  f1_score))
                                        # 创建或打开CSV文件以追加内容
                                        with open(csv_file, 'a') as f:
                                            f.write("{},{},{},{},{},{},{},{},{},{},{}\n".format(msw, tt, back_distance, min_mapping_threshold, min_clip, read_length, min_non_overlap, discordant_z, f1_score, precision, recall))

                                        # ...
                                        current_iteration += 1
                                        sys.stdout.write("\r进度: {}/{} ({:.2f}%)".format(current_iteration, total_iterations, (float(current_iteration)/total_iterations)*100))
                                        sys.stdout.flush()

    print("Sample {} processed.".format(sample_name))

if __name__ == "__main__":

    log_folder = "/home/cloudam/NA12878/log"
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)

    args = parse_args()
    if not os.path.exists("/home/cloudam/NA12878/vcf"):
        os.makedirs("/home/cloudam/NA12878/vcf")
    
    sample_names = read_and_clean_csv(args.csv_path)
    sample_datas = [row[0] for row in sample_names]  # 提取每一行的第一个元素，即样本名
   # sample_names =[str(i) for i in range(1, 4801)]
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

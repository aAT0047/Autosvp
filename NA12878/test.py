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

output_dir = "/home/cloudam/simulat_2/sample_bam"
vcf_dir = "/home/cloudam/simulat_2/sample_vcf"
output_dirt = "/home/cloudam/simulat_2/sample_bam_test"
def process_parameters(args):
    
    # 在这里，你可以对这些参数进行你想要的处理
    sample_name, msw, tt, back_distance, min_mapping_threshold, min_clip, read_length, min_non_overlap, discordant_z ,mean ,stdev= args
#    print sample_name, msw, tt, back_distance, min_mapping_threshold, min_clip, read_length, min_non_overlap, discordant_z ,mean ,stdev
    output_vcf_path = os.path.join(vcf_dir, "{}.vcf".format(sample_name))
    w_values = [10000]  # Modified this line to be a list
    msw_values =[ msw]
    tt_values =[tt]
    back_distance_values =[ back_distance]
    min_mapping_threshold_values =[ min_mapping_threshold]
    min_clip_values =[ min_clip]
    read_length_values =[ read_length]
    min_non_overlap_values =[ min_non_overlap]
    discordant_z_values = [discordant_z] 
   # mean = [mean]
   # stdev = [stdev]
    sample_name =str( sample_name  )
    # 初始化CSV文件
    csv_path = "/home/cloudam/simulat_2/csvtest"
    if not os.path.exists(csv_path):
        os.makedirs(csv_directory)
    csv_file = "{}/{}_params.csv".format(csv_path, sample_name)
    # with open(csv_file, 'wb') as f:
    #     f.write("msw,tt,back_distance,min_mapping_threshold,min_clip,read_length,min_non_overlap,discordant_z,f1_score\n")
    with open(csv_file, 'wb') as f:
        f.write("sample_name,f1_score,precision,recall\n")

    total_iterations = len(w_values) * len(msw_values) * len(tt_values) * len(back_distance_values) * len(min_mapping_threshold_values) * len(min_clip_values) * len(read_length_values) * len(min_non_overlap_values) * len(discordant_z_values)
    current_iteration = 0
    # 创建新文件夹
    new_folder = os.path.join(output_dirt, sample_name)
    if not os.path.exists(new_folder):
        os.makedirs(new_folder)
    unique_output = "/home/cloudam/simulat_2/unique_output"
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
                                        output_vcf = "{}/{}_variants_w{}_msw{}_tt{}_bd{}_mmt{}_mc{}_mean{}_stdev{}_rl{}_mno{}_dz{}.vcf".format(output_dirt, sample_name, w, msw, tt, back_distance, min_mapping_threshold, min_clip, mean, stdev, read_length, min_non_overlap, discordant_z)

                                        # 运行LUMPY并记录到VCF文件
                                        # 构建lumpy命令
                                        cmd = [
                                            "lumpy",
                                            "-w", str(w),
                                            "-msw", str(msw),
                                            "-tt", str(tt),
                                            "-pe", 'id:{},read_group:readgroup1,bam_file:{}/{},histo_file:{}/{},mean:{},stdev:{},read_length:{},min_non_overlap:{},discordant_z:{},back_distance:{},weight:1,min_mapping_threshold:{}'.format(sample_name, output_dir, sample_name + ".discordants.bam", output_dir, sample_name + ".lib1.histo", mean, stdev, read_length, min_non_overlap, discordant_z, back_distance, min_mapping_threshold),
                                            "-sr", 'id:{},bam_file:{}/{},back_distance:{},weight:1,min_mapping_threshold:{}'.format(sample_name, output_dir, sample_name + ".splitters.bam", back_distance, min_mapping_threshold),
                                            ">", output_vcf,
                                            "2>&1"
                                        ]
                       
                                        proc = subprocess.Popen(" ".join(cmd), shell=True)
                                        proc.communicate()

                                        # unique_output_filename = "output_f1_score_{}.txt".format(sample_name)  # 使用样本名称作为文件名的一部分
                                        unique_output_filename = "output_f1_score_{}_iter{}.txt".format(sample_name, current_iteration)
                                       # unique_output = os.path.join(log_folder, unique_output_filename)
                                        subprocess.call("evaf1.py {} {} {}".format(output_vcf_path, output_vcf, unique_output), shell=True)
                                        
                                        with open(unique_output, 'r') as f:
                                            lines = f.readlines()
                                            f1_score = lines[0].strip()
                                            precision = lines[1].strip()
                                            recall = lines[2].strip()

                                        # 指定生成的VCF文件和目标文件夹
                                        generated_vcf = output_vcf  # 这里应该是生成的VCF文件的完整路径
                                        target_folder = '/home/cloudam/simulat/vcf/'

                                        # 使用mv命令移动VCF文件到目标文件夹
                                        mv_cmd = 'mv {} {}'.format(generated_vcf, target_folder)
                                        subprocess.call(mv_cmd, shell=True)
                                        # 创建或打开CSV文件以追加内容
                                        # with open(csv_file, 'a') as f:
                                        #     f.write("{},{},{},{},{},{},{},{},{}\n".format( msw, tt, back_distance, min_mapping_threshold, min_clip, read_length, min_non_overlap, discordant_z,  f1_score))
                                        # 创建或打开CSV文件以追加内容
                                        with open(csv_file, 'a') as f:
                                            f.write("{},{},{},{}\n".format(sample_name,f1_score, precision, recall))

                                        # ...
                                        current_iteration += 1
                                        sys.stdout.write("\r进度: {}/{} ({:.2f}%)".format(current_iteration, total_iterations, (float(current_iteration)/total_iterations)*100))
                                        sys.stdout.flush()

#    print("Sample {} processed.".format(sample_name))

def parse_args():
    parser = argparse.ArgumentParser(description="从CSV文件读取样本的数据")
    parser.add_argument("csv_path", type=str, help="CSV文件的路径")
    parser.add_argument("output", type=str, help="输出文件名 (默认为 merged_lumpytestdata.csv)")
    return parser.parse_args()

def main():
    args = parse_args()

    tasks = []

    with open(args.csv_path, 'rb') as csvfile: 
        reader = csv.DictReader(csvfile)
        for row in reader:
            task = (row['sample_name'], row['msw'], row['tt'], row['back_distance'], row['min_mapping_threshold'], row['min_clip'], row['read_length'], row['min_non_overlap'], row['discordant_z'],row['mean'],row['stdev'])
            tasks.append(task)
   # print(tasks)
    # Create a pool of worker processes
    pool = Pool(processes=cpu_count())

    # Use tqdm with the pool's map function
    for _ in tqdm(pool.imap_unordered(process_parameters, tasks), total=len(tasks)):
        pass

    pool.close()
    pool.join()



if __name__ == '__main__':
    main()
    print('merged into test_data start')
    # Specify the path
    path = "/home/cloudam/simulat_2/csvtest"
    paths = "/home/cloudam/simulat_2"
    # Generate the list of filenames with the full path
    file_names = os.listdir(path)
    # 按照字母顺序排序
    sorted_file_names = sorted(file_names)
    # 提取文件名中的数字部分并保存到列表中
    numbers = [name.split("_")[0] for name in sorted_file_names]
    files = ["{}/{}_params.csv".format(path, i) for i in numbers]

    # Check if the first file exists, and if so, read its content
    if os.path.exists(files[0]):
        merged_df = pd.read_csv(files[0])
    else:
        raise ValueError("File {} not found.".format(files[0]))

    # From the second file onwards, only read the content and skip the header
    for file in tqdm(files[1:], desc="Merging files", unit="file"):
        if os.path.exists(file):
            df = pd.read_csv(file)
            merged_df = pd.concat([merged_df, df], ignore_index=True)
            # 删除原始的CSV文件
            os.remove(file)
    # Save the combined data to a new CSV file
    args = parse_args()
    output_file = args.output if args.output else "merged_lumpytestdata.csv"
    # 如果没有提供 --output 参数，将使用默认文件名 "merged_lumpytestdata.csv"
    # 保存到CSV文件
    merged_df.to_csv(output_file, index=False)
        # 构建匹配模式来查找所有带有"_optparams"的CSV文件
    file_pattern = os.path.join(path, '_params.csv')
    import glob
    # 查找匹配的文件
    matching_files = glob.glob(file_pattern)

    # 遍历并删除匹配的文件
    for file_path in matching_files:
        os.remove(file_path)

    print("All CSV files merged into test_data.csv")

#!/usr/bin/env python
# -*- coding: utf-8 -*-
from multiprocessing import Pool, cpu_count
import os
import subprocess

def execute_command(cmd):
    return subprocess.check_output(cmd, shell=True).strip()

output_dir = "/home/cloudam/simulat_2/sample_bam"
vcf_dir = "/home/cloudam/simulat_2/sample_vcf"
# 确保文件夹存在
if not os.path.exists("/home/cloudam/simulat_2/vcf"):
    os.makedirs("/home/cloudam/simulat_2/vcf")



# 循环处理每个BAM文件
for i in range(116, 117):  # Assuming you have 1.bam to 4800.bam
    sample_name = str(i)
    output_vcf_path = os.path.join(vcf_dir, "{}.vcf".format(sample_name))
    
    # # 提取discordant配对末端比对数据
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
    histo_path = os.path.join(output_dir, "{}.lib1.histo".format(sample_name))


    # 首先，生成 BAM 文件中每个库的经验插入大小统计数据
    cmd1 = 'samtools view {} | tail -n+10 | /home/cloudam/.conda/envs/py2env/share/lumpy-sv-0.2.14a-2/scripts/pairend_distro.py -r150  -X 4 -N 100 -o make.lib1.histo'.format(bam_file)
    output = execute_command(cmd1)
    mean = execute_command('echo "{}" | grep -oP "mean:\\K[^,\t]+" | head -n 1'.format(output))
    stdev = execute_command('echo "{}" | grep -oP "stdev:\\K[^,\t]+" | head -n 1'.format(output))



    # 定义参数范围
    w_values = [1000000]  # Modified this line to be a list
    msw_values = xrange(1,5, 1)
    tt_values = xrange(0, 100, 10)
    back_distance_values = xrange(10, 100000, 100)
    min_mapping_threshold_values = xrange(1,100, 200)
    min_clip_values = xrange(0, 1000, 10)
    read_length_values = [150]
    min_non_overlap_values = xrange(148, 155, 4)
    discordant_z_values = xrange(2, 20, 5)

    # 创建日志文件
    log_file = "{}_log.txt".format(sample_name)
    with open(log_file, 'wb') as f:
        f.write("开始处理样本 {}".format(sample_name))

    # 初始化CSV文件
    csv_path = "/home/cloudam/simulat_2/csvpro"
    if not os.path.exists(csv_path):
        os.mkdir(csv_path)
    csv_file = "{}/{}_params.csv".format(csv_path, sample_name)
    with open(csv_file, 'wb') as f:
        f.write("w,msw,tt,back_distance,min_mapping_threshold,min_clip,read_length,min_non_overlap,discordant_z,TP,FP,FN,Precision,Recall,f1_score\n")

    total_iterations = len(w_values) * len(msw_values) * len(tt_values) * len(back_distance_values) * len(min_mapping_threshold_values) * len(min_clip_values) * len(read_length_values) * len(min_non_overlap_values) * len(discordant_z_values)
    current_iteration = 0

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
                                        with open(log_file, 'ab') as f:
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
                                            ">", output_vcf,
                                            "2>&1"
                                        ]
                                        # if proc.returncode != 0:
                                        #     print("Error executing lumpy command.")
                                        # else:
                                        #     print("lumpy command executed successfully.")
                                        # 使用subprocess执行命令
                                        proc = subprocess.Popen(" ".join(cmd), shell=True)
                                        proc.communicate()

                                        subprocess.call("evaf1.py {} {} output_f1_score.txt".format(output_vcf_path, output_vcf), shell=True)
                                        with open('output_f1_score.txt', 'r') as f:
    
                                            lines = f.readlines()
                                            TP = lines[0].strip()
                                            FP = lines[1].strip()
                                            FN = lines[2].strip()
                                            f1_score = lines[3].strip()
                                            precision = lines[4].strip()
                                            recall = lines[5].strip()
                                        # 指定生成的VCF文件和目标文件夹
                                        generated_vcf = output_vcf  # 这里应该是生成的VCF文件的完整路径
                                        target_folder = '/home/cloudam/simulat_2/vcf/'

                                        # 使用mv命令移动VCF文件到目标文件夹
                                        mv_cmd = 'mv {} {}'.format(generated_vcf, target_folder)
                                        subprocess.call(mv_cmd, shell=True)
                                        # 将参数配置保存到CSV文件
                                        # 假设所有的变量（如w, msw, tt等）已经在Python代码中定义了

                                        # 创建或打开CSV文件以追加内容
                                        with open(csv_file, 'a') as f:
                                            f.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(w, msw, tt, back_distance, min_mapping_threshold, min_clip, read_length, min_non_overlap, discordant_z, TP,FP,FN,Precision,Recall,f1_score))

                                        current_iteration += 1
                                        print current_iteration

    print "\n所有参数设置的运行已完成。"
    with open(log_file, 'ab') as f:
        f.write("所有参数设置的运行已完成。\n")



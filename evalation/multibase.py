# -*- coding: utf-8 -*-
import subprocess
import os
import itertools
import multiprocessing

# 参数验证
if len(os.sys.argv) != 3:
    print ("用法: script.py <样本名称> <output.vcf路径>")
    os.sys.exit(1)

# 提取命令行参数中的样本名称
sample_name = os.sys.argv[1]
output_vcf_path = os.sys.argv[2]

# 创建输出目录
output_dir = "{}_output".format(sample_name)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 定义文件名参数
reference_genome = "/home/cloudam/my_folder_graph/ref.fasta"
input_fq1 = "{}1.fq".format(sample_name)
input_fq2 = "{}2.fq".format(sample_name)

# 获取mean和stdev
command = "samtools view {}/{}.bam | tail -n+10 | /home/cloudam/.conda/envs/py2env/share/lumpy-sv-0.2.14a-2/scripts/pairend_distro.py -r 75 -X 4 -N 100000 -o make.lib1.histo".format(output_dir, sample_name)
output = subprocess.check_output(command, shell=True)
mean = output.split('mean:')[1].split(',')[0]
stdev = output.split('stdev:')[1].split(',')[0]

# 定义参数范围
w_values = [100, 5000, 10000]
msw_values = [2, 3, 4, 5]
tt_values = range(0, 11, 5)
back_distance_values = range(10, 51, 10)
min_mapping_threshold_values = range(0, 3, 1)
min_clip_values = range(0, 31, 10)
read_length_values = range(70, 75, 2)
min_non_overlap_values = range(70, 75, 4)
discordant_z_values = range(2, 7, 2)

# 创建日志文件和CSV文件
with open("{}_log.txt".format(sample_name), 'w') as log_file, open("{}_params.csv".format(sample_name), 'w') as csv_file:
    log_file.write("开始处理样本 {}\n".format(sample_name))
    csv_file.write("w,msw,tt,back_distance,min_mapping_threshold,min_clip,read_length,min_non_overlap,discordant_z,output_vcf,f1_score\n")

    # 生成所有可能的参数组合
    param_combinations = list(itertools.product(w_values, msw_values, tt_values, back_distance_values, min_mapping_threshold_values, min_clip_values, read_length_values, min_non_overlap_values, discordant_z_values))

    # 使用多进程并行处理参数组合
    num_processes = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)

    def process_param_combination(params):
        w, msw, tt, back_distance, min_mapping_threshold, min_clip, read_length, min_non_overlap, discordant_z = params
        output_vcf = "{}/{}_variants_w{}_msw{}_tt{}_bd{}_mmt{}_mc{}_mean{}_stdev{}_rl{}_mno{}_dz{}.vcf".format(output_dir, sample_name, w, msw, tt, back_distance, min_mapping_threshold, min_clip, mean, stdev, read_length, min_non_overlap, discordant_z)

        # 记录到日志文件
        log_file.write("正在运行参数设置:w={}, msw={}, tt={}, back_distance={}, min_mapping_threshold={}, min_clip={}, mean={}, stdev={}, read_length={}, min_non_overlap={}, discordant_z={}\n".format(w, msw, tt, back_distance, min_mapping_threshold, min_clip, mean, stdev, read_length, min_non_overlap, discordant_z))

        # 运行LUMPY并记录到VCF文件
        lumpy_command = """lumpy -w {} -msw {} -tt {} \
        -pe id:{},read_group:readgroup1,bam_file:{}/{}.discordants.bam,histo_file:{}/{}.lib1.histo,mean:{},stdev:{},read_length:{},min_non_overlap:{},discordant_z:{},back_distance:{},weight:1,min_mapping_threshold:{} \
        -sr id:{},bam_file:{}/{}.splitters.bam,back_distance:{},weight:1,min_mapping_threshold:{} \
        > {} 2>&1""".format(w, msw, tt, sample_name, output_dir, sample_name, output_dir, sample_name, mean, stdev, read_length, min_non_overlap, discordant_z, back_distance, min_mapping_threshold, sample_name, output_dir, sample_name, back_distance, min_mapping_threshold, output_vcf)

        import subprocess
        subprocess.call(lumpy_command, shell=True)


        # 计算F1分数
        f1_command = "evaf1.py {} {} output_f1_score.txt".format(output_vcf_path, output_vcf)
        subprocess.call(f1_command, shell=True)
        with open("output_f1_score.txt", 'r') as f1_file:
            f1_score = f1_file.read().strip()

        # 将参数配置保存到CSV文件
        csv_file.write("{},{},{},{},{},{},{},{},{},{},{},{}\n".format(w, msw, tt, back_distance, min_mapping_threshold, min_clip, mean, stdev, read_length, min_non_overlap, discordant_z, output_vcf, f1_score))

    pool.map(process_param_combination, param_combinations)
    pool.close()
    pool.join()

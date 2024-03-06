import subprocess
import os
import itertools

# 参数验证
if len(os.sys.argv) != 3:
    print("用法: script.py <样本名称> <output.vcf路径>")
    os.sys.exit(1)

# 提取命令行参数中的样本名称
sample_name = os.sys.argv[1]
output_vcf_path = os.sys.argv[2]

# 创建输出目录
output_dir = f"{sample_name}_output"
os.makedirs(output_dir, exist_ok=True)

# 定义文件名参数
reference_genome = "/home/cloudam/my_folder_graph/ref.fasta"
input_fq1 = f"{sample_name}1.fq"
input_fq2 = f"{sample_name}2.fq"

# 获取mean和stdev
command = f"samtools view {output_dir}/{sample_name}.bam | tail -n+10 | /home/cloudam/.conda/envs/py2env/share/lumpy-sv-0.2.14a-2/scripts/pairend_distro.py -r 75 -X 4 -N 100000 -o make.lib1.histo"
output = subprocess.getoutput(command)
mean = output.split('mean:')[1].split(',')[0]
stdev = output.split('stdev:')[1].split(',')[0]

# 定义参数范围
w_values = [100, 5000, 10000]
msw_values = [2, 3, 4, 5]
tt_values = list(range(0, 11, 5))
back_distance_values = list(range(10, 51, 10))
min_mapping_threshold_values = list(range(0, 31, 10))
min_clip_values = list(range(0, 31, 10))
read_length_values = list(range(70, 75, 2))
min_non_overlap_values = list(range(70, 75, 4))
discordant_z_values = list(range(2, 7, 2))

# 创建日志文件和CSV文件
with open(f"{sample_name}_log.txt", 'w') as log_file, open(f"{sample_name}_params.csv", 'w') as csv_file:
    log_file.write(f"开始处理样本 {sample_name}\n")
    csv_file.write("w,msw,tt,back_distance,min_mapping_threshold,min_clip,read_length,min_non_overlap,discordant_z,output_vcf,f1_score\n")

    # 生成所有可能的参数组合
    total_iterations = len(w_values) * len(msw_values) * len(tt_values) * len(back_distance_values) * len(min_mapping_threshold_values) * len(min_clip_values) * len(read_length_values) * len(min_non_overlap_values) * len(discordant_z_values)
    current_iteration = 0

    for w, msw, tt, back_distance, min_mapping_threshold, min_clip, read_length, min_non_overlap, discordant_z in itertools.product(w_values, msw_values, tt_values, back_distance_values, min_mapping_threshold_values, min_clip_values, read_length_values, min_non_overlap_values, discordant_z_values):
        output_vcf = f"{output_dir}/{sample_name}_variants_w{w}_msw{msw}_tt{tt}_bd{back_distance}_mmt{min_mapping_threshold}_mc{min_clip}_mean{mean}_stdev{stdev}_rl{read_length}_mno{min_non_overlap}_dz{discordant_z}.vcf"

        # 记录到日志文件
        log_file.write(f"正在运行参数设置：w={w}, msw={msw}, tt={tt}, back_distance={back_distance}, min_mapping_threshold={min_mapping_threshold}, min_clip={min_clip}, mean={mean}, stdev={stdev}, read_length={read_length}, min_non_overlap={min_non_overlap}, discordant_z={discordant_z}\n")

        # 运行LUMPY并记录到VCF文件
        lumpy_command = f"""lumpy -w {w} -msw {msw} -tt {tt} \
        -pe id:{sample_name},read_group:readgroup1,bam_file:{output_dir}/{sample_name}.discordants.bam,histo_file:{output_dir}/{sample_name}.lib1.histo,mean:{mean},stdev:{stdev},read_length:{read_length},min_non_overlap:{min_non_overlap},discordant_z:{discordant_z},back_distance:{back_distance},weight:1,min_mapping_threshold:{min_mapping_threshold} \
        -sr id:{sample_name},bam_file:{output_dir}/{sample_name}.splitters.bam,back_distance:{back_distance},weight:1,min_mapping_threshold:{min_mapping_threshold} \
        > {output_vcf} 2>&1"""
        subprocess.run(lumpy_command, shell=True)

        # 计算F1分数
        f1_command = f"evaf1.py {output_vcf_path} {output_vcf} output_f1_score.txt"
        subprocess.run(f1_command, shell=True)
        with open("output_f1_score.txt", 'r') as f1_file:
            f1_score = f1_file.read().strip()

        # 将参数配置保存到CSV文件
        csv_file.write(f"{w},{msw},{tt},{back_distance},{min_mapping_threshold},{min_clip},{mean},{stdev},{read_length},{min_non_overlap},{discordant_z},{output_vcf},{f1_score}\n")
        
        current_iteration += 1
        print(f"Current iteration: {current_iteration}")
        percentage = int((current_iteration / total_iterations) * 100)
        print(f"Progress: [{'#' * (percentage // 2):<50}] {percentage}%")

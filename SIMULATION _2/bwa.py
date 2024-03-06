import subprocess
from tqdm import tqdm

# 您的参考基因组路径
REFERENCE = "/home/cloudam/my_folder_graph/ref.fasta"

# 设定您的数据文件夹的基路径
BASE_PATH = "/home/cloudam/simulat"

def run_command(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print("\nError: '{}'".format(err))  # 添加换行以避免与进度条冲突
        exit(p.returncode)
    return out

# 创建参考基因组的索引
print("正在为参考基因组创建索引...")
run_command('bwa index "{}"'.format(REFERENCE))
print("索引创建完成！")

# 遍历所有文件夹，并为每对文件执行比对
for FOLDER_NUM in tqdm(range(1, 49), desc="处理进度", ncols=100):
    FOLDER = "{}/my_folder_{}".format(BASE_PATH, FOLDER_NUM)
    R1 = "{}/make1.fq".format(FOLDER)
    R2 = "{}/make2.fq".format(FOLDER)
    OUT_BAM = "{}/make.bam".format(FOLDER)
    SAMPLE_NAME = "Sample_{}".format(FOLDER_NUM)

    cmd = ('bwa mem -R "@RG\tID:id\tSM:{0}\tLB:lib" "{1}" "{2}" "{3}" | '
           'samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 | '
           'samtools view -S -b - > "{4}"').format(SAMPLE_NAME, REFERENCE, R1, R2, OUT_BAM)

    run_command(cmd)

print("\n处理完成")

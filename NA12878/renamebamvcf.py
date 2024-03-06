import os
import shutil
from tqdm import tqdm

# 创建目标文件夹
os.makedirs("/home/cloudam/NA12878/sample_vcf", exist_ok=True)
os.makedirs("/home/cloudam/NA12878/sample_bam", exist_ok=True)

# 定义全局计数器
count = 1201

# 初始化tqdm进度条
pbar = tqdm(total=len(list2 )*200)  # 48文件夹，每个文件夹有200个文件 (100 bam + 100 vcf)
list2 = [str(i) for i in range(13, 23)] + ['X']
# 循环遍历每个文件夹中的 BAM 和 VCF 文件
for folder_num in list2 :
    # 在这里执行您的操作，使用 folder_num 变量来引用当前的文件夹编号

    src_folder = f"/home/cloudam/NA12878/{folder_num}/split_files"

    # 复制并重命名 BAM 文件

    for i in range(1,101 ):
        src_bam = f"{src_folder}/split_{i}.bam"
        dest_bam = f"/home/cloudam/NA12878/sample_bam/{count}.bam"
        if os.path.exists(src_bam):
            shutil.copy(src_bam, dest_bam)
            pbar.update(1)
            pbar.set_description(f"Copying {src_bam} to {dest_bam}")
            #count += 1  # 在成功拷贝后增加计数器
        else:
            pbar.write(f"Warning: {src_bam} does not exist!")
            pbar.update(1)

    # 复制并重命名 VCF 文
        src_vcf = f"{src_folder}/split_{i}.vcf"
        dest_vcf = f"/home/cloudam/NA12878/sample_vcf/{count}.vcf"
        if os.path.exists(src_vcf):
            shutil.copy(src_vcf, dest_vcf)
            pbar.update(1)
            pbar.set_description(f"Copying {src_vcf} to {dest_vcf}")
          
        else:
            pbar.write(f"Warning: {src_vcf} does not exist!")
            pbar.update(1)
        count +=1
pbar.close()

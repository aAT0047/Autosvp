import os
import subprocess
from tqdm import tqdm

BASE_PATH = "/home/cloudam/simulat"

for i in tqdm(range(1, 49)):
    # 设置文件夹路径
    folder_path = os.path.join(BASE_PATH, "my_folder_{}".format(i))
    
    # 设置.sh文件的路径
    sh_file = os.path.join("/home/cloudam/simulat/sh_files", "A_stableCallerPaperSimFlowShell_{}.sh".format(i))

    # 在该目录下执行shinvcf.py命令来生成base.vcf文件
    cmd = "shinvcf.py {} base.vcf".format(sh_file)
    subprocess.call(cmd, cwd=folder_path, shell=True)

print ("处理完成！")

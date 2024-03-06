import os
import subprocess
import argparse
from tqdm import tqdm

def main(args):
    vcf_file = args.vcf_file
    sorted_vcf_file = "NA12878_sorted.vcf.gz"
    chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']

    # 排序VCF文件
    cmd_sort = f"bcftools sort {vcf_file} -O z -o {sorted_vcf_file} "
    subprocess.check_call(cmd_sort, shell=True)

    # 为排序后的VCF创建索引
    cmd_index = f"bcftools index {sorted_vcf_file}"
    subprocess.check_call(cmd_index, shell=True)

    # 按染色体拆分并重命名
    for chr_name in tqdm(chromosomes, desc="Processing chromosomes"):
        folder_name = chr_name
        output_vcf = os.path.join(folder_name, "base.vcf")
        cmd_split = f"bcftools view -r {chr_name} {sorted_vcf_file} -O z -o {output_vcf}.gz "
        subprocess.check_call(cmd_split, shell=True)

        # 创建索引文件（如果您需要的话）
        cmd_index_per_chr = f"bcftools index {output_vcf}.gz"
        subprocess.check_call(cmd_index_per_chr, shell=True)

    print("VCF file has been split by chromosomes!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split VCF file by chromosomes.")
    parser.add_argument("vcf_file", help="Path to the VCF file.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads for bcftools.")
    args = parser.parse_args()
    main(args)

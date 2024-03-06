import os
import subprocess
from multiprocessing import Pool

# 1. 定义 execute_command 函数
def execute_command(cmd):
    """Execute the command and handle errors."""
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError:
        print(f"Error executing command: {cmd}")

def get_vcf_positions(vcf_file):
    cmd = ['bcftools', 'query', '-f', '%POS\n', vcf_file]
    try:
        output = subprocess.check_output(cmd).decode('utf-8')
        return [int(pos) for pos in output.splitlines()]
    except subprocess.CalledProcessError:
        print("Error executing bcftools command.")
        return []

def split_segment(args):
    bam_file, vcf_file, work_dir, chrom, start, end, i = args
    bam_output = os.path.join(work_dir, "split_files", f"split_{i + 1}.bam")
    vcf_output = os.path.join(work_dir, "split_files", f"split_{i + 1}.vcf")
    region = f"{chrom}:{start}-{end}"
    regionbam = f"{chrom[-1]}:{start}-{end}"
    cmd_bam = ['samtools', 'view', '-b', bam_file, regionbam, '-o', bam_output]
    subprocess.check_call(cmd_bam)

    cmd_vcf = ['bcftools', 'view', vcf_file, region, '-o', vcf_output]
    subprocess.check_call(cmd_vcf)

def split_bam_and_vcf(bam_file, vcf_file, num_splits, work_dir, chrom):
    positions = get_vcf_positions(vcf_file)
    if not positions:
        print(f"No variants found in VCF.")
        return

    variants_per_split = len(positions) // num_splits
    split_points = []
    for i in range(num_splits):
        start = positions[i * variants_per_split]
        if i == num_splits - 1:
            end = positions[-1]
        else:
            end = positions[(i + 1) * variants_per_split - 1]
        split_points.append((start, end))

    split_dir = os.path.join(work_dir, "split_files")
    os.makedirs(split_dir, exist_ok=True)

    with Pool(32) as pool:
        pool.map(split_segment, [(bam_file, vcf_file, work_dir, chrom, start, end, i) for i, (start, end) in enumerate(split_points)])

    return split_points

def process_folder(folder_num):
    folder_path = os.path.join(BASE_DIR, folder_num)

    # Process make.bam
    bam_file = os.path.join(folder_path, f"{folder_num}.bam")
    sorted_bam_file = os.path.join(folder_path, "make_sorted.bam")
    execute_command(f"samtools sort -@ 32 {bam_file} -o {sorted_bam_file}")
    execute_command(f"samtools index {sorted_bam_file}")

    # Process base.vcf
    # vcf_file = os.path.join(folder_path, "base.vcf")
    gz_vcf_file = os.path.join(folder_path, "base.vcf.gz")
    sorted_gz_vcf_file = os.path.join(folder_path, "base_sorted.vcf.gz")

    # Compress the VCF using bgzip
    # execute_command(f"bgzip -c {vcf_file} > {gz_vcf_file}")

    # Sort the compressed VCF using bcftools
    execute_command(f"bcftools sort {gz_vcf_file} -Oz -o {sorted_gz_vcf_file}")
    
    # Index the sorted compressed VCF using tabix
    execute_command(f"tabix -p vcf {sorted_gz_vcf_file}")

list2 = [str(i) for i in range(1, 23)] + ['X']

BASE_DIR = "/home/cloudam/NA12878"
for item in list2:
    # Sort, index, and compress bam and vcf
    process_folder(item)

    work_dir = os.path.join(BASE_DIR, item)
    bam_file = os.path.join(work_dir, "make_sorted.bam")
    vcf_file = os.path.join(work_dir, "base_sorted.vcf.gz")
    chrom = "chr" + item

    split_points = split_bam_and_vcf(bam_file, vcf_file,100 , work_dir, chrom)
   # print(split_points)
    if split_points:
        # Save split points to the new directory
        with open(os.path.join(work_dir, "split_files", "split_points.txt"), "w") as f:
            for start, end in split_points:
                f.write(f"{start}-{end}\n")

    print(f"Finished processing and splitting BAM and VCF for {work_dir}!")
# 这是一个完整的处理流程。请确保您的环境中已经安装了所有的命令行工具（例如samtools, bcftools, bgzip, 和 tabix）。此外，您的原始代码中的split_bam_and_vcf函数需要完整地被插入到上述代码中的相应位置。








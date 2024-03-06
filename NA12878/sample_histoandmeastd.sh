# 指定BAM文件和输入/输出目录
input_dir="/home/cloudam/NA12878/sample_bam"
output_dir="/home/cloudam/NA12878/sample"
# 创建输出目录
mkdir -p "${output_dir}"
# 创建CSV文件并写入标题行
csv_file="/home/cloudam/NA12878/meanstdev.csv"
echo "Sample,Mean,Stdev" > "${csv_file}"

# 循环处理每个样本
for sample_name in {1201..2300}; do
    # 生成经验插入大小统计数据
    samtools view "${input_dir}/${sample_name}.bam" \
        | tail -n+10 \
        | /home/cloudam/.conda/envs/py2env/share/lumpy-sv-0.2.14a-2/scripts/pairend_distro.py \
        -r 150 \
        -X 4 \
        -N 1000 \
        -o "${output_dir}/${sample_name}.lib1.histo"

    # 提取mean和stdev并写入CSV文件
    output=$(samtools view "${input_dir}/${sample_name}.bam" | tail -n+10 | /home/cloudam/.conda/envs/py2env/share/lumpy-sv-0.2.14a-2/scripts/pairend_distro.py -r 150 -X 4 -N 1000 -o "${output_dir}/${sample_name}.lib1.histo")
    mean=$(echo "$output" | grep -oP 'mean:\K[^,\t]+' | head -n 1)
    stdev=$(echo "$output" | grep -oP 'stdev:\K[^,\t]+' | head -n 1)
    # 移动文件回到input_dir
    mv "${output_dir}/${sample_name}.lib1.histo" "${input_dir}/${sample_name}.lib1.histo"
    # 写入CSV文件
    echo "${sample_name},${mean},${stdev}" >> "${csv_file}"

    echo "Sample ${sample_name} processed."
done

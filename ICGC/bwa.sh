#/bin/bash

# 您的参考基因组路径
REFERENCE="/home/cloudam/my_folder_graph/ref.fasta"

# 设定您的数据文件夹的基路径
BASE_PATH="/home/cloudam/ICGCsimulat_2"

# 创建参考基因组的索引
bwa index "$REFERENCE"

# 遍历所有文件夹，并为每对文件执行比对
for FOLDER_NUM in {1,18}; do
    FOLDER="$BASE_PATH/my_folder_$FOLDER_NUM"
    R1="$FOLDER/make1.fq"
    R2="$FOLDER/make2.fq"
    OUT_BAM="$FOLDER/make.bam"
    SAMPLE_NAME="Sample_$FOLDER_NUM"

    # 显示当前的进度
    echo "正在处理：$FOLDER_NUM / 2"

    # 使用 BWA MEM 进行比对，并使用 samblaster 标记 PCR 重复，最后转换为 BAM 格式
    bwa mem -t 32 -R "@RG\tID:id\tSM:${SAMPLE_NAME}\tLB:lib" "$REFERENCE" "$R1" "$R2" \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -S -b -@ 32 - \
    > "$OUT_BAM"
done

echo "处理完成！"

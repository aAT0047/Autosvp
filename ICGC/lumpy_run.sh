#!/bin/bash

# 参数验证
if [ "$#" -ne 2 ]; then
    echo "用法: $0 <样本名称> <output.vcf路径>"
    exit 1
fi

# 提取命令行参数中的样本名称
sample_name="$1"
output_vcf_path="$2"


# 创建输出目录
sudo mkdir -p /home/cloudam/ICGCsimulat_2/split_files


# 定义文件名参数

output_dir="/home/cloudam/ICGCsimulat_2/split_files"

# 提取discordant配对末端比对数据
samtools view -b -F 1294 "${output_dir}/${sample_name}.bam" > "${output_dir}/${sample_name}.discordants.unsorted.bam"

# Extract the split-read alignments
samtools view -h "${output_dir}/${sample_name}.bam" \
    | /home/cloudam/.conda/envs/py2env/share/lumpy-sv-0.2.14a-2/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > "${output_dir}/${sample_name}.splitters.unsorted.bam"


# 对discordant和split-read比对数据进行排序
samtools sort -o "${output_dir}/${sample_name}.discordants.bam" "${output_dir}/${sample_name}.discordants.unsorted.bam"
samtools sort -o "${output_dir}/${sample_name}.splitters.bam" "${output_dir}/${sample_name}.splitters.unsorted.bam"


# 清理未排序的文件（可选）
rm "${output_dir}/${sample_name}.discordants.unsorted.bam"
rm "${output_dir}/${sample_name}.splitters.unsorted.bam"

#首先，生成 BAM 文件中每个库的经验插入大小统计数据
samtools view  ${output_dir}/${sample_name}.bam \
    | tail -n+10 \
    | /home/cloudam/.conda/envs/py2env/share/lumpy-sv-0.2.14a-2/scripts/pairend_distro.py\
    -r 150\
    -X 4 \
    -N 1000 \
    -o $output_dir/$sample_name.lib1.histo


output=$(samtools view $output_dir/$sample_name.bam | tail -n+10 | /home/cloudam/.conda/envs/py2env/share/lumpy-sv-0.2.14a-2/scripts/pairend_distro.py -r 75 -X 4 -N 100000 -o make.lib1.histo)
mean=$(echo "$output" | grep -oP 'mean:\K[^,\t]+' | head -n 1)
stdev=$(echo "$output" | grep -oP 'stdev:\K[^,\t]+' | head -n 1)

# 定义参数范围
# 定义参数范围
w_values=(  10000000 )
msw_values=(2)
tt_values=( 0  )
back_distance_values=($(seq 20 200 1000))
min_mapping_threshold_values=(10)
min_clip_values=(20)
read_length_values=(  150 )
min_non_overlap_values=( 150  )
discordant_z_values=( 4  )


# 创建日志文件
log_file="${sample_name}_log.txt"
echo "开始处理样本 $sample_name" > "$log_file"

# 初始化CSV文件

csv_file="${sample_name}_params.csv"
echo "w,msw,tt,back_distance,min_mapping_threshold,min_clip,read_length,min_non_overlap,discordant_z,mean,stdev,f1_score,precision,recall" > "$csv_file"


total_iterations=$((${#w_values[@]} * ${#msw_values[@]} * ${#tt_values[@]} * ${#back_distance_values[@]} * ${#min_mapping_threshold_values[@]} * ${#min_clip_values[@]} * ${#read_length_values[@]} * ${#min_non_overlap_values[@]} * ${#discordant_z_values[@]}))
echo $total_iterations
current_iteration=0

# 循环运行LUMPY，并记录参数和VCF文件
for w in "${w_values[@]}"; do
    for msw in "${msw_values[@]}"; do
        for tt in "${tt_values[@]}"; do
            for back_distance in "${back_distance_values[@]}"; do
                for min_mapping_threshold in "${min_mapping_threshold_values[@]}"; do
                    for min_clip in "${min_clip_values[@]}"; do
                        for read_length in "${read_length_values[@]}"; do
                            for min_non_overlap in "${min_non_overlap_values[@]}"; do
                                for discordant_z in "${discordant_z_values[@]}"; do
                                            # 定义当前运行的结果文件名
                                    output_vcf="$output_dir/${sample_name}_variants_w${w}_msw${msw}_tt${tt}_bd${back_distance}_mmt${min_mapping_threshold}_mc${min_clip}_mean${mean}_stdev${stdev}_rl${read_length}_mno${min_non_overlap}_dz${discordant_z}.vcf"
                                            
                                            # 打印当前参数设置到日志文件
                                    echo "正在运行参数设置：w=$w, msw=$msw, tt=$tt, back_distance=$back_distance, min_mapping_threshold=$min_mapping_threshold, min_clip=$min_clip, mean=$mean, stdev=$stdev, read_length=$read_length, min_non_overlap=$min_non_overlap, discordant_z=$discordant_z" >> "$log_file"

                                            # 运行LUMPY并记录到VCF文件
									
									lumpy -w "$w" -msw "$msw" -tt "$tt" \
									-pe id:"$sample_name",read_group:readgroup1,bam_file:"$output_dir/$sample_name.discordants.bam",histo_file:"$output_dir/$sample_name.lib1.histo",mean:"$mean",stdev:"$stdev",read_length:"$read_length",min_non_overlap:"$min_non_overlap",discordant_z:"$discordant_z",back_distance:"$back_distance",weight:1,min_mapping_threshold:"$min_mapping_threshold" \
									-sr id:"$sample_name",bam_file:"$output_dir/$sample_name.splitters.bam",back_distance:"$back_distance",weight:1,min_mapping_threshold:"$min_mapping_threshold" > "$output_vcf" 2>&1

									evaf1.py "$output_vcf_path" "$output_vcf" output_f1_score.txt
									f1_score=$(sed -n '4p' output_f1_score.txt | tr -d '\n')
									precision=$(sed -n '5p' output_f1_score.txt | tr -d '\n')
									recall=$(sed -n '6p' output_f1_score.txt | tr -d '\n')

                                            # 将参数配置保存到CSV文件
									echo "$w,$msw,$tt,$back_distance,$min_mapping_threshold,$min_clip,$read_length,$min_non_overlap,$discordant_z,$mean,$stdev,$f1_score,$precision,$recall" >> "$csv_file"
									let current_iteration++
									percentage=$((100 * current_iteration / total_iterations))
									printf "\rProgress: [%-50s] %d%%" $(printf '%.0s#' $(seq 1 $((percentage/2)))) $percentage

                                done
                            done
                        done 
                    done
                done
            done
        done
    done
done
echo ""  # 这会在进度条后面打印一个新行
mv "${sample_name}_params.csv" /home/cloudam/simulat_2/csvpro/
echo "所有参数设置的运行已完成。" >> "$log_file"

import pandas as pd

# 替换以下路径为您的 CSV 文件的实际路径
file_path = r'C:\Users\1\Desktop\code_ppg\ICGC\1.csv'
# 使用 Pandas 读取 CSV 文件
df = pd.read_csv(file_path)

# 假设您的文件有两列，名为 'A' 和 'B'，下面的代码将找到这两列中相同的元素
common_elements_indices = df.index[df['icgc_donor_id'].isin(df['Sample'])].tolist()

print("第一列和第二列中相同元素对应的索引（编号）:")
print(common_elements_indices)
common_elements_indices = df.index[df['Sample'].isin(df['icgc_donor_id'])].tolist()

print("第一列和第二列中相同元素对应的索引（编号）:")
print(common_elements_indices)
df1 = pd.read_csv(r"C:\Users\1\Desktop\code_ppg\ICGC\2.csv")
new_df = df1.loc[common_elements_indices]

print("基于 common_elements_indices 的新 DataFrame:")
print(new_df)
# # 指定您想要保存新 CSV 文件的路径和文件名
# output_file_path = r'C:\Users\1\Desktop\code_ppg\ICGC\filtered_data.csv'

# # 将新的 DataFrame 保存为 CSV 文件
# new_df.to_csv(output_file_path, index=False)

# print(f"文件已成功保存到 {output_file_path}")

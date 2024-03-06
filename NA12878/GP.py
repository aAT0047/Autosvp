import numpy as np
import GPy
import pandas as pd
from scipy.optimize import minimize
from skopt import gp_minimize
from skopt.space import Real
from sklearn.model_selection import train_test_split
import random
from tqdm import tqdm
import concurrent.futures
import os
import re
from functools import partial
class BayesianOptimizer_f1:
    def __init__(self, data_path):
        self.data_path = data_path
        self.epsilon = 1e-5

    def custom_loss(self, y_true, y_pred, theta=0.05, A=1.5, B=2):
        diff = y_true - y_pred
        abs_diff = np.abs(diff)
        w = 1 + (A - 1) * (1 / (1 + np.exp(-10 * (abs_diff - theta))))
        lw = w * diff**2
        capped_abs_diff = np.clip(abs_diff, -1/B, 1/B)
        ls = np.exp(B * capped_abs_diff) * lw
        return np.mean(ls)

    def objective(self, params):
        self.model[:] = params
        y_pred, _ = self.model.predict(self.X_train)
        loss = self.custom_loss(self.y_train, y_pred)
        return loss

    def adjusted_bounds(self, col):
        col_min = self.X[col].min()
        col_max = self.X[col].max()
        if col_min == col_max:
            col_min -= self.epsilon
            col_max += self.epsilon
        return col_min, col_max

    def objective_fun(self, params):
        param_dict = {name: value for name, value in zip(self.X.columns, params)}
        param_df = pd.DataFrame([param_dict])
        y_pred, _ = self.model.predict(param_df.values)
        return -y_pred[0][0]

    def load_data(self):
        self.data = pd.read_csv(self.data_path)
        columns_to_drop = ["precision", "recall", "f1_score"]
        self.X = self.data.drop(columns=columns_to_drop)
        self.y = self.data['f1_score']
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X, self.y, test_size=0.2, random_state=42)
        self.X_train = self.X_train.values
        self.y_train = self.y_train.values.reshape(-1, 1)

    def train_model(self):
        input_dim = self.X_train.shape[1]
        kernel = GPy.kern.Matern52(input_dim=input_dim)+ GPy.kern.White(input_dim=input_dim, variance=1e-5)
        self.model = GPy.models.GPRegression(self.X_train, self.y_train, kernel)
        result = minimize(self.objective, self.model.param_array, method='L-BFGS-B')
        self.model[:] = result.x

    def optimize(self):
        space = [Real(*self.adjusted_bounds(col), name=col) for col in self.X.columns]
        self.result = gp_minimize(self.objective_fun, space, n_calls=100, random_state=0, verbose=True)
        sorted_indices = np.argsort(self.result.func_vals)
        # n_recommendations = random.randint(5, 10)
        n_recommendations =1
        # 创建一个 DataFrame 以收集所有的 recommended_params
        params_list = []

        for idx in sorted_indices[:n_recommendations]:
            recommended_params = {name: value for name, value in zip(self.X.columns, self.result.x_iters[idx])}
            print(f"推荐的参数组合 {idx+1}: {recommended_params}")
            params_list.append(recommended_params)

        df_opt_params = pd.DataFrame(params_list)


    # 保存 DataFrame 到 CSV
    # df_opt_params.to_csv(r"C:\Users\1\Desktop\som_ppg\csv\1_optparams.csv", index=False)

    def run(self):
        self.load_data()
        self.train_model()
        self.optimize()

    def save_optimized_parameters(self, save_path, sample_id):
        sorted_indices = np.argsort(self.result.func_vals)
        n_recommendations = random.randint(5, 10)
        n_recommendations =1
        # 创建一个 DataFrame 以收集所有的 recommended_params
        df_opt_params = pd.DataFrame()

        for idx in sorted_indices[:n_recommendations]:
            recommended_params = {name: value for name, value in zip(self.X.columns, self.result.x_iters[idx])}
            df_opt_params = df_opt_params.append(recommended_params, ignore_index=True)
        # 在DataFrame的最前面插入一列名为 'sample ID'，所有行都具有相同的 sample ID 值
        df_opt_params.insert(0, 'sample_name', sample_id)
        # 保存 DataFrame 到 CSV
        df_opt_params.to_csv(save_path, index=False)

class BayesianOptimizer_precision(BayesianOptimizer_f1):
    

    def load_data(self):
        self.data = pd.read_csv(self.data_path)
        columns_to_drop = ["precision", "recall","f1_score"]
        self.X = self.data.drop(columns=columns_to_drop)
        self.y = self.data['precision']
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X, self.y, test_size=0.2, random_state=42)
        self.X_train = self.X_train.values
        self.y_train = self.y_train.values.reshape(-1, 1)


class BayesianOptimizer_recall(BayesianOptimizer_f1):

    def load_data(self):
        self.data = pd.read_csv(self.data_path)
        columns_to_drop = ["precision", "recall","f1_score"]
        self.X = self.data.drop(columns=columns_to_drop)
        self.y = self.data['recall']
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X, self.y, test_size=0.2, random_state=42)
        self.X_train = self.X_train.values
        self.y_train = self.y_train.values.reshape(-1, 1)

def optimize_for_file(i, input_base_path, output_base_path, optimization_type='f1'):
    # 定义输入和输出文件的路径
    input_file_name = f"{i}_params.csv"
    input_path = os.path.join(input_base_path, input_file_name)
    output_file_name = f"{i}_optparams.csv"
    output_path = os.path.join(output_base_path, output_file_name)
    
    # 根据所选择的优化类型初始化不同的优化器
    if optimization_type == 'f1':
        optimizer = BayesianOptimizer_f1(data_path=input_path)
    elif optimization_type == 'precision':
        optimizer = BayesianOptimizer_precision(data_path=input_path)
    elif optimization_type == 'recall':
        optimizer = BayesianOptimizer_recall(data_path=input_path)
    else:
        raise ValueError(f"Unsupported optimization type: {optimization_type}")
    
    print(f"Optimizing for {input_file_name} using {optimization_type}...")
    
    try:
        optimizer.run()
        optimizer.save_optimized_parameters(output_path, i)
        print(f"Completed optimization for {input_file_name}!")
    except np.linalg.LinAlgError:
        print(f"Error occurred with {input_file_name}. Skipping this sample.")


# def optimize_for_filef1(i, input_base_path, output_base_path):
#     # 定义输入和输出文件的路径
#     input_file_name = f"{i}_params.csv"
#     input_path = os.path.join(input_base_path, input_file_name)
#     output_file_name = f"{i}_optparams.csv"
#     output_path = os.path.join(output_base_path, output_file_name)
#     print(f"Optimizing for {input_file_name}...")
#     optimizer = BayesianOptimizer_f1(data_path=input_path)
#     optimizer.run()
#     optimizer.save_optimized_parameters(output_path,i)
#     print(f"Completed optimization for {input_file_name}!")

# def optimize_for_fileprecision(i, input_base_path, output_base_path):
#     # 定义输入和输出文件的路径
#     input_file_name = f"{i}_params.csv"
#     input_path = os.path.join(input_base_path, input_file_name)
#     output_file_name = f"{i}_optparams.csv"
#     output_path = os.path.join(output_base_path, output_file_name)
#     print(f"Optimizing for {input_file_name}...")
#     optimizer = BayesianOptimizer_precision(data_path=input_path)
#     optimizer.run()
#     optimizer.save_optimized_parameters(output_path,i)
#     print(f"Completed optimization for {input_file_name}!")

# def optimize_for_filerecall(i, input_base_path, output_base_path):
#     # 定义输入和输出文件的路径
#     input_file_name = f"{i}_params.csv"
#     input_path = os.path.join(input_base_path, input_file_name)
#     output_file_name = f"{i}_optparams.csv"
#     output_path = os.path.join(output_base_path, output_file_name)
#     print(f"Optimizing for {input_file_name}...")
#     optimizer = BayesianOptimizer_recall(data_path=input_path)
#     optimizer.run()
#     optimizer.save_optimized_parameters(output_path,i)
#     print(f"Completed optimization for {input_file_name}!")

def main():
    # 获取文件夹中的所有文件名
    folder_path = '/home/cloudam/NA12878/csv'  # 请替换为您的文件夹路径
    file_names = os.listdir(folder_path)
    # 按照字母顺序排序
    sorted_file_names = sorted(file_names)
    # 提取文件名中的数字部分并保存到列表中
    numbers = [name.split("_")[0] for name in sorted_file_names]
    input_base_path = '/home/cloudam/NA12878/csv'
    output_base_path = '/home/cloudam/NA12878/csvtestdata'
    path = '/home/cloudam/NA12878/csvtestdata/'
    # 检查并创建输入目录
    if not os.path.exists(input_base_path):
        os.makedirs(input_base_path)
    # 检查并创建输出目录
    if not os.path.exists(output_base_path):
        os.makedirs(output_base_path)
    # 使用全部CPU核心进行并行运算。你可以通过修改max_workers参数来限制核心数量。
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # 使用 functools.partial 来预先设置 optimize_for_file 的 optimization_type 参数值
        func = partial(optimize_for_file, optimization_type='f1')
        
        # 使用tqdm提供一个进度条
        list(tqdm(executor.map(func, numbers, [input_base_path]*len(numbers), [output_base_path]*len(numbers)), total=len(numbers)))
        
    print("All optimizations completed!")
    # Generate the list of filenames with the full path
    files = [f"{output_base_path}/{i}_optparams.csv" for i in  numbers]

    # Check if the first file exists, and if so, read its content
    if os.path.exists(files[0]):
        merged_df = pd.read_csv(files[0])
    else:
        raise ValueError(f"File {files[0]} not found.")

    # From the second file onwards, only read the content and skip the header
    for file in tqdm(files[1:], desc="Merging files", unit="file"):
        if os.path.exists(file):
            df = pd.read_csv(file)
            merged_df = pd.concat([merged_df, df], ignore_index=True)
            # 删除原始的CSV文件
            os.remove(file)
    # Save the combined data to a new CSV file
    merged_df.to_csv(f"{path}/merged_dataf1.csv", index=False)

    print("All CSV files merged into merged_dataf1.csv")
        # 使用全部CPU核心进行并行运算。你可以通过修改max_workers参数来限制核心数量。
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # 使用 functools.partial 来预先设置 optimize_for_file 的 optimization_type 参数值为 'precision'
        func = partial(optimize_for_file, optimization_type='precision')
        
        # 使用tqdm提供一个进度条
        list(tqdm(executor.map(func, numbers, [input_base_path]*len(numbers), [output_base_path]*len(numbers)), total=len(numbers)))
        
    print("All optimizations completed!")


    # Generate the list of filenames with the full path
    files = [f"{output_base_path}/{i}_optparams.csv" for i in  numbers]

    # Check if the first file exists, and if so, read its content
    if os.path.exists(files[0]):
        merged_df = pd.read_csv(files[0])
    else:
        raise ValueError(f"File {files[0]} not found.")

    # From the second file onwards, only read the content and skip the header
    for file in tqdm(files[1:], desc="Merging files", unit="file"):
        if os.path.exists(file):
            df = pd.read_csv(file)
            merged_df = pd.concat([merged_df, df], ignore_index=True)
            # 删除原始的CSV文件
            os.remove(file)
    # Save the combined data to a new CSV file
    merged_df.to_csv(f"{path}/merged_datapre.csv", index=False)

    print("All CSV files merged into merged_datapre.csv")
    with concurrent.futures.ProcessPoolExecutor(max_workers=64) as executor:
        # 使用 functools.partial 来预先设置 optimize_for_file 的 optimization_type 参数值为 'recall'
        func = partial(optimize_for_file, optimization_type='recall')
        
        # 使用tqdm提供一个进度条
        list(tqdm(executor.map(func, numbers, [input_base_path]*len(numbers), [output_base_path]*len(numbers)), total=len(numbers)))
        
    print("All optimizations completed!")


    # Generate the list of filenames with the full path
    files = [f"{output_base_path}/{i}_optparams.csv" for i in  numbers]

    # Check if the first file exists, and if so, read its content
    if os.path.exists(files[0]):
        merged_df = pd.read_csv(files[0])
    else:
        raise ValueError(f"File {files[0]} not found.")

    # From the second file onwards, only read the content and skip the header
    for file in tqdm(files[1:], desc="Merging files", unit="file"):
        if os.path.exists(file):
            df = pd.read_csv(file)
            merged_df = pd.concat([merged_df, df], ignore_index=True)
            # 删除原始的CSV文件
            os.remove(file)
    # Save the combined data to a new CSV file
    merged_df.to_csv(f"{path}/merged_datarecall.csv", index=False)

    print("All CSV files merged into merged_datarec.csv")

if __name__ == "__main__":
    main()
import random
import numpy as np
import re
def find_longest_consecutive_sequences(original_list):
    nums = [x for x in original_list if isinstance(x, int)]
    all_longest_sequences = []
    current_sequence = []

    for i in range(len(nums)):
        # 如果是列表的第一个数字或者是前一个数字的连续值
        if not current_sequence or nums[i] == current_sequence[-1] + 1:
            current_sequence.append(nums[i])
        else:
            # 更新最长的连续子序列列表
            if not all_longest_sequences or len(current_sequence) > len(all_longest_sequences[0]):
                all_longest_sequences = [current_sequence]
            elif len(current_sequence) == len(all_longest_sequences[0]):
                all_longest_sequences.append(current_sequence)

            # 重置当前序列
            current_sequence = [nums[i]]

    # 检查最后一个序列
    if not all_longest_sequences or len(current_sequence) > len(all_longest_sequences[0]):
        all_longest_sequences = [current_sequence]
    elif len(current_sequence) == len(all_longest_sequences[0]):
        all_longest_sequences.append(current_sequence)

    return all_longest_sequences


# 提取子列表
def extract_sublists_by_indices(matrix, index_list):
    selected_sublists = [sublist for sublist in matrix if sublist[-1] in index_list]
    return selected_sublists

# 提取MAP边标签
def get_last_numbers(matrix):
    last_nums = [edge[-1] for edge in matrix]
    return last_nums

# 生成主链kemers n=k-1
def sliding_window(string, n):
    result = []
    for i in range(len(string) - n + 1):
        substring = string[i:i + n].upper()
        result.append(substring)
    return result

# 定位reads在染色体的位置
def find_substring_indices(main_string, data):
    first_elements = [item[0][0] for item in data]
    end_elements =  [item[-1] for item in data][0]
    sub_string = ''.join(first_elements)
    start_index = main_string.find(sub_string) + 1
    end_index = start_index + len(sub_string) - 1
    chr_start_index = start_index
    return chr_start_index,start_index-1,end_index 

def find_indel_indices(main_string, data):
    
    sub_string = data[0][1::]
    start_index = main_string.find(sub_string[0]) + 1
    end_index = start_index + len(sub_string) - 1
    # chr_start_index = end_elements-start_index+1
    chr_start_index = start_index
    return chr_start_index,start_index-1,end_index 


def get_map_edge_labels(edge_l, nodes):
    edge_labels = []
    for i in range(len(nodes) - 1):
        src_node = nodes[i]
        dst_node = nodes[i + 1]
        
        for edge_info in edge_l:
            current_edge_key, label = edge_info
            current_src_node, current_dst_node, current_key = current_edge_key
            
            if current_src_node == src_node and current_dst_node == dst_node:
                edge_labels.extend([[src_node, dst_node, current_key, label]])
                
    return edge_labels

# 统计list中 invert 和dele 个数

def count_elements(data):
    # 创建一个空字典来存储统计结果
    invert_count = {}
    delete_count = {}

    # 遍历列表元素
    i = 0
    while i < len(data):
        element = data[i]
        if isinstance(element, str):
            if element.startswith('invert_'):
                # 提取 'invert_' 后面的数字
                invert_number = int(element.split('_')[1])

                # 检查该数字是否存在于列表中
                if invert_number in data:
                    # 如果存在，统计数量
                    if invert_number not in invert_count:
                        invert_count[invert_number] = 2
                    else:
                        invert_count[invert_number] += 1
                # 跳过处理过的元素
                i += 1
            elif element.startswith('delete_'):
                # 提取 'delete_' 后面的两个数字
                delete_numbers = element.split('_')[1:]
                delete_numbers = [int(num) for num in delete_numbers]

                # 检查两个数字中的任何一个是否存在于列表中
                if any(num in data for num in delete_numbers):
                    # 统计以 'delete_2_7' 开头的元素数量
                    delete_prefix = f'delete_{delete_numbers[0]}_{delete_numbers[1]}'
                    delete_count[delete_prefix] = delete_count.get(delete_prefix, 0) + 1
               
                
                # 跳过处理过的元素
                i += 1
            else:
                i += 1
        else:
            i += 1
    for key in delete_count:
        delete_count[key] += 1
    combined_dict = {**invert_count, **delete_count}

    # 获取合并后字典的最大值
    if combined_dict:
        max_value = max(combined_dict.values())

        # 找到最大值对应的键
        max_keys = [key for key, value in combined_dict.items() if value == max_value][0]
        
        if isinstance( max_keys, int) is False:
            key = [key.split('_')[-1] for key in max_keys][-1]
            max_keys = key
    
    
        return max_keys, max_value

    else:
        # 处理字典为空的情况
        return 0, 0
def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-1):
    # 创建一个矩阵用于存储得分，以及一个矩阵用于存储回溯指针
    rows, cols = len(seq1) + 1, len(seq2) + 1
    scores = np.zeros((rows, cols))
    pointers = np.zeros((rows, cols), dtype=int)

    # 填充矩阵
    for i in range(1, rows):
        for j in range(1, cols):
            diag = scores[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            up = scores[i - 1][j] + gap
            left = scores[i][j - 1] + gap
            scores[i][j] = max(diag, up, left, 0)
            if scores[i][j] == diag:
                pointers[i][j] = 1  # 对角线
            elif scores[i][j] == up:
                pointers[i][j] = 2  # 上方
            elif scores[i][j] == left:
                pointers[i][j] = 3  # 左侧

    # 找到最高分及其位置
    max_score = 0
    max_pos = (0, 0)
    for i in range(1, rows):
        for j in range(1, cols):
            if scores[i][j] > max_score:
                max_score = scores[i][j]
                max_pos = (i, j)

    # 回溯以找到比对
    i, j = max_pos
    alignment_seq1 = ''
    alignment_seq2 = ''
    while i > 0 and j > 0 and scores[i][j] > 0:
        if pointers[i][j] == 1:
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            i -= 1
            j -= 1
        elif pointers[i][j] == 2:
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = '-' + alignment_seq2
            i -= 1
        elif pointers[i][j] == 3:
            alignment_seq1 = '-' + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            j -= 1
        else:
            break

    return alignment_seq1, alignment_seq2, max_score

def generate_unique_random_numbers(count, range_start, range_end):
    numbers = list(range(range_start, range_end + 1))
    random.shuffle(numbers)
    return numbers[:count]

# 生成两组不相同的 100 个随机数
random_numbers1 = generate_unique_random_numbers(10000, 0, 20000)
random_numbers2 = generate_unique_random_numbers(20000, 0, 20000)
# import random

def generate_random_sequence(length):
    bases = ['A', 'T', 'G', 'C']
    sequence = [random.choice(bases) for _ in range(length)]
    return ''.join(sequence)

# 生成 100 个不同长度的随机序列
random_sequences = [generate_random_sequence(random.randint(1, 50)) for _ in range(1000)]



def generate_random_sequence(length):
    bases = ['A', 'T', 'G', 'C']
    sequence = ''.join(random.choice(bases) for _ in range(length))
    return sequence


# 组成reads子图
import networkx as nx

def create_labeled_multidigraph(target_sequence, k, edge_index, edge_labels):
    # 创建一个 MultiDiGraph
    read = nx.MultiDiGraph()

    # 将节点添加到图中
    read_main_nodes = sliding_window(target_sequence, k - 1)
    for node in read_main_nodes:
        read.add_node(node)

    # 将边添加到图中
    for i in range(len(read_main_nodes) - 1):
        read.add_edge(read_main_nodes[i], read_main_nodes[i + 1])

    # 创建一个字典来存储边标签
    read_edge_labels = {}
    
    for i in edge_index:
        # 给定值
        target_value = i
        # 遍历字典，查找包含给定值的键
        matching_keys = [(u, v, key) for (u, v, key), value in edge_labels.items() if value == target_value]
        if matching_keys:
            # 获取第一个匹配键的源节点、目标节点和键（边标签）
            u, v, key = matching_keys[0]
            # 将边标签添加到字典中
            read_edge_labels[(u, v, key)] = i

    # 使用 set_edge_attributes 函数将边标签添加到图中
    nx.set_edge_attributes(read, read_edge_labels, "label")

    return read,read_edge_labels,read.edges,read_main_nodes 


# # # 寻找最长连续子序列的前后延长序列
# def find_min_key_for_value(dictionary, target_value):
#     def closest_value_to_target(lst, target):
#         numeric_values = [x for x in lst if isinstance(x, (int, float))]
#         if not numeric_values:
#             return float('inf')
#         return min(numeric_values, key=lambda x: abs(x - target))
    
#     closest_key = None
#     closest_distance = float('inf')

#     for k, v in dictionary.items():
#         if isinstance(v, list):
#             closest = closest_value_to_target(v, target_value)
#             distance = abs(closest - target_value)
#             if distance < closest_distance:
#                 closest_distance = distance
#                 closest_key = k

#     return closest_key

# def keys_info_containing_smaller_value(d, value, max_k):
#     """Return all keys in the dictionary (with key < max_k) that contain values smaller than the provided value."""
#     # return [k for k, v in d.items() if any(i < value for i in v if isinstance(i, int)) and k < max_k]
#     return [k for k, v in d.items() if value is not None and any(i <= value for i in v if isinstance(i, int)) and (max_k is not None and k <= max_k)]

# def get_sorted_infer_result(array_dict, initial_value):
#     array_dict_copy = array_dict.copy()
#     current_value = initial_value
#     last_key = find_min_key_for_value(array_dict_copy , current_value)

#     result = {}

#     while True:
#         keys = keys_info_containing_smaller_value(array_dict_copy, current_value, last_key)
    
#         # For each key, find the largest value that's less than current_value
#         # potential_pairs = [(k, max([v for v in array_dict_copy[k] if isinstance(v, int) and v < current_value])) for k in keys]
#         potential_pairs = [(k, max([v for v in array_dict_copy[k] if v is not None and isinstance(v, int) and v <= current_value], default=None))\
#             for k in keys if any(v is not None and isinstance(v, int) and v <= current_value for v in array_dict_copy[k])]
#         # Sort by key in descending order and then by value in descending order
#         potential_pairs = sorted(potential_pairs, key=lambda x: (-x[0], -x[1]))

#         if potential_pairs:
#             chosen_key, chosen_value = potential_pairs[0]
            
#             # If we encounter a problematic value, we remove the key and continue
#             if not result:
#                 pass
#             else:
#                 if 'invert_' + str(result[min(result.keys())][-1]-1) in array_dict.get(chosen_key, []):
#                     del array_dict_copy[chosen_key]
#                     continue
                
#                 search_pattern = 'delete_' + str(result[min(result.keys())][-1])
#                 contains_pattern = any(re.search(search_pattern, value) for value in array_dict.get(chosen_key, []) if isinstance(value, str))
#                 if contains_pattern:
#                     del array_dict_copy[chosen_key]
#                     continue
            
#             result[chosen_key] = [chosen_value]
#             del array_dict_copy[chosen_key]
#             last_key = chosen_key
#             current_value = chosen_value
#         else:
#             # If there are no suitable keys, we stop
#             break

#     return dict(sorted(result.items()))

# def keys_after_containing_greater_value(d, value, min_k):
#     """Return all keys in the dictionary (with key > min_k) that contain values greater than the provided value."""
#     return [k for k, v in d.items() if value is not None and any(i >= value for i in v if isinstance(i, int)) and (min_k is not None and k >= min_k)]



# def get_sorted_after_result(array_dict, initial_value):
#     array_dict_copy = array_dict.copy()
#     current_value = initial_value
#     last_key = find_min_key_for_value(array_dict_copy, current_value)

#     result = {}

#     while True:
#         keys = keys_after_containing_greater_value(array_dict_copy, current_value, last_key)
#         potential_pairs = [(k, min([v for v in array_dict_copy[k] if v is not None and isinstance(v, int) and v >= current_value], default=None))\
#             for k in keys if any(v is not None and isinstance(v, int) and v >= current_value for v in array_dict_copy[k])]
#         # potential_pairs = [(k, min([v for v in array_dict_copy[k] if isinstance(v, int) and v > current_value])) for k in keys]
#         potential_pairs = sorted(potential_pairs, key=lambda x: (x[0], x[1]))  # sort by key first, then by value

#         if potential_pairs:
#             chosen_key, chosen_value = potential_pairs[0]
#             if not result:
#                 pass
#             else:
#                 if 'invert_' + str(result[max(result.keys())][-1]) in array_dict.get(chosen_key, []):
#                     del array_dict_copy[chosen_key]
#                     continue
                    
#                 search_pattern = 'delete_' + str(result[max(result.keys())][-1])
#                 contains_pattern = any(re.search(search_pattern, value) for value in array_dict.get(chosen_key, []) if isinstance(value, str))
#                 if contains_pattern:
#                     del array_dict_copy[chosen_key]
#                     continue
#             result[chosen_key] = [chosen_value]
#             del array_dict_copy[chosen_key]
#             last_key = chosen_key
#             current_value = chosen_value
#         else:
#             # If there are no suitable keys, we stop
#             break

#     return dict(sorted(result.items()))

# # 生成最长匹配列表字典
# def filter_edges(read_edges, edge_labels):
#     result = {index + 1: sub_tuple for index, sub_tuple in enumerate(read_edges)}
#     # print(result)
#     read_dict = {}
#     for key, value in result.items():
#         # read_dict[key] = [edge_labels[k] for k in edge_labels if value == k[:2]]
#         read_dict[key] = []
#         for k in edge_labels:
#             if value == k[:2]:
#                 read_dict[key].append(edge_labels[k])
#     filtered_dict = {}

#     return read_dict

# def extract_values_in_range(input_dict, start_range, end_range):
#     # 创建一个新的字典来存储值在指定范围内的项
#     result_dict = {}

#     # 遍历原始字典的键值对
#     for key, value_list in input_dict.items():
#         # 检查值是否在指定的范围内
#         filtered_values = [v for v in value_list if start_range <= v <= end_range]

#         # 如果有符合条件的值，将其存储在新的字典中
#         if filtered_values:
#             result_dict[key] = filtered_values

#     return result_dict
# # 生成连续的连接边
# def create_string_pairs_with_repeat(strings):
#     # 使用循环生成字符串对，每个字符串都在两个子列表中出现
#     string_pairs = []

#     for i in range(1, len(strings)):
#         pair1 = (strings[i-1], strings[i])
#         string_pairs.append(pair1)

#     return string_pairs

# def find_matching_subdicts(input_dict, values_to_find):
#     matching_subdicts = {}
#     for key, value in input_dict.items():
#         if value in values_to_find:
#             matching_subdicts[key] = value
#     return matching_subdicts

# def find_middle_reads(my_dict,first_index,last_index):
#     keys_with_value_1 = [key for key, value in my_dict.items() if value == [first_index]]
    
#     for start_key in keys_with_value_1:
#         current_key = start_key
#         for i in range(first_index+1, last_index+1):  # 从2开始，因为1已经是起点了
#             current_key += 1
#             if my_dict.get(current_key) != [i]:
#                 break
#             # return {}    
#         else:
#             # 如果for循环没有提前中断，说明找到了连续的keys
#             keys_range = list(range(start_key, current_key + 1))
#             return {key: my_dict[key] for key in keys_range}
    




# # 获得最终read
# def process_read_data(read_main_nodes, edge_labels, edge_index):
#     result_read = {}
#     read_edges = create_string_pairs_with_repeat(read_main_nodes)
#     pre_midd_r = filter_edges(read_edges,find_matching_subdicts(edge_labels,edge_index))
    
#     max_read_dict = filter_edges(read_edges, edge_labels)

#     middle_read = find_middle_reads(pre_midd_r, edge_index[0], edge_index[-1])
#     first_key, first_value = next(iter(middle_read.items()))
#     last_key, last_value = max(middle_read.items())
#     af = {key: value for key, value in max_read_dict.items() if  key > last_key }
#     after_read = get_sorted_after_result(af, last_value[0]+1)
#     inf = {key: value for key, value in max_read_dict.items() if  key <first_key}
#     infer_read = get_sorted_infer_result(inf \
#         , first_value[0]-1)
#     result_read.update(infer_read)
#     result_read.update(middle_read)
#     result_read.update(after_read)
#     output_dict = {}
#     for key, value in result_read.items():
#         if value not in output_dict.values():
#             output_dict[key] = value

#     return output_dict,read_edges


import numpy as np
import re
import multiprocessing

def find_min_key_for_value(dictionary, target_value):
    numeric_values = [v for v in dictionary.values() if isinstance(v, (int, float))]
    
    if not numeric_values:
        return None

    closest_key = min(dictionary, key=lambda k: abs(dictionary[k] - target_value))
    return closest_key

def keys_info_containing_smaller_value(d, value, max_k):
    return [k for k, v in d.items() if v is not None and any(i <= value for i in v if isinstance(i, int)) and (max_k is None or k <= max_k)]

def get_sorted_infer_result(array_dict, initial_value):
    array_dict_copy = array_dict.copy()
    current_value = initial_value
    result = {}

    while True:
        keys = keys_info_containing_smaller_value(array_dict_copy, current_value, None)
        
        if not keys:
            break

        potential_pairs = [(k, max(v, default=None)) for k, v in array_dict_copy.items() if k in keys]
        potential_pairs = sorted(potential_pairs, key=lambda x: (-x[0], -x[1]))

        if potential_pairs:
            chosen_key, chosen_value = potential_pairs[0]

            if not result:
                pass
            else:
                if 'invert_' + str(result[min(result.keys())][-1]-1) in array_dict.get(chosen_key, []):
                    del array_dict_copy[chosen_key]
                    continue
                
                search_pattern = 'delete_' + str(result[min(result.keys())][-1])
                contains_pattern = any(re.search(search_pattern, value) for value in array_dict.get(chosen_key, []) if isinstance(value, str))
                if contains_pattern:
                    del array_dict_copy[chosen_key]
                    continue
            
            result[chosen_key] = [chosen_value]
            del array_dict_copy[chosen_key]
            current_value = chosen_value
        else:
            break

    return dict(sorted(result.items()))

def keys_after_containing_greater_value(d, value, min_k):
    return [k for k, v in d.items() if v is not None and any(i >= value for i in v if isinstance(i, int)) and (min_k is None or k >= min_k)]

def get_sorted_after_result(array_dict, initial_value):
    array_dict_copy = array_dict.copy()
    current_value = initial_value
    result = {}

    while True:
        keys = keys_after_containing_greater_value(array_dict_copy, current_value, None)
        potential_pairs = [(k, min(v, default=None)) for k, v in array_dict_copy.items() if k in keys]
        potential_pairs = sorted(potential_pairs, key=lambda x: (x[0], x[1]))

        if potential_pairs:
            chosen_key, chosen_value = potential_pairs[0]

            if not result:
                pass
            else:
                if 'invert_' + str(result[max(result.keys())][-1]) in array_dict.get(chosen_key, []):
                    del array_dict_copy[chosen_key]
                    continue
                    
                search_pattern = 'delete_' + str(result[max(result.keys())][-1])
                contains_pattern = any(re.search(search_pattern, value) for value in array_dict.get(chosen_key, []) if isinstance(value, str))
                if contains_pattern:
                    del array_dict_copy[chosen_key]
                    continue

            result[chosen_key] = [chosen_value]
            del array_dict_copy[chosen_key]
            current_value = chosen_value
        else:
            break

    return dict(sorted(result.items()))

def filter_edges(read_edges, edge_labels):
    result_dict = {}
    for index, sub_tuple in enumerate(read_edges, start=1):
        result_dict[index] = [edge_labels[k] for k in edge_labels if sub_tuple == k[:2]]
    return result_dict

def extract_values_in_range(input_dict, start_range, end_range):
    result_dict = {}
    for key, value_list in input_dict.items():
        filtered_values = [v for v in value_list if start_range <= v <= end_range]
        if filtered_values:
            result_dict[key] = filtered_values
    return result_dict

def create_string_pairs_with_repeat(strings):
    string_pairs = [(strings[i-1], strings[i]) for i in range(1, len(strings))]
    return string_pairs

def find_matching_subdicts(input_dict, values_to_find):
    matching_subdicts = {key: value for key, value in input_dict.items() if value in values_to_find}
    return matching_subdicts

def find_middle_reads(my_dict, first_index, last_index):
    keys_with_value_1 = [key for key, value in my_dict.items() if value == [first_index]]

    for start_key in keys_with_value_1:
        current_key = start_key
        for i in range(first_index+1, last_index+1):
            current_key += 1
            if my_dict.get(current_key) != [i]:
                break
        else:
            keys_range = list(range(start_key, current_key + 1))
            return {key: my_dict[key] for key in keys_range}

def process_read_data(read_main_nodes, edge_labels, edge_index):
    result_read = {}
    read_edges = create_string_pairs_with_repeat(read_main_nodes)
    pre_midd_r = filter_edges(read_edges, find_matching_subdicts(edge_labels, edge_index))
    
    max_read_dict = filter_edges(read_edges, edge_labels)

    middle_read = find_middle_reads(pre_midd_r, edge_index[0], edge_index[-1])
    first_key, first_value = next(iter(middle_read.items()))
    last_key, last_value = max(middle_read.items())
    af = {key: value for key, value in max_read_dict.items() if key > last_key}
    after_read = get_sorted_after_result(af, last_value[0] + 1)
    inf = {key: value for key, value in max_read_dict.items() if key < first_key}
    infer_read = get_sorted_infer_result(inf, first_value[0] - 1)
    result_read.update(infer_read)
    result_read.update(middle_read)
    result_read.update(after_read)
    output_dict = {}
    for key, value in result_read.items():
        if value not in output_dict.values():
            output_dict[key] = value

    return output_dict, read_edges




# 读取fq文件
def read_reference_genome(reference_genome_file):
    genome_sequence = ""

    with open(reference_genome_file, 'r') as file:
        # 跳过第一行描述信息
        file.readline()
        
        for line in file:
            # 移除行末尾的换行符
            line = line.strip()
            genome_sequence += line

    return genome_sequence
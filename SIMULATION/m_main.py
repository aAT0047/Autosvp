import hashlib
import networkx as nx
from usingtools import *
from PAN_GENE_GRAPH import DeBruijnGraph
from NWAligner import *
from sv import *
import time

def process_sequence_matching(target_sequence, sequence, lennode):
    chr_number = None  # 初始化 chr_number
    chr_start_index = None  # 初始化 chr_start_index
    edge_index = None  # 初始化 edge_index
    matching_paths=[]
    for chromosome_name, chromosome in chromosomes.items():
        start_time = time.time()
        matching_path = chromosome.find_sequence_location(target_sequence, lennode=k-1)[0]
            # 记录结束时间
        end_time = time.time()
        
        # 计算运行时间
        elapsed_time = end_time - start_time
        # print(f" 所用时间：{elapsed_time} 秒")
        matching_paths.append(matching_path )
        max_length = max(len(sublist) for sublist in matching_paths)
        indices_of_longest = [i for i, sublist in enumerate(matching_paths) if len(sublist) == max_length]
        # matching_paths = [ matching_paths[i] for i in indices_of_longest]
        chr_index = []
    for i in indices_of_longest:
        chr_key=[]
        edge_labels = chromosomes[list(chromosomes.keys())[i]].get_edge_labels()
        if len(matching_paths[i]) > (len(target_sequence)+ k-1) * 0.006:
            edge_label_list = get_last_numbers(get_map_edge_labels(edge_labels.items(), matching_paths[i]))
            max_keys, max_value = count_elements(edge_label_list)
            edge_index = find_longest_consecutive_sequences(edge_label_list)
            
            if max_value/2 < len(edge_index[0]):
                edgepath = extract_sublists_by_indices(get_map_edge_labels(edge_labels.items(), matching_paths[i]), edge_index[0])
                chr_start_index_main, sta_main, end_main = find_substring_indices(target_sequence, edgepath)
                # chr_number = i+1
                aligned_seq1, aligned_seq2, score_main = smith_waterman(
                    sequence[i][0][ chr_start_index_main-1: chr_start_index_main-1+sta_main]+\
                        sequence[i][0][chr_start_index_main-1-end_main:len(target_sequence)],\
                          target_sequence[0:sta_main]+target_sequence[end_main::])
                chr_key.insert(0,i+1)
                chr_key.insert(1,chr_start_index_main)
                chr_key.insert(2,score_main) 
            elif max_value/2 > len(edge_index[0]):
                edgepath_indel = extract_sublists_by_indices(get_map_edge_labels(edge_labels.items(), matching_paths[i]), [max_keys])
                chr_start_index_indel, sta_indel, end_indel = find_indel_indices(target_sequence, edgepath_indel)
                if sta_indel < 0.5 * len(target_sequence):
                    _, _, score_indel = smith_waterman(
                                                    sequence[i][0][chr_start_index_indel - 1:chr_start_index_indel - 1 + sta_indel]\
                                                        ,target_sequence[0:sta_indel])
                else:
                    _, _, score_indel = smith_waterman(
                                                    sequence[i][0][chr_start_index_indel - 1 + end_indel:len(target_sequence)]\
                                                        ,target_sequence[end_indel::])
                chr_key.insert(0,i+1)
                chr_key.insert(1,chr_start_index_indel)
                chr_key.insert(2,score_indel)

            else:
                edgepath_main = extract_sublists_by_indices(get_map_edge_labels(edge_labels.items(), matching_paths[i]), edge_index[0])
                chr_start_index_main, sta_main, end_main = find_substring_indices(target_sequence, edgepath)
                edgepath_indel = extract_sublists_by_indices(get_map_edge_labels(edge_labels.items(), matching_paths[i]), [max_keys])
                chr_start_index_indel, sta_indel, end_indel = find_indel_indices(target_sequence, edgepath)

                aligned_seq1, aligned_seq2, score_main = smith_waterman(
                    sequence[i][0][chr_start_index_main - 1: chr_start_index_main - 1 + sta_main] \
                        + sequence[i][0][chr_start_index_main - 1 - end_main:len(target_sequence)],\
                          target_sequence[0:sta_main]+target_sequence[end_main::])
                
                if sta_indel < 0.5 * len(target_sequence):
                    _, _, score_indel = smith_waterman(
                                                    sequence[i][0][chr_start_index_indel - 1:chr_start_index_indel - 1 + sta_indel]\
                                                        ,target_sequence[0:sta_indel])
                else:
                    _, _, score_indel = smith_waterman(
                                                    sequence[i][0][chr_start_index_indel - 1 + end_indel:len(target_sequence)]\
                                                        ,target_sequence[end_indel::])
                
                if score_indel > score_main:
                    chr_start_index = chr_start_index_indel
                    key = chr_start_index
                    chr_key.insert(0,i+1)
                    chr_key.insert(1,key)
                    chr_key.insert(2,score_indel)
                else:
                    chr_start_index = chr_start_index_main
                    key = chr_start_index
                    chr_key.insert(0,i+1)
                    chr_key.insert(1,key)
                    chr_key.insert(2,score_main)
            chr_index.append( chr_key)
            csi =  max(chr_index, key=lambda x: x[-1])   
            chr_number = csi[0]
            chr_start_index = csi[1]
        
    
            # edge_labels = chromosomes[list(chromosomes.keys())[chr_number-1]].get_edge_labels()
            edge_label_list = get_last_numbers(get_map_edge_labels(edge_labels.items(), matching_paths[chr_number-1]))
            edge_index = find_longest_consecutive_sequences(edge_label_list)[0]

        else:
                None   
    return chr_number, chr_start_index ,edge_index ,matching_paths
import time
import numpy as np

# def process_sequence_matching(target_sequence, sequence, lennode, k, chromosomes):
#     matching_paths = []
    
#     # 计算 target_sequence 的长度，以减少重复计算
#     target_len = len(target_sequence)
    
#     # 提前获取 chromosomes 的键和边缘标签
#     chromosome_keys = list(chromosomes.keys())
#     edge_labels = chromosomes[chromosome_keys[0]].get_edge_labels()
    
#     # 计算 edge_labels 中的值，以减少重复计算
#     edge_label_list = get_last_numbers(get_map_edge_labels(edge_labels.items(), sequence))
    
#     for i, chromosome_name in enumerate(chromosome_keys):
#         start_time = time.time()
        
#         matching_path = chromosomes[chromosome_name].find_sequence_location(target_sequence, lennode=k-1)[0]
#         matching_paths.append(matching_path)
        
#         # 记录结束时间
#         end_time = time.time()
        
#         # 计算运行时间
#         elapsed_time = end_time - start_time
#         print(f"所用时间：{elapsed_time}秒")
    
#     # 找到最长匹配路径的索引
#     max_length = max(len(sublist) for sublist in matching_paths)
#     indices_of_longest = [i for i, sublist in enumerate(matching_paths) if len(sublist) == max_length]
    
#     chr_index = []
#     for i in indices_of_longest:
#         chr_key = []
#         edge_label_list_i = get_last_numbers(get_map_edge_labels(edge_labels.items(), matching_paths[i]))
#         max_keys, max_value = count_elements(edge_label_list_i)
#         edge_index = find_longest_consecutive_sequences(edge_label_list_i)
        
#         if max_value / 2 < len(edge_index[0]):
#             edgepath = extract_sublists_by_indices(get_map_edge_labels(edge_labels.items(), matching_paths[i]), edge_index[0])
#             chr_start_index_main, sta_main, end_main = find_substring_indices(target_sequence, edgepath)
#             aligned_seq1, aligned_seq2, score_main = smith_waterman(
#                 sequence[i][0][chr_start_index_main-1:chr_start_index_main-1+sta_main] +
#                 sequence[i][0][chr_start_index_main-1-end_main:chr_start_index_main-1+target_len],
#                 target_sequence[0:sta_main] + target_sequence[end_main:])
            
#             chr_key.extend([i+1, chr_start_index_main, score_main])
#         elif max_value / 2 > len(edge_index[0]):
#             edgepath_indel = extract_sublists_by_indices(get_map_edge_labels(edge_labels.items(), matching_paths[i]), [max_keys])
#             chr_start_index_indel, sta_indel, end_indel = find_indel_indices(target_sequence, edgepath_indel)
            
#             if sta_indel < 0.5 * target_len:
#                 _, _, score_indel = smith_waterman(
#                     sequence[i][0][chr_start_index_indel-1:chr_start_index_indel-1+sta_indel],
#                     target_sequence[0:sta_indel])
#             else:
#                 _, _, score_indel = smith_waterman(
#                     sequence[i][0][chr_start_index_indel-1+end_indel:chr_start_index_indel-1+target_len],
#                     target_sequence[end_indel:])
            
#             chr_key.extend([i+1, chr_start_index_indel, score_indel])
#         else:
#             edgepath_main = extract_sublists_by_indices(get_map_edge_labels(edge_labels.items(), matching_paths[i]), edge_index[0])
#             chr_start_index_main, sta_main, end_main = find_substring_indices(target_sequence, edgepath_main)
#             edgepath_indel = extract_sublists_by_indices(get_map_edge_labels(edge_labels.items(), matching_paths[i]), [max_keys])
#             chr_start_index_indel, sta_indel, end_indel = find_indel_indices(target_sequence, edgepath_indel)
            
#             aligned_seq1, aligned_seq2, score_main = smith_waterman(
#                 sequence[i][0][chr_start_index_main-1:chr_start_index_main-1+sta_main] +
#                 sequence[i][0][chr_start_index_main-1-end_main:chr_start_index_main-1+target_len],
#                 target_sequence[0:sta_main] + target_sequence[end_main:])
            
#             if sta_indel < 0.5 * target_len:
#                 _, _, score_indel = smith_waterman(
#                     sequence[i][0][chr_start_index_indel-1:chr_start_index_indel-1+sta_indel],
#                     target_sequence[0:sta_indel])
#             else:
#                 _, _, score_indel = smith_waterman(
#                     sequence[i][0][chr_start_index_indel-1+end_indel:chr_start_index_indel-1+target_len],
#                     target_sequence[end_indel:])
            
#             if score_indel > score_main:
#                 chr_start_index = chr_start_index_indel
#                 chr_key.extend([i+1, chr_start_index, score_indel])
#             else:
#                 chr_start_index = chr_start_index_main
#                 chr_key.extend([i+1, chr_start_index, score_main])
        
#         chr_index.append(chr_key)
    
#     # 从 chr_index 中找到具有最高分数的染色体
#     chr_index = np.array(chr_index)
#     max_score_index = np.argmax(chr_index[:, 2].astype(float))
    
#     chr_number = chr_index[max_score_index, 0]
#     chr_start_index = int(chr_index[max_score_index, 1])
    
#     return chr_number, chr_start_index, edge_index, matching_paths[chr_number - 1]

# # 假设你已经有了所需的辅助函数和数据结构
# # 请注意，还需要传递参数 k 和 chromosomes 到函数中


def process_translocation(read_edges,result_read,target_sequence, sequence, k, chromosomes, chr_number):
    length  = len(read_edges) - len(result_read)
    start = list(result_read.values())[-1]
    tra_dict = {}
    
    sec_chr_number, _, sec_edge_index, sec_matching_paths = process_sequence_matching(target_sequence[start[0]+1::], sequence, lennode=k-1)
    sec_edge_labels = chromosomes[list(chromosomes.keys())[sec_chr_number-1]].get_edge_labels()
    
    sec_read, _, _, sec_read_main_nodes = create_labeled_multidigraph(target_sequence[start[0]+1::], k, sec_edge_index, sec_edge_labels)
    sec_result_read, sec_read_edges = process_read_data(sec_read_main_nodes, sec_edge_labels, sec_edge_index)
    end = next(iter(sec_result_read.values()))
    
    tra_dict['chr'+str(chr_number)] = (start[0]+1, str(']'+str(sec_chr_number)+':'+str(end[0])+']'))
    tra_dict['chr'+str(sec_chr_number)] = (end[0], str('['+str(chr_number)+':'+str(start[0]+1)+'['))

    return tra_dict

# def find_inverted_sequences(k, readdictt, previous_dict, next_dict, chromosomes, chr_number):
    """
    Given some parameters, this function attempts to identify and return inverted sequences.
    
    Args:
    - k (int): Presumably the k-mer length.
    - readdictt (dict): Dictionary related to sequence reads.
    - previous_dict (dict): Dictionary representing previous sequences or nodes.
    - next_dict (dict): Dictionary representing next sequences or nodes.
    - chromosomes (dict): Dictionary containing chromosome data.
    - chr_number (int): The chromosome number to consider.
    
    Returns:
    - dict: Dictionary containing the inverted sequences found.
    """
    
    undefined_seq = undefined_dicts(chr_number,k, readdictt, previous_dict, next_dict)
    invert_dict = {}
    
    items_previous = list(previous_dict.items())
    items_next = list(next_dict.items())
    
    for i in range(len(items_previous)):
        key_previous, value_previous = items_previous[i]
        key_next, value_next = items_next[i]

        reversed_sequence = undefined_seq[i][::-1]
        mat = chromosomes[list(chromosomes.keys())[chr_number-1]].find_sequence_location(reversed_sequence, lennode=k-1)[0]
        
        if len(mat) >= (len(reversed_sequence)-k+1)*0.8:
            result = 'invert'+str(value_previous[0]+k)
            count = value_next[0]-1-(value_previous[0]+k)+1
            invert_dict[result] = (chr_number,value_previous[0]+k, value_next[0]-1,count)
    
    return invert_dict

def process_fasta_file(fasta_filename):
    global  sequence,  k ,chromosomes
    fragments = []
    result_dict = {}
    
    for record in SeqIO.parse(fasta_filename, "fasta"):
        fragments.append(record.seq)
    
    for i, fragment in tqdm(enumerate(fragments[i], start=0), total=len(fragments), desc="Processing Fragments"):
        target_sequence = fragments[i].upper()
        start_time = time.time()
        chr_number, _, edge_index, matching_paths = process_sequence_matching(target_sequence, sequence, k-1)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(elapsed_time)
        
        edge_labels = chromosomes[list(chromosomes.keys())[chr_number-1]].get_edge_labels()
        
        read_main_nodes = sliding_window(target_sequence, k - 1)
        start_time = time.time()
        result_read, read_edges = process_read_data(read_main_nodes, edge_labels, edge_index)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(elapsed_time)
        readdictt = {i+1: item for i, item in enumerate(read_edges)}
        
        previous_dict, next_dict = find_discontinuous_keys(k, result_read, 100)
        insert_dict = insert_dicts(chr_number, k, readdictt, previous_dict, next_dict)
        print(insert_dict)
        # previous_dict, next_dict = find_repeat_keys(k, result_read, 50)
        # repeat_dict = repeat_dicts(chr_number, k, readdictt, previous_dict, next_dict)
        repeat_dict ={}
        previous_dict, next_dict = find_delete_keys(k, result_read)
        delet_dict = dele_dicts(chr_number, k, readdictt, previous_dict, next_dict)
        
        invert_dict = None
        # previous_dict, next_dict = find_invert_keys(k, result_read)
        
        # if previous_dict:
        #     invert_dict = find_inverted_sequences(k, readdictt, previous_dict, next_dict, chromosomes, chr_number)
        
        tra_dict = None
        # # if list(readdictt.keys())[-1] - list(result_read.keys())[-1] > 2 * k:
        # if readdictt and result_read and (list(readdictt.keys())[-1] - list(result_read.keys())[-1] > 2 * k):

        #     tra_dict = process_translocation(read_edges, result_read, target_sequence, sequence, k, chromosomes, chr_number)
        
        if insert_dict is not None:
            result_dict.update(insert_dict)
        if repeat_dict is not None:
            result_dict.update(repeat_dict)
        if delet_dict is not None:
            result_dict.update(delet_dict)
        if invert_dict is not None:
            result_dict.update(invert_dict)
        if tra_dict is not None:
            result_dict.update(tra_dict)
    
    return result_dict
# 调用函数
if __name__ == "__main__":
    # sequence = [["gtatttcatgATTTGCTCCAtttatgcaatgtatcctagtactagcttACTactccttaaatctatcactatcgaatc"],\
    #     ['axojompmihsuhdlndwogwjonbhfyryuabvcghxgfdyyvyssswwfefefefe']]

    reference_genome_file =r"C:\Users\1\Desktop\code_ppg\output.fasta" # 将文件路径替换为您的参考基因组文件路径
    genome_sequence = [read_reference_genome(reference_genome_file)]
    # reference_genome_file2 = r"C:\Users\1\Desktop\code_ppg\frandom_sequence_2.fasta"
    # genome_sequence1 = [read_reference_genome(reference_genome_file2)]
    sequence= []
    sequence.append(genome_sequence)
    # sequence.append(genome_sequence1)
    # 打印前100个碱基
    # print(res)
    chromosomes = {}
    for chri in range(1,2):
    # 生成对象名称，例如 "chr1", "chr2", ...
        chromosome_name = "chr" + str(chri)
        
        # 创建对象，并将其存储在字典中
        chromosome = DeBruijnGraph()  # 用您的类初始化对象
        k =11
        main_chain_nodes = sliding_window(sequence[chri-1][0], k-1)
        nodes_to_add = main_chain_nodes
        node_iterator = (new_node for new_node in main_chain_nodes)
        chromosome.add_nodes_from(node_iterator)

        for m in range (len(main_chain_nodes)-1):
        # 添加带有标签的边，连接相同内容的不同节点
            chromosome.add_main_edge(main_chain_nodes[m], main_chain_nodes[m+1], label=m+1)
        # chromosome.insert_subgraph(k,1, sequence[chri-1][0], ['aaaa'][0])
        # chromosome.insert_subgraph(k,1, sequence[chri-1][0], ['CTATAGTAGGCCTTAGTGTGGATG'][0])
        for i in [90065056, 90203042, 90328858, 90495928, 90629094, 91580048, 91696541, 91920177, 92350829, 92500019, 92762507, 92902006, 93206786, 93568072, 93917624, 93940849, 94271745, 94671938, 94965974, 95568587, 95676546, 96115250, 96361021, 97082882, 97681657, 97822521, 98453305, 98664839, 99220143, 99276658, 99757679, 99878789]:
            chromosome.insert_subgraph(k,i, sequence[chri-1][0], ['CTAAT'][0])
        # chromosome.delete(k,sequence[chri-1][0],10,21)
       
        chromosomes[chromosome_name] = chromosome
        # chromosome.save_as_gfa(r"C:\Users\1\Desktop\code_ppg\sample.fasta")
    # graph\
   
    from tqdm import tqdm
    import networkx as nx
    import matplotlib.pyplot as plt 
    from Bio import SeqIO 
    #region
    # fasta_filename = r"C:\Users\1\Desktop\code_ppg\fragments.fasta"
    # fragments = []
    # result_dict = {}
    # for record in SeqIO.parse(fasta_filename, "fasta"):
    #     fragments.append(record.seq)
    # for i, fragment in tqdm(enumerate(fragments, start=0), total=len(fragments), desc="Processing Fragments"):
    #     target_sequence = fragments[i].upper()
    #     chr_number,_,edge_index ,matching_paths = process_sequence_matching(target_sequence, sequence, lennode =k-1)
    #     edge_labels = chromosomes[list(chromosomes.keys())[chr_number-1]].get_edge_labels()
    #     read ,_,_,read_main_nodes = create_labeled_multidigraph(target_sequence,k,edge_index, edge_labels)
    #     result_read,read_edges = process_read_data(read_main_nodes, edge_labels, edge_index)
    #     # print( result_read)
    #     readdictt ={i+1: item for i, item in enumerate(read_edges)}
    #     # insert
    #     previous_dict,next_dict = find_discontinuous_keys( k,result_read,100)
    #     insert_dict = insert_dicts(chr_number,k,readdictt,previous_dict,next_dict)
    #     # print(insert_dict)
    #     # repeat
    #     previous_dict,next_dict = find_repeat_keys(k,result_read,50)
    #     repeat_dict = repeat_dicts(chr_number,k,readdictt,previous_dict,next_dict)
    #     # print(repeat_dict)
    #     # delete
    #     previous_dict,next_dict = find_delete_keys(k,result_read)
    #     delet_dict = dele_dicts(chr_number,k,readdictt,previous_dict,next_dict)
    #     # print( delet_dict)
    #     # invert
    #     invert_dict=None
    #     previous_dict,next_dict = find_invert_keys(k,result_read)
    #     # print( previous_dict,next_dict)
    #     if previous_dict:
        
    #         invert_dict =find_inverted_sequences(k, readdictt, previous_dict, next_dict, chromosomes, chr_number)
    #     # print( invert_dict)

    #     # TRA
    #     tra_dict = None
    #     if list(readdictt.keys())[-1] - list(result_read.keys())[-1] > 2 * k:
        
    #         tra_dict = process_translocation(read_edges,result_read,target_sequence, sequence, k, chromosomes, chr_number)
    #         # print(tra_dict)
       
    #     if insert_dict is not None :
    #         result_dict.update(insert_dict)
    #         # print(result_dict)
    #     if repeat_dict is not None :
    #         result_dict.update(repeat_dict)
    #     if delet_dict is not None :
    #         result_dict.update(delet_dict)
    #     if invert_dict is not None :
    #         result_dict.update(invert_dict)
    #     if tra_dict is not None:
    #         result_dict.update(tra_dict)
    # print(result_dict)
#endregion

    fasta_filename = r"C:\Users\1\Desktop\code_ppg\sample.fasta"
    result_dict = process_fasta_file(fasta_filename)

    # 打印或处理result_dict
    print(result_dict)


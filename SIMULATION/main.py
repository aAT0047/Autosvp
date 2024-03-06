import hashlib
import networkx as nx
from usingtools import *
from PAN_GENE_GRAPH import DeBruijnGraph
from NWAligner import *
if __name__ == "__main__":
    de_bruijn_graph =  DeBruijnGraph()
    # sequences = [generate_random_sequence(2000) for _ in range(2000)]
    # print(sequences[0])
    sequences11 ="CCGAGTTGGTGCGGATTATTCAGGAGTCAGGGGCCTATAGTAGGCCTTAGTGTGGATGAGGTGAGGACGCGGAA"
    sequence =  sequences11
    # print(sequence)
    k =7
    main_chain_nodes = sliding_window(sequence, k-1)


    nodes_to_add = main_chain_nodes
    for i in range(len(main_chain_nodes)):
        
       de_bruijn_graph.add_nodes_from(nodes_to_add)
    for i in range (len(main_chain_nodes)-1):
    # 添加带有标签的边，连接相同内容的不同节点
        de_bruijn_graph.add_main_edge(main_chain_nodes[i], main_chain_nodes[i+1], label=i+1)
    # de_bruijn_graph.add_main_edge(main_chain_nodes[-2], main_chain_nodes[-1], label=len(main_chain_nodes))
    # for i in range(1000):
    #     de_bruijn_graph.insert_subgraph(k,random_numbers1[i], sequences11, random_sequences[i])
   
    # for i in range(2,2000,2):
    #     num =sorted( [random_numbers2[i],random_numbers2[i+1]])
    #     de_bruijn_graph.add_missing_edges(k,sequences11,num[0],num[1])




    de_bruijn_graph.insert_subgraph(k,1, sequences11, ['agactatatat'][0])
    # de_bruijn_graph.add_missing_edges(k,sequences11,2,7)
    edge_labels = de_bruijn_graph.get_edge_labels()
    print(edge_labels.items() )
    print("\n边和标签:")
    for (src, dst, key), label in edge_labels.items():
        print("源节点:", src, " 目标节点:", dst, " 边标签:", label)
    de_bruijn_graph.print_node_hash_map()
    # 获取所有节点
    for node, node_hash in de_bruijn_graph.node_hash_map:
        print(f"Node: {node}, Hash: {node_hash}")


    # target_sequence = "ACTAGTATTTCTG.TTTCAAATACTTTTTGGCTCAGCGGAAAGAC"
    # k =31
    # # 在图中查找目标序列的位置

    # matching_paths= de_bruijn_graph.find_sequence_location(target_sequence, lennode=k-1)
    # print(matching_paths)
    # # de_bruijn_graph.has_edge('GGCG','GCGA')
    # # print(all_matching_paths)
    # # 获得路径
    # edge_label_list=get_last_numbers(get_map_edge_labels(edge_labels.items(),matching_paths[0]))
    # print( edge_label_list)
    # max_keys, max_value = count_elements( edge_label_list)
    # edge_index =find_longest_consecutive_sequences (  edge_label_list)
    # if max_value <= len(edge_index[0]):
    #     edgepath = extract_sublists_by_indices( get_map_edge_labels(edge_labels.items(),matching_paths[0]), edge_index[0])
    #     chr_start_index= find_substring_indices(target_sequence,edgepath)
    # else:
    #     edgepath = extract_sublists_by_indices( get_map_edge_labels(edge_labels.items(),matching_paths[0]), [max_keys])
    #     chr_start_index= find_indel_indices(target_sequence,edgepath)

    # print(edgepath) 
    # print(chr_start_index) 
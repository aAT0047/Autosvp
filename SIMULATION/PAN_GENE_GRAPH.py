import networkx as nx
import matplotlib.pyplot as plt
import random
import json
import hashlib
from usingtools import *
class DeBruijnGraph:
    def __init__(self):
         self.graph = nx.MultiDiGraph()
         self.node_hash_map = {}
    def print_node_hash_map(self):
        print(self.node_hash_map)
    def add_invert_edge(self, node, neighbor,label):  # Add label parameter
        self.graph.add_edge(node, neighbor, key=None, label=label)
    def add_delete_edge(self, node, neighbor,label):  # Add label parameter
        self.graph.add_edge(node, neighbor, key=None, label=label)    
    def add_main_edge(self, node, neighbor, label):  # Add label parameter
        self.graph.add_edge(node, neighbor, key=None, label=label)
        
    def add_label(self, node,label):
        if node in self.graph:
            self.graph[node].append(label)
        else:
            self.graph[node] = [label]

    def visualize_graph(self):
            # 使用 networkx 创建图
            G = nx.DiGraph()

            # 添加图中的节点和边
            for node, neighbors in self.graph.items():
                G.add_node(node)
                for neighbor in neighbors:
                    G.add_edge(node, neighbor)

            # 绘制图形并使用 graphviz_layout 进行布局
            pos = nx.drawing.nx_pydot.graphviz_layout(G, prog='dot')
            nx.draw(G, pos, with_labels=True, node_size=1000, node_color='skyblue', font_size=8)
            plt.show()
    def get_nodes(self):
        return list(self.graph.nodes())

    def get_edge_labels(self):
        edge_labels = nx.get_edge_attributes(self.graph, 'label')
        return edge_labels

    def get_edges(self):
        return list(self.graph.edges())
        return edges
    def insert_subgraph(self,k,insert_index,seq,insert_seq):
        if len(seq[0:insert_index])<k:
            sequence = seq[0:insert_index]+insert_seq+seq[insert_index:insert_index+k]
        elif len(seq[insert_index::])<k:
            sequence = seq[insert_index-k+1:insert_index]+insert_seq+seq[insert_index::]
            # print(seq[insert_index+1+k::])
        else:
            # 为了防止存在超过1个位置的节点相似，在3k之外开始插入
            
            
           
            sequence = seq[insert_index-k+1:insert_index]+insert_seq+seq[insert_index:insert_index+k-1]
            
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            nodes=[kmer[:-1],kmer[1:]]
            self.add_nodes_from( nodes)
            
            self.add_invert_edge(kmer[:-1], kmer[1:],label =f"invert_{insert_index}")
        # return  de_bruijn_graph

    def delete(self,k,seq,delete_start_index,delete_end_index):
        if len(seq[0:delete_start_index])<k:
            sequence = seq[0:delete_start_index]+seq[delete_end_index-1:delete_end_index+k-1]
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                # print(kmer)
                self.add_delete_edge(kmer[:-1], kmer[1:],label =f"delete_{delete_start_index}_{delete_end_index}")
        elif len(seq[delete_end_index-1::])<k:
            sequence = seq[delete_start_index-k:delete_start_index]+seq[delete_end_index-1::]
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                # print(kmer)
                self.add_delete_edge(kmer[:-1], kmer[1:],label =f"delete_{delete_start_index}_{delete_end_index}")
            # print(seq[insert_index+1+k::])
        else:
           
            sequence = seq[delete_start_index-k:delete_start_index]+seq[delete_end_index-1:delete_end_index+k-1]
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                # print(kmer)
                self.add_delete_edge(kmer[:-1], kmer[1:],label =f"delete_{delete_start_index}_{delete_end_index}")
       
    def count_nodes(self):
        return len(self.graph)
    
    def to_dict(self):
        return self.graph

    def save_as_gfa(self, filename):
        with open(filename, "w") as f:
            f.write("H\tVN:Z:1.0\n")

            for node in self.graph.nodes():
                label = self.graph.nodes[node].get('label', "")
                f.write(f"S\t{node}\t{label}\n")

            for u, v, label in self.graph.edges(data='label'):
                f.write(f"L\t{u}\t+\t{v}\t+\t{label}\n")

    def save_as_json(self, filename):
        data = {
            "nodes": list(self.graph.keys()),
            "edges": [(node, neighbor) for node, neighbors in self.graph.items() for neighbor in neighbors]
        }
        with open(filename, "w") as f:
            json.dump(data, f, indent=4)

    
    def generate_hash(self, node):
        hasher = hashlib.sha256()
        node = str(node) 
        hasher.update(node.encode('utf-8'))
        return hasher.hexdigest()

    def build_node_hash_map(self):
        for node in self.graph:
            node_hash = self.generate_hash(node)
            self.node_hash_map[node_hash] = node
            # print(self.node_hash_map)

    def has_edge(self, node1, node2):
        return self.graph.has_edge(node1, node2)

    def add_nodes_from(self, nodes):
        self.graph.add_nodes_from(nodes)

    def find_sequence_location(self, target_sequence, lennode):
        matching_nodes = []
        unmatching_nodes = []
        self.build_node_hash_map()

        prev_matching_node = None  # 用于跟踪最近的匹配节点

        for i in range(len(target_sequence) - lennode + 1):
            segment = target_sequence[i:i + lennode]
            segment_hash = self.generate_hash(segment)

            if segment_hash in self.node_hash_map:
                node = self.node_hash_map[segment_hash]
                matching_nodes.append(node)
                
                if prev_matching_node is not None:  # 如果存在最近的匹配节点
                    if self.has_edge(prev_matching_node, node):
                        prev_matching_node = node
                    else:
                        matching_nodes.pop()  # 移除刚刚添加的节点
                        break
                else:
                    prev_matching_node = node
            else:
                remaining_sequence = segment
                unmatching_nodes.append(segment)
                prev_matching_node = None  # 重置最近的匹配节点

        return matching_nodes, unmatching_nodes

    def find_multiple_sequences(self, target_sequence, lennode):
        all_matching_paths = []

        self.build_node_hash_map()

        for i in range(len(target_sequence) - lennode + 1):
            matching_nodes = []
            unmatching_nodes = []
            self.build_node_hash_map()

            prev_matching_node = None  # 用于跟踪最近的匹配节点

            for i in range(len(target_sequence) - lennode + 1):
                segment = target_sequence[i:i + lennode]
                segment_hash = self.generate_hash(segment)

                if segment_hash in self.node_hash_map:
                    node = self.node_hash_map[segment_hash]
                    matching_nodes.append(node)

                    if prev_matching_node is not None:  # 如果存在最近的匹配节点
                        if self.has_edge(prev_matching_node, node):
                            prev_matching_node = node
                        else:
                            matching_nodes.pop()  # 移除刚刚添加的节点
                            break
                    else:
                        prev_matching_node = node
                else:
                    remaining_sequence = segment
                    unmatching_nodes.append(segment)
                    prev_matching_node = None  # 重置最近的匹配节点

            all_matching_paths.append((matching_nodes, unmatching_nodes))

        # 从所有路径中选择匹配节点最多的路径
        max_matching_nodes = 0
        best_matching_path = None

        for matching_nodes, _ in all_matching_paths:
            if len(matching_nodes) > max_matching_nodes:
                max_matching_nodes = len(matching_nodes)
                best_matching_path = matching_nodes

        return best_matching_path, all_matching_paths
    
    def find_paths_with_label_to_target(self, label_to_find, target_node):
        paths = []
        edge_labels = nx.get_edge_attributes(self.graph, 'label')

        for edge, edge_label in edge_labels.items():
            if edge_label == label_to_find and edge[0] == target_node:
                self._find_paths_recursive(paths, edge, edge_labels, target_node)

        return paths

    def _find_paths_recursive(self, paths, current_edge, edge_labels, target_node, current_path=[]):
        current_path = current_path + [current_edge[0]]

        if current_edge[1] == target_node:
            paths.append(current_path + [target_node])
            return

        successors = self.graph.successors(current_edge[1])

        for next_node in successors:
            next_edge = (current_edge[1], next_node)
            if next_edge in edge_labels and edge_labels[next_edge] == 1:
                self._find_paths_recursive(paths, next_edge, edge_labels, target_node, current_path)



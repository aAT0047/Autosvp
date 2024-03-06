# class NWAligner:
#     def __init__(self, maxIndel, M, I, G):
#         # 初始化NWAligner类，设置最大插入/删除限制（maxIndel）和得分参数（M、I、G）
#         # 初始化一些属性，包括currSize、maxIndel、M、I、G和matrix

#     def reserveBanded(self, l1, l2):
#         # 根据输入的s1和s2的长度，为对角线对齐算法保留带状矩阵的内存空间

#     def S(self, a, b):
#         # 计算两个字符a和b之间的得分，用于匹配或替代

#     def alnGlobFreeEndGap(self, s1, s2):
#         # 在全局自由端缺口对齐中对两个序列s1和s2进行对齐
#         # 移除尾部的缺口并返回对齐结果

#     def align(self, s1, s2):
#         # 在全局端到端对齐中对两个序列s1和s2进行对齐
#         # 返回对齐的结果，包括对齐长度和得分

#     def trimTrailingGaps(self, alnRes):
#         # 去除对齐结果中的尾部缺口

#     def alignBanded(self, s1, s2):
#         # 在带状对齐中对两个序列s1和s2进行对齐

# class AlnRes:
#     def __init__(self, s1len, s2len, score):
#         # 初始化AlnRes类，表示对齐的结果，包括对齐长度和得分
import numpy as np

def smith_waterman(seq1, seq2, match=2, mismatch=-2, gap=-1):
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

if __name__ == "__main__":
    seq1 = "AGTACGTACGTACTCTACC"
    seq2 = "ATACGTACGTACTCTACC"

    aligned_seq1, aligned_seq2, score = smith_waterman(seq1, seq2)

    print("比对后序列 1:", aligned_seq1)
    print("比对后序列 2:", aligned_seq2)
    print("比对得分:", score)

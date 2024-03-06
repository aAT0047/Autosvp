import torch
import torch.nn.functional as F
# from sklearn.datasets import make_multilabel_classification
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv
from torch_geometric.utils import from_scipy_sparse_matrix
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from torch_geometric.utils import from_scipy_sparse_matrix



# 设置随机种子以保证结果可重复
np.random.seed(42)

# 定义数据数量和特征数量
n_samples = 1000
n_features = 10
n_targets = 3

# 模拟特征数据
X = np.random.randn(n_samples, n_features)

# 模拟目标数据（这里我们简单地使每个目标成为输入特征的线性组合，加上一些噪声）
weights = np.array([[1.5, 2.0, 0.5],
                    [-0.5, 1.0, 2.0],
                    [2.0, -1.0, 1.5],
                    [0.5, 1.5, -1.0],
                    [-1.0, 2.5, 0.5],
                    [1.0, -0.5, 2.0],
                    [0.5, 1.0, -2.0],
                    [2.0, -1.0, 0.5],
                    [-2.0, 0.5, 1.0],
                    [1.0, -2.0, -0.5]])
Y = X.dot(weights) + np.random.normal(0, 0.1, size=(n_samples, n_targets))

# 将数据转换为DataFrame格式（仅为了展示，下面的模型部分不需要这个格式）
df_features = pd.DataFrame(X, columns=[f'feature_{i}' for i in range(n_features)])
df_targets = pd.DataFrame(Y, columns=[f'target_{i}' for i in range(n_targets)])
df = pd.concat([df_features, df_targets], axis=1)

# print(df.head())
# 将df拆分为特征和目标
X = df.drop(columns=['target_0', 'target_1', 'target_2'])  # 假设目标列名为target_0, target_1, target_2
y = df[['target_0', 'target_1', 'target_2']]

from sklearn.neighbors import kneighbors_graph

# 使用KNN构建邻接矩阵
A = kneighbors_graph(X, n_neighbors=5, mode='connectivity')  # 例如, 选择5个最近邻
A = A.toarray()

# 转换为PyTorch张量
A = torch.tensor(A, dtype=torch.float)

# 标准化数据
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
# 在拆分数据集之前构建edge_index
X_tensor = torch.tensor(X_scaled, dtype=torch.float)
y_tensor = torch.tensor(y.values, dtype=torch.float)
A_sparse = csr_matrix(A)
edge_index, _ = from_scipy_sparse_matrix(A_sparse)
# print(edge_index)
# 使用train_test_split按节点划分数据集
train_indices, test_indices, y_train, y_test = train_test_split(
    range(len(X_tensor)), y_tensor, test_size=0.2, random_state=42
)

train_mask = torch.zeros(len(X_tensor), dtype=torch.bool)
train_mask[train_indices] = True

test_mask = torch.zeros(len(X_tensor), dtype=torch.bool)
test_mask[test_indices] = True

# 至此，我们已经为每个节点设置了一个mask，指示它是训练节点还是测试节点。
# 在训练过程中，您可以只更新和计算train_mask=True的节点的损失。
# 在测试过程中，可以单独为每个test_mask=True的节点进行预测。

# 示范：定义一个简单的GCN模型
class GCN(torch.nn.Module):
    def __init__(self, in_channels, out_channels):
        super(GCN, self).__init__()
        self.conv1 = GCNConv(in_channels, 64)
        self.conv2 = GCNConv(64, out_channels)

    def forward(self, x, edge_index):
        x = F.relu(self.conv1(x, edge_index))
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.conv2(x, edge_index)
        return x

# 初始化模型和优化器
model = GCN(X_tensor.size(1), y_tensor.size(1))
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

# 训练
model.train()
for epoch in range(500):
    optimizer.zero_grad()
    out = model(X_tensor, edge_index)
    loss = F.mse_loss(out[train_mask], y_train)
    loss.backward()
    optimizer.step()
    if epoch % 10 == 0:
        print(f"Epoch: {epoch}, Loss: {loss.item()}")

# 测试
model.eval()
with torch.no_grad():
    predictions = model(X_tensor, edge_index)[test_mask]
    test_loss = F.mse_loss(predictions, y_test)
    print(f"Test Loss: {test_loss.item()}")
    print (predictions)
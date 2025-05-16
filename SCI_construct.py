import pandas as pd
import numpy as np
import math

def compute_SCI(Cij,disease_patient_num,N,dic_cols,prevalence):
    """
        计算 Salton Cosine Index（SCI）并基于 phi 显著性筛选边。

        参数:
            Cij: 共病矩阵 (2D numpy array, 上三角)
            disease_patient_num: 每个疾病的患者数（DataFrame）
            N: 总样本人群数
            dic_cols: 疾病 id 对应名称的字典
            prevalence: 各疾病流行率

        返回:
            edge_list_new: 筛选后的 SCI 边表（DataFrame）
        """

    phi_node_set=set()
    phi_edge_list=[]
    SCI_edge_list = []

    #1. 对每个共病对，统计具有显著性意义的φ边
    for i in range(Cij.shape[0]):
        for j in range(i + 1, Cij.shape[0]):
            Nij = Cij[i][j]
            if Nij == 0:
                continue

            N1 = disease_patient_num.loc[i, 'patient_number']
            N2 = disease_patient_num.loc[j, 'patient_number']

            # ---------- φ检验 ----------
            a = (0.1 * N1) * (0.1 * N2) * (0.1 * (N - N1)) * (0.1 * (N - N2))
            if (a) <= 0:
                raise Exception("phi 溢出", a)
            phi_ij = ((0.1 * Cij[i][j]) * (0.1 * N) - (0.1 * N1) * (0.1 * N2)) / np.sqrt(a)
            t = 0

            n = 0
            if abs(phi_ij) < 1:  # phi=1时，会发生除零错误,|phi|>1时，会发生计算错误
                n = max(N1, N2)
                # n=N     # 注意测试一下
                t = (phi_ij * math.sqrt(n - 2)) / np.sqrt(1 - (phi_ij ** 2))
            elif phi_ij > 1 or phi_ij < -1:  # 不会大于1
                print("有phi大于1 或者小于-1 ，考虑截断,phi值为：", phi_ij)
                # 若phi=1，只能是这种情况：A病和B病必定同时出现，且A病和B病不单独出现，这时的phi=1；因为前面步骤去除了流行度小于1%的疾病，所以这种情况基本不会发生吧
                t = 0
            else:
                t = 2.77
                n = max(N1, N2)
                raise Exception("有phi等于-1、1 ，n = max(prevalence1, prevalence2)值为：", n)
            if ((n > 1000 and phi_ij > 0 and t >= 2.58) or (n > 500 and phi_ij > 0 and t >= 2.59) or (
                    n > 200 and phi_ij > 0 and t >= 2.60) or (n > 90 and phi_ij > 0 and t >= 2.63) or (
                    n > 80 and phi_ij > 0 and t >= 2.64) or (n > 70 and phi_ij > 0 and t >= 2.65) or (
                    n > 60 and phi_ij > 0 and t >= 2.66) or (n > 50 and phi_ij > 0 and t >= 2.68) or (
                    n > 40 and phi_ij > 0 and t >= 2.70) or (n > 38 and phi_ij > 0 and t >= 2.71) or (
                    n > 35 and phi_ij > 0 and t >= 2.72) or (n > 33 and phi_ij > 0 and t >= 2.73) or (
                    n > 31 and phi_ij > 0 and t >= 2.74) or (n > 30 and phi_ij > 0 and t >= 2.75) or (
                    n > 28 and phi_ij > 0 and t >= 2.76) or (
                    n > 27 and phi_ij > 0 and t >= 2.77)):  # 这里只考虑了两个节点联系比随机情况下更强的情况
                phi_edge_list.append([i, j, N1, N2, phi_ij, t, -999, Cij[i][j]])
                phi_node_set.update([i, j])

            # ---------- 计算 SCI ----------
            try:
                SCI_ij = (0.1 * Cij[i][j]) / np.sqrt(0.1 * N1 * 0.1 * N2)
                SCI_edge_list.append([i, j, N1, N2, SCI_ij, -999, Cij[i][j]])
            except ZeroDivisionError:
                continue

    #2. 计算nab_minimun
    sum=0

    phi_node_list = list(phi_node_set)
    if len(phi_node_list) < 2:
        print("phi_node_list 长度小于2，跳过筛选，返回0条边")
        q = 0
    else:
        total_cooccurrence = sum(Cij[i][j] for i in range(len(phi_node_list)) for j in range(i + 1, len(phi_node_list)))
        total_pairs = len(phi_node_list) * (len(phi_node_list) - 1) / 2
        nab_minimum = total_cooccurrence * 2 / total_pairs if total_pairs > 0 else 0

        #3. 统计满足 nab > nab_min 的边数 q
        q = sum(1 for i in range(Cij.shape[0]) for j in range(i + 1, Cij.shape[1]) if Cij[i][j] > nab_minimum)

    #4. 筛选q个有意义的SCI边
    SCI_edge_list = pd.DataFrame(SCI_edge_list,columns=['node1', 'node2', 'prevalence1', 'prevalence2', 'SCI', '无意义位', 'Nij'])
    edge_list_new = SCI_edge_list.sort_values(by='SCI', ascending=False).iloc[0:q]

    return edge_list_new

if __name__ == '__main__':
    N = 10000  # 总人数
    disease_names = ['A', 'B', 'C', 'D']
    disease_patient_num = pd.DataFrame({
        'disease_id': list(range(4)),
        'patient_number': [1200, 1500, 800, 1000]
    })
    prevalence = disease_patient_num['patient_number'] / N
    dic_cols = {'disease_name': {0: 'A', 1: 'B', 2: 'C', 3: 'D'}}

    Cij = np.array([
        [0, 100, 50, 80],
        [0, 0, 60, 70],
        [0, 0, 0, 30],
        [0, 0, 0, 0]
    ])

    #共病网络边表
    edges = compute_SCI(Cij, disease_patient_num, N, dic_cols, prevalence)
    print(edges)



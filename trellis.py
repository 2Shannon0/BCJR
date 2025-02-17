from util import *
import argparse
import GFn
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-input_matrix')
args = parser.parse_args()

def symbol_all(nbit=1, b=1):
    ret_list = []
    if nbit == 1:
        for a0 in range(0, 2**b):
            ret_list.append(GFn.GFn([a0], 1))
        return ret_list
    if nbit == 2:
        for a0 in range(0, 2):
            for a1 in range(0, 2**b):
                ret_list.append(GFn.GFn([a0, a1], 2))
        return ret_list
    raise ValueError(f"No symbol for nbit = {nbit}")

# Функция построения и отображения графа (решетки) по слоям
def plot_sections(vex, edg, position):
    start, end = position # Определяем границы слоев для отрисовки
    edges = []  # Список рёбер
    pos_dict = {} # Словарь с координатами узлов

    # Обход вершин по слоям
    for num_layer, v_layer in zip(range(start, end+1), vex[start:end+1]):
        for v in v_layer:
            v_name = v2str(v, num_layer)
            pos_dict[v_name] = np.array([num_layer, -arr2int(v)])
    
    # Обход рёбер по слоям
    for num_layer, e_layer in zip(range(start, end), edg[start:end]):
        for e in e_layer:
            v0, a, v1 = e
            edges.append((v2str(v0, num_layer), v2str(v1, num_layer+1), {'weight': int(a)}))

    G = nx.Graph()
    G.add_edges_from(edges)
    
    pos = nx.spring_layout(G) # Вычисление позиций узлов
    nx.draw_networkx_nodes(G, pos_dict, node_size=1200) # Рисуем узлы
    nx.draw_networkx_labels(G, pos_dict, font_size=12, font_family='sans-serif') # Подписываем узлы
    nx.draw_networkx_edges(G, pos_dict, edgelist=edges, width=6) # Отображаем рёбра
    labels = nx.get_edge_attributes(G, 'weight') # Получаем веса рёбер
    nx.draw_networkx_edge_labels(G, pos_dict, edge_labels=labels) # Добавляем подписи рёбер
    plt.show()

# Функция удаления рёбер, ведущих в ненулевые состояния (оптимизация решетки) OLD version
def remove_nonzero(edg, vex):
    edg.append([(vex[0][0], 0, vex[0][0])])
    for i in reversed(range(len(edg)-1)):
        vex_next = vex_connected(vex[i+1], edg[i+1])
        edg[i] = [edge for edge in edg[i] if in_list(vex_next, edge[2])]

p_mat = read_mat(args.input_matrix)

# Определяем начальное (нулевое) состояние как массив нулей
zero_state = np.array([GFn.GFn(0, 1)] * int(p_mat.shape[1]))
vex = [[zero_state]] # Начальная вершина решетки
edg = [] # Список рёбер

# Создаём массив символов на основе входной матрицы
symbol_np_arr = np.empty([int(p_mat.shape[0]), int(p_mat.shape[1])], dtype=GFn.GFn)
for x in range(p_mat.shape[0]):
    for y in range(p_mat.shape[1]):
        symbol_np_arr[x][y] = GFn.GFn(p_mat[x][y], 1)
print('symbol_np_arr:\n')
print_matrix(symbol_np_arr)

# Построение решетки
for layer in range(p_mat.shape[0]):
    vex_new = []
    edg_new = []
    symbol_layer = symbol_np_arr[layer]
    for v_last in vex[-1]:
        for symbol in symbol_all():
            add_v = symbol * symbol_layer + v_last
            edge = (v_last, symbol, add_v)
            edg_new.append(edge)
            if not in_list(vex_new, add_v):
                vex_new.append(add_v)
    vex.append(vex_new)
    edg.append(edg_new)

# Отобразить структуру
# for layer, states in enumerate(vex):
#     print(f"Слой {layer}:")
#     for state in states:
#         print(f"  {state}")

# for layer, edges in enumerate(edg):
#     print(f"Edges at Layer {layer}:")
#     for v_last, symbol, add_v in edges:
#         print(f"  {v_last} --({symbol})--> {add_v}")


remove_nonzero(edg, vex)
plot_sections(vex, edg, [0, p_mat.shape[0]])

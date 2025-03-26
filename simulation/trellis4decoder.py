from util import *
import numpy as np
import GFn
import networkx as nx
import matplotlib.pyplot as plt

class Trellis:
    def __init__(self, matrix_file = None, vex = None, edg = None):
        self.file = matrix_file
        self.p_mat = matrix_file and read_mat(matrix_file)
        self.vex = vex
        self.edg = edg

    def symbol_all(self, nbit=1, b=1):
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
    
    def remove_nonzero(self, edg, vex, vex_set):
        edg.append([(vex[0][0], 0, vex[0][0])])
        print("removing nonzero")
        for i in reversed(range(len(edg)-1)):
            print(f"process {len(edg) - 1 - i} / {len(edg) - 1}")
            vex_next = vex_connected(vex_set[i+1], edg[i+1])
            edg[i] = [edge for edge in edg[i] if gfn_array_to_str(edge[2]) in vex_next]

    def build_trellis(self):
        # p_mat = read_mat(self.file)

        # Определяем начальное (нулевое) состояние как массив нулей
        zero_state = np.array([GFn.GFn(0, 1)] * int(self.p_mat.shape[1]))
        self.vex = [[zero_state]] # Начальная вершина решетки
        vex_set = [set(gfn_array_to_str(zero_state))] # Начальная вершина решетки
        self.edg = [] # Список рёбер

        # Создаём массив символов на основе входной матрицы
        symbol_np_arr = np.empty([int(self.p_mat.shape[0]), int(self.p_mat.shape[1])], dtype=GFn.GFn)
        for x in range(self.p_mat.shape[0]):
            for y in range(self.p_mat.shape[1]):
                symbol_np_arr[x][y] = GFn.GFn(self.p_mat[x][y], 1)

        symbol_all = self.symbol_all()

        # Построение решетки
        print("building trellis")
        for layer in range(self.p_mat.shape[0]):
            vex_new = []
            vex_new_set = set()
            edg_new = []
            symbol_layer = symbol_np_arr[layer]
            print(f"process {layer + 1} / {self.p_mat.shape[0]}. vex to pass {len(self.vex[-1])}")
            for v_last in self.vex[-1]:
                for symbol in symbol_all:
                    add_v = symbol * symbol_layer + v_last
                    edge = (v_last, symbol, add_v)
                    edg_new.append(edge)

                    add_v_str = gfn_array_to_str(add_v)
                    if add_v_str not in vex_new_set:
                        vex_new.append(add_v)
                        vex_new_set.add(add_v_str)
            vex_set.append(vex_new_set)
            self.vex.append(vex_new)
            self.edg.append(edg_new)

        self.remove_nonzero(self.edg, self.vex, vex_set)
    
    # Функция построения и отображения графа (решетки) по слоям
    def plot_sections(self, title=None, save_name=None):
        # Определяем границы слоев для отрисовки
        start, end = [0, len(self.vex)-1]

        edges = []  # Список рёбер
        pos_dict = {} # Словарь с координатами узлов

        # Обход вершин по слоям
        for num_layer, v_layer in zip(range(start, end+1), self.vex[start:end+1]):
            for v in v_layer:
                v_name = v2str(v, num_layer)
                pos_dict[v_name] = np.array([num_layer, -arr2int(v)])
        
        # Обход рёбер по слоям
        for num_layer, e_layer in zip(range(start, end), self.edg[start:end]):
            for e in e_layer:
                v0, a, v1 = e
                edges.append((v2str(v0, num_layer), v2str(v1, num_layer+1), {'weight': int(a)}))
        plt.figure(figsize=(16, 9))
        G = nx.Graph()
        G.add_edges_from(edges)
        
        pos = nx.spring_layout(G) # Вычисление позиций узлов
        nx.draw_networkx_nodes(G, pos_dict, node_size=1200) # Рисуем узлы
        nx.draw_networkx_labels(G, pos_dict, font_size=12, font_family='sans-serif') # Подписываем узлы
        nx.draw_networkx_edges(G, pos_dict, edgelist=edges, width=6) # Отображаем рёбра
        labels = nx.get_edge_attributes(G, 'weight') # Получаем веса рёбер
        nx.draw_networkx_edge_labels(G, pos_dict, edge_labels=labels) # Добавляем подписи рёбер

        if title is not None:
            plt.title(title)

        if save_name is None:
            plt.show()
        else: 
            plt.savefig(f'process_in_pictures/{save_name}.png', dpi=300, bbox_inches='tight')
    
    def plot_sections_float(self, title=None, save_name=None):
        # Определяем границы слоев для отрисовки
        start, end = [0, len(self.vex)-1]

        edges = []  # Список рёбер
        pos_dict = {} # Словарь с координатами узлов

        # Обход вершин по слоям
        for num_layer, v_layer in zip(range(start, end+1), self.vex[start:end+1]):
            for v in v_layer:
                v_name = v2str(v, num_layer)
                pos_dict[v_name] = np.array([num_layer, -arr2int(v)])
        
        # Обход рёбер по слоям
        for num_layer, e_layer in zip(range(start, end), self.edg[start:end]):
            for e in e_layer:
                v0, a, v1 = e
                edges.append((v2str(v0, num_layer), v2str(v1, num_layer+1), {'weight': round(a, 3)}))
        plt.figure(figsize=(16, 9))
        G = nx.Graph()
        G.add_edges_from(edges)

        pos = nx.spring_layout(G) # Вычисление позиций узлов
        nx.draw_networkx_nodes(G, pos_dict, node_size=1200) # Рисуем узлы
        nx.draw_networkx_labels(G, pos_dict, font_size=12, font_family='sans-serif') # Подписываем узлы
        nx.draw_networkx_edges(G, pos_dict, edgelist=edges, width=6) # Отображаем рёбра
        labels = nx.get_edge_attributes(G, 'weight') # Получаем веса рёбер
        nx.draw_networkx_edge_labels(G, pos_dict, edge_labels=labels) # Добавляем подписи рёбер
        if title is not None:
            plt.title(title)

        if save_name is None:
            plt.show()
        else:
            plt.savefig(f'process_in_pictures/{save_name}.png', dpi=300, bbox_inches='tight')

    def plot_sections_float_alpha(self, title=None, save_name=None):
        start, end = [0, len(self.vex) - 1]
        edges = []  # Список рёбер
        pos_dict = {}  # Словарь с координатами узлов
        node_values = {}  # Подписи узлов
        # Обход вершин по слоям
        for num_layer, v_layer in zip(range(start, end + 1), self.vex[start:end + 1]):
            for key in v_layer:
                node_name = f"{num_layer}_{key}"  # Приводим к формату v2str
                pos_dict[node_name] = np.array([num_layer, -arr2int(key)])  # Позиция узла
                node_values[node_name] = f"{key}\n{round(v_layer[key], 3)}"
        
        # Обход рёбер по слоям
        for num_layer, e_layer in zip(range(start, end), self.edg[start:end]):
            for e in e_layer:
                v0, a, v1 = e
                edges.append((v2str(v0, num_layer), v2str(v1, num_layer+1), {'weight': round(a, 3)}))
        plt.figure(figsize=(16, 9))
        G = nx.Graph()
        G.add_edges_from(edges)
        
        # Вычисление позиций узлов
        nx.draw_networkx_nodes(G, pos_dict, node_size=1200)
        nx.draw_networkx_labels(
            G, pos_dict,
            font_size=14,
            font_family='sans-serif',
            labels=node_values,
            # verticalalignment='bottom',
            font_color="red"
        )
        nx.draw_networkx_edges(G, pos_dict, edgelist=edges, width=6)
        labels = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edge_labels(G, pos_dict, edge_labels=labels)
        
        if title is not None:
            plt.title(title)
        
        if save_name is None:
            plt.show()
        else:
            plt.savefig(f'process_in_pictures/{save_name}.png', dpi=300, bbox_inches='tight')
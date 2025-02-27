from copy import deepcopy
import numpy as np
from scipy.special import logsumexp
from trellis4decoder import Trellis

def gfn_array_to_str(gfns: list) -> str:
    result = ""
    for gfn in gfns:
        result += str(gfn)
    return result

class BCJRDecoder:
    def __init__(self, vex, edg, H, sigma=1.0):
        self.vex = vex
        self.edg = edg
        self.H = H
        self.sigma = sigma
        self.snr = 1 / (2 * sigma**2)

    def decode(self, llr, sigma):

        a_priori = 0.5
        dispersion = sigma**2

        Trellis(vex = self.vex, edg = self.edg).plot_sections('Исходная решетка', '1.Исходная решетка')

        self.make_edges_with_bpsk()

        Trellis(vex = self.vex, edg = self.edg).plot_sections('Решетка после BPSK', '2.Решетка после BPSK')

        # Берем структуру edg, меняя значения в ребре на гамму
        gamma = deepcopy(self.edg)
        # последнее ребро не нужно, ибо оно по факту ни о чем. Кол-во ярусов должно быть N, а в массиве их N+1
        del gamma[-1]

        print("gamma")
        constant_coef = a_priori * (1 / np.sqrt(2 * np.pi * dispersion)) # ВОПРОС. В матлабе эта величина еще возводится в степень N
        for i in range(len(gamma)):
            for j in range(len(gamma[i])):
                diff = (llr[i] - self.edg[i][j][1])**2 / (2 * dispersion)
                gamma[i][j] = (
                    gfn_array_to_str(gamma[i][j][0]),
                    constant_coef * np.exp(-diff),
                    gfn_array_to_str(gamma[i][j][2])
                )
                print(gamma[i][j])
            print()
        Trellis(vex = self.vex, edg = gamma).plot_sections_float('Поменяли значения в ребрах на расчитанные gamma', '3.Gamma')

        print("alpha")
        # Сохраняем альфы в виде списка из словарей. Каждый i-ый словарь - это, по сути, ярус
        # Ключ словаря - синдром, он на каждом ярусе уникален. Значение - посчитанная альфа
        # Такой подход обусловлен тем, что пробегаясь по ребрам, мы могли быстро извлечь альфу следующей вершины
        alpha = [{} for _ in range(len(gamma) + 1)]
        alpha[0][gamma[0][0][0]] = 1
        print(alpha[0])
        for i in range(len(gamma)):
            for j in range(len(gamma[i])):
                cur_gamma = gamma[i][j][1]
                next_vex = gamma[i][j][2]
                prev_vex = gamma[i][j][0]

                new_alpha = cur_gamma * alpha[i][prev_vex]
                # Если вершина с синдромом уже была рассмотрена, то считаем альфу
                # как сумму раннее посчитанного значения и текущего.
                # Это случай когда в вершину введет 2 ребра
                if next_vex in alpha[i + 1]:
                    alpha[i + 1][next_vex] += new_alpha
                # В противном случае просто записываем значение
                else:
                    alpha[i + 1][next_vex] = new_alpha

            # Normalization
            summa = sum(alpha[i + 1].values())
            for key in alpha[i + 1].keys():
                alpha[i + 1][key] /= summa

            print(alpha[i + 1])
        print()
        Trellis(vex = alpha, edg = gamma).plot_sections_float('Поменяли значения в вершинах на расчитанные alpha', '4.Alpha')

        print("beta")
        beta = [{} for _ in range(len(gamma) + 1)]
        beta[-1][gamma[0][0][0]] = 1 # в матлабе последний бета равен последнему альфа, а не харкод 1. Однако, в таком подходе в результате нормализации послдений альфа всегда будет 1
        print(beta[-1])
        for i in range(len(gamma) - 1, -1, -1):
            for j in range(len(gamma[i])):
                cur_gamma = gamma[i][j][1]
                # обратный проход - next и prev меняются местами
                next_vex = gamma[i][j][0]
                prev_vex = gamma[i][j][2]

                new_beta = cur_gamma * beta[i + 1][prev_vex]
                if next_vex in beta[i]:
                    beta[i][next_vex] += new_beta
                else:
                    beta[i][next_vex] = new_beta

            # Normalization
            summa = sum(beta[i].values())
            for key in beta[i].keys():
                beta[i][key] /= summa

            print(beta[i])
        print()
        Trellis(vex = beta, edg = gamma).plot_sections_float('Поменяли значения в вершинах на расчитанные beta', '5.Beta')

        print("sigma")
        sigma = deepcopy(gamma)
        for i in range(len(gamma)):
            for j in range(len(gamma[i])):
                cur_gamma = gamma[i][j][1]
                next_vex = gamma[i][j][2]
                prev_vex = gamma[i][j][0]

                # Берем альфу, которая слева от текущего ребра
                cur_alpha = alpha[i][prev_vex]
                # Берем бету, которая справа от текущего ребра
                cur_beta = beta[i + 1][next_vex]

                sigma[i][j] = (
                    prev_vex,
                    cur_gamma * cur_alpha * cur_beta,
                    next_vex
                )

                print(sigma[i][j])
            print()

        print("LLR")
        out_llr = []
        for i in range(len(sigma)):
            up, down = 0, 0
            for j in range(len(sigma[i])):
                # не забываем про БПСК, по факту ищем ребра с 0
                if self.edg[i][j][1] == 1:
                    up += sigma[i][j][1]
                else:
                    down += sigma[i][j][1]

            # такое вообще возможно?
            if down != 0:
                out_llr.append(np.log(up / down))
            else:
                out_llr.append(9999)
        print(out_llr)

    def make_edges_with_bpsk(self):
        for i in range(len(self.edg)):
            for j in range(len(self.edg[i])):
                self.edg[i][j] = (self.edg[i][j][0], -2 * int(self.edg[i][j][1]) + 1, self.edg[i][j][2])

import bcjr_decoder
from copy import deepcopy
import numpy as np
import mpmath as mp
from util import gfn_array_to_str


class BCJRDecoder:
    def __init__(self, edg: list[list[list[tuple[str, int, str]]]]):
        """
        Инициализация BCJR декодера для линейных блочных кодов на основе решетки.
        :param edg: Массив с рёбрами (список кортежей (синдром1, цифра ребра, синдром2))
        """
        self.edg = edg

        del self.edg[-1]

        self.edg_bpsk = deepcopy(self.edg)
        self.make_edges_with_bpsk()

        for i in range(len(self.edg)):
            for j in range(len(self.edg[i])):
                prev_vex = gfn_array_to_str(self.edg[i][j][0])
                next_vex = gfn_array_to_str(self.edg[i][j][2])

                self.edg[i][j] = (prev_vex, int(self.edg[i][j][1]), next_vex)
                self.edg_bpsk[i][j] = (prev_vex, self.edg_bpsk[i][j][1], next_vex)

    def decode_cpp(self, llr_in, sigma2):
        def bcjr_decode_python_callback(llr):
            return self.decode(llr, sigma2)

        return bcjr_decoder.decode(bcjr_decode_python_callback, self.edg, self.edg_bpsk, llr_in, sigma2)

    def decode(self, llr_in, sigma2):
        a_priori = 0.5

        '''
        FORWARD
        '''
        # Берем структуру edg, меняя значения в ребре на гамму
        gammas = deepcopy(self.edg_bpsk)

        # Сохраняем альфы в виде списка из словарей. Каждый i-ый словарь - это, по сути, ярус
        # Ключ словаря - синдром, он на каждом ярусе уникален. Значение - посчитанная альфа
        # Такой подход обусловлен тем, что пробегаясь по ребрам, мы могли быстро извлечь альфу следующей вершины
        alphas = [{} for _ in range(len(gammas) + 1)]
        alphas[0][gfn_array_to_str(gammas[0][0][0])] = 1

        constant_coef = a_priori * (1 / (2 * mp.pi * sigma2))  # ВОПРОС. В матлабе эта величина еще возводится в степень N
        for i in range(len(gammas)):
            for j in range(len(gammas[i])):
                '''
                GAMMA SECTION
                '''
                diff = (llr_in[i] - self.edg_bpsk[i][j][1]) ** 2 / (2 * sigma2)
                cur_gamma = constant_coef * mp.exp(-diff)
                prev_vex = gfn_array_to_str(gammas[i][j][0])
                next_vex = gfn_array_to_str(gammas[i][j][2])

                gammas[i][j] = (prev_vex, cur_gamma, next_vex)

                # print("layer", i, "/ gamma =", cur_gamma)

                '''
                ALPHA SECTION
                '''
                new_alpha = cur_gamma * alphas[i][prev_vex]
                # Если вершина с синдромом уже была рассмотрена, то считаем альфу
                # как сумму раннее посчитанного значения и текущего.
                # Это случай когда в вершину введет 2 ребра
                if next_vex in alphas[i + 1]:
                    alphas[i + 1][next_vex] += new_alpha
                # В противном случае просто записываем значение
                else:
                    alphas[i + 1][next_vex] = new_alpha

            # Normalization
            summa = sum(alphas[i + 1].values())
            if summa != 0:
                for key in alphas[i + 1].keys():
                    alphas[i + 1][key] /= summa

        '''
        BACKWARD
        '''
        betas = [{} for _ in range(len(gammas) + 1)]
        betas[-1][gammas[0][0][0]] = 1  # в матлабе последний бета равен последнему альфа, а не харкод 1. Однако, в таком подходе в результате нормализации послдений альфа всегда будет 1

        llr_out = [0] * len(llr_in)

        for i in range(len(gammas) - 1, -1, -1):
            up, down = 0, 0  # числитель и знаменатель для рассчета llr
            for j in range(len(gammas[i])):
                '''
                BETA SECTION
                '''
                cur_gamma = gammas[i][j][1]
                # обратный проход - next и prev меняются местами
                next_vex = gammas[i][j][0]
                prev_vex = gammas[i][j][2]

                new_beta = cur_gamma * betas[i + 1][prev_vex]
                if next_vex in betas[i]:
                    betas[i][next_vex] += new_beta
                else:
                    betas[i][next_vex] = new_beta

                '''
                SIGMA SECTION
                '''
                # Берем альфу, которая слева от текущего ребра
                cur_alpha = alphas[i][next_vex]
                # Берем бету, которая справа от текущего ребра
                cur_beta = betas[i + 1][prev_vex]

                cur_sigma = cur_gamma * cur_alpha * cur_beta

                '''
                LLR SECTION
                '''
                # не забываем про БПСК, по факту ищем ребра с 0
                if self.edg_bpsk[i][j][1] == 1:
                    up += cur_sigma
                else:
                    down += cur_sigma

            # print("up=", up, " down=", down, " log=", np.log(up / down), end=" // ", sep="")
            # llr_out[i] = np.log(up / down)
            # такое вообще возможно?
            if down == 0:
                llr_out[i] = 9999
            else:
                llr_out[i] = mp.ln(up / down)

            # Normalization beta
            summa = sum(betas[i].values())
            if summa != 0:
                for key in betas[i].keys():
                    betas[i][key] /= summa

        return llr_out

    def make_edges_with_bpsk(self):
        for i in range(len(self.edg_bpsk)):
            for j in range(len(self.edg_bpsk[i])):
                self.edg_bpsk[i][j] = (self.edg_bpsk[i][j][0], -2 * int(self.edg_bpsk[i][j][1]) + 1, self.edg_bpsk[i][j][2])

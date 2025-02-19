import numpy as np
from scipy.special import logsumexp

class BCJRDecoder:
    def __init__(self, vex, edg, H, sigma=1.0):
        """
        Инициализация BCJR декодера для линейных блочных кодов на основе решетки.
        :param vex: Массив с синдромами (список списков символов)
        :param edg: Массив с рёбрами (список кортежей (синдром1, цифра ребра, синдром2))
        :param H: Проверочная матрица (numpy.ndarray)
        :param sigma: Стандартное отклонение шума (по умолчанию 1.0)
        """
        self.vex = vex
        self.edg = edg
        self.H = H
        self.sigma = sigma
        self.snr = 1 / (2 * sigma**2)

    def compute_llr(self, y):
        """
        Вычисляет LLR (Log Likelihood Ratio) для входных данных.
        :param y: Принятый вектор (numpy.ndarray)
        :return: Вектор LLR (numpy.ndarray)
        """
        return 2 * self.snr * y

    def decode(self, y):
        """
        Выполняет декодирование BCJR для решётки.
        :param y: Принятый вектор (numpy.ndarray)
        :return: Декодированный кодовый вектор (numpy.ndarray)
        """
        llr = self.compute_llr(y)
        
        # Инициализация сообщений по решетке
        alpha = {v: -np.inf for v in self.vex}  # Прямой проход
        # alpha = {tuple(map(int, v)): -np.inf for v in self.vex}

        beta = {v: -np.inf for v in self.vex}   # Обратный проход
        
        alpha[self.vex[0]] = 0  # Начальное состояние
        
        # Прямой проход (форвард)
        for synd1, bit, synd2 in self.edg:
            alpha[synd2] = logsumexp([alpha[synd2], alpha[synd1] + bit * llr[bit]])
        
        beta[self.vex[-1]] = 0  # Конечное состояние
        
        # Обратный проход (бэкуорд)
        for synd1, bit, synd2 in reversed(self.edg):
            beta[synd1] = logsumexp([beta[synd1], beta[synd2] + bit * llr[bit]])
        
        # Вычисление апостериорных LLR
        L = np.array([alpha[s1] + beta[s2] for s1, bit, s2 in self.edg])
        
        return (L < 0).astype(int)  # Преобразование в бинарный код

# Пример использования:
if __name__ == "__main__":
    # Пример решетки
    vex = [[0], [1]]  # Синдромы
    edg = [([0], 0, [1]), ([1], 1, [0])]  # Рёбра (синдром1, цифра ребра, синдром2)
    
    # Проверочная матрица для (7,4) кода Хэмминга
    H = np.array([[1, 1, 0, 1, 1, 0, 0],
                  [1, 0, 1, 1, 0, 1, 0],
                  [0, 1, 1, 1, 0, 0, 1]])
    
    decoder = BCJRDecoder(vex, edg, H, sigma=0.5)
    received = np.array([1.2, -0.8, 0.5, -1.1, 0.9, -0.7, 0.3])  # Пример принятых значений
    decoded = decoder.decode(received)
    print("Декодированное слово:", decoded)

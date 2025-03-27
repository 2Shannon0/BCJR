from BCJR import BCJRDecoder
import bcjr_decoder
from trellis4decoder import Trellis
from trellis_repo import get_trellis
from bpsk import bpsk_modulation, bpsk_demodulation
from awgn import awgn_llr
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import time
# def awgn(c, sigma=0.5):
#     # rewrite
#     # 1 -> -1
#     # 0 -> +1
#     c1 = c * -2 + 1

#     #AWGN
#     loc = 0,
#     sigma = np.sqrt(sigma)
#     size = len(c1)
#     e = np.array(np.random.normal(loc, sigma, size)) # 2y / sigma^2

#     c2 = c1+e

#     # rewrite vector with errors in -1 +1
#     c3 = []

#     for i in c2:
#         if i < 0: i = -1
#         if i > 0: i = 1
#         c3.append(i)
#     c3 = np.array(c3)

#     #determining the number of errors
#     # def findNumberOfE(a,b):
#     # k = 0
#     # for i in c1:
#     #   if c1[i] == c3[i]: k += 1
#     errorIndexcies = []
#     for i in range(c3.shape[0]):
#         if c3[i] != c1[i]: errorIndexcies.append(i)

#     cLLR =[]
#     for i in c2:
#         i = 2*i/(sigma**2)
#         cLLR.append(i)
#     cLLR = np.array(cLLR)

#     return cLLR, errorIndexcies


# Раскоментить, если нет закэшированной решетки
# trellis = Trellis("/home/i17m5/BCJR/matricies/BCH_MATRIX_N_15_K_7_DEFAULT.csv")
# trellis.build_trellis()

trellis = get_trellis("BCH_MATRIX_N_63_K_51_DEFAULT")

N = len(trellis.vex) - 1

# Задаем кодовое слово
# codeword_initial = [1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0] # BCH(15, 7)
# codeword_initial = [1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1] # BCH(15, 7)
codeword_initial = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# codeword_initial = [0] * N
# codeword_initial = [0, 0, 0, 0, 0, 0, 0] 
print(f'Кодовое слово: {codeword_initial}\n\n')

codeword_modulated = bpsk_modulation(codeword_initial)
print(f'Кодовое слово BPSK: {codeword_modulated}\n\n')

# Инициализируем декодер
start_time = time.time()
python_decoder = BCJRDecoder(trellis.edg)
end_time = time.time()
print(f"Инициализаци python-декодера: {end_time-start_time:.6f} секунд")

# Для заданного отношения сигнал-шум считаем llr
llr_in, sigma2 = awgn_llr(codeword_modulated, 0)
print(f'Значение sigma2: {sigma2}\n\n')
print(f'Входные LLR: {llr_in}\n\n')

# llr после декодирования
start_time = time.time()
llr_out = python_decoder.decode(llr_in, sigma2)
end_time = time.time()
print(f"Python декодирование: {end_time-start_time:.6f} секунд")

start_time = time.time()
result_cpp = bcjr_decoder.decode(python_decoder.edg, python_decoder.edg_bpsk, llr_in, sigma2)
end_time = time.time()
print(f"CPP декодирование: {end_time-start_time:.6f} секунд")

print(f'Выходные LLR: {llr_out}\n\n')
print(f'result_cpp: {result_cpp}\n\n')


# Декодированное кодовое слово в бинарном виде
codeword_result = bpsk_demodulation(llr_out)

# считаем кол-во ошибок
errors = 0
for j in range(N):
    if codeword_result[j] != codeword_initial[j]:
        errors += 1

print(f'Количество ошибок: {errors}\n\n')
print(f'Декодированный вектор:\n{codeword_result}')
print(f'{codeword_initial}\nкодовый вектор^^^')

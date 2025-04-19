import ctypes
import numpy as np
import os
from trellis4decoder import Trellis
from trellis_repo import get_trellis
from bpsk import bpsk_modulation, bpsk_demodulation
from awgn import awgn_llr, awgn_llr_complex
from BCJR import BCJRDecoder
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import time

start = time.time()

ESNO_START = 0
ESNO_END = 10
ESNO_STEP = 0.2
WRONG_DECODING_NUMBER = 120

# Раскоментить, если нет закэшированной решетки
# trellis = Trellis("/Users/aleksejbandukov/Documents/python/BCJR_Project/matricies/BCH_MATRIX_N_15_K_11_DEFAULT.csv")
# trellis.build_trellis()
trellis_name = 'BCH_MATRIX_N_31_K_16_DEFAULT'
trellis = get_trellis(f'./trellis_binaries/{trellis_name}')

N = len(trellis.vex) - 1

TITLE = f'Decoding BCJR, trellis: {trellis_name} WRONG_DECODING_NUMBER = {WRONG_DECODING_NUMBER}, ESNO_END = {ESNO_END}'
print('\n',TITLE,'\n')

# Задаем кодовое слово
# codeword_initial = [0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0] # BCH(15, 5)
# codeword_initial = [1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0] # BCH(15, 7)
# codeword_initial = [1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0] # BCH(15, 11)
codeword_initial = [1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0]   # BCH(31, 16)
# codeword_initial = [0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0] # BCH(31, 26)
# codeword_initial = [0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0] # BCH(31, 21)

# codeword_initial = [0] * N

codeword_modulated = bpsk_modulation(codeword_initial)

# Задаем список EsNo
esno_array = []
value = ESNO_START
while round(value, 2) <= ESNO_END:
    esno_array.append(round(value, 2))
    value += ESNO_STEP

# Создаем fer & ber
fer = [0] * len(esno_array)
ber = [0] * len(esno_array)

# Инициализируем декодер
decoder = BCJRDecoder(trellis.edg)

# # Запускаем моделирование
for (i, esno) in enumerate(esno_array):
    tests_passed, wrong_decoding, errors_at_all = 0, 0, 0

    print(f"\n-------------------- EsNo = {esno} --------------------")

    while wrong_decoding < WRONG_DECODING_NUMBER:
        tests_passed += 1

        # Для заданного отношения сигнал-шум считаем llr
        llr_in, sigma2 = awgn_llr_complex(codeword_modulated, esno)

        # llr после декодирования
        llr_out = decoder.decode_cpp(llr_in, sigma2)

        # Декодированное кодовое слово в бинарном виде
        codeword_result = bpsk_demodulation(llr_out)

        # считаем кол-во ошибок
        errors = 0
        for j in range(N):
            if codeword_result[j] != codeword_initial[j]:
                errors += 1

        # если ошибки есть, то считаем fer & ber
        if errors != 0:
            wrong_decoding += 1
            errors_at_all += errors

            fer[i] = wrong_decoding / tests_passed
            ber[i] = errors_at_all / N / tests_passed

            print(f"fer = {fer[i]}, ber = {ber[i]}, tests_passed = {tests_passed}")

    print("\nRESULTS")
    print(esno_array)
    print(fer)
    print(ber)

end = time.time()

print(f"Время выполнения: {end - start:.4f} секунд")

fer_smooth = gaussian_filter1d(fer, sigma=1).tolist() # Параметр sigma овечает за то, насколько сильно сглаживать график. При 2 выглядит оптимально

plt.plot(esno_array, fer, label="Original", alpha=0.5, linewidth=1)
plt.plot(esno_array, fer_smooth, label="Smoothed", linewidth=2)
plt.yscale("log")  # Логарифмическая шкала по Y
plt.xlabel("EsNo")
plt.ylabel("FER")
plt.legend()
plt.grid(True, which="both", linestyle="--")
# plt.show()
# plt.savefig(f'../modeling_results/BCJR_{trellis_name}_from_{ESNO_START}_to_{ESNO_END}.png', dpi=300, bbox_inches='tight')


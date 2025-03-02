from BCJR import BCJRDecoder
from simulation.trellis4decoder import Trellis
from simulation.bpsk import bpsk_modulation, bpsk_demodulation
from simulation.awgn import awgn_llr
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d


ESNO_START = -12
ESNO_END = -2
ESNO_STEP = 0.4
WRONG_DECODING_NUMBER = 5
SUPERCODE_ITERATIONS = 5


trellis1 = Trellis("../matricies/BCH_MATRIX_N_15_K_7_HALF_1.csv")
trellis1.build_trellis()
trellis2 = Trellis("../matricies/BCH_MATRIX_N_15_K_7_HALF_2.csv")
trellis2.build_trellis()

N = len(trellis1.vex) - 1

# Задаем нулевое кодовое слово
codeword_initial = [0] * N
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

decoder1 = BCJRDecoder(trellis1.vex, trellis1.edg)
decoder2 = BCJRDecoder(trellis2.vex, trellis2.edg)

for (i, esno) in enumerate(esno_array):
    tests_passed, wrong_decoding, errors_at_all = 0, 0, 0

    print(f"\n-------------------- EsNo = {esno} --------------------")

    while wrong_decoding < WRONG_DECODING_NUMBER:
        tests_passed += 1

        # Для заданного отношения сигнал-шум считаем llr
        llr, sigma2 = awgn_llr(codeword_modulated, esno)

        llr_prev_super_1 = [0] * N
        llr_prev_super_2 = [0] * N

        for _ in range(SUPERCODE_ITERATIONS):

            '''
            1-st supercode
            '''
            llr_to_decode = [0] * N

            for j in range(N):
                llr_to_decode[j] = llr[j] - llr_prev_super_1[j]

            llr_out_1 = decoder1.decode(llr_to_decode, sigma2)

            for j in range(N):
                llr_out_1[j] -= llr[j] - llr_to_decode[j]
                llr_prev_super_1[j] = llr_out_1[j]

            '''
            2-nd supercode
            '''
            llr_to_decode = [0] * N

            for j in range(N):
                llr_to_decode[j] = llr[j] - llr_prev_super_2[j]

            llr_out_2 = decoder2.decode(llr_to_decode, sigma2)

            for j in range(N):
                llr_out_2[j] -= llr[j] - llr_to_decode[j]
                llr_prev_super_2[j] = llr_out_2[j]


            for j in range(N):
                llr[j] += llr_out_1[j] + llr_out_2[j]


        # Декодированное кодовое слово в бинарном виде
        codeword_result = bpsk_demodulation(llr)

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
            ber[i] = wrong_decoding / N / tests_passed

            print(f"fer = {fer[i]}, ber = {ber[i]}, tests_passed = {tests_passed}")

print("\nRESULTS")
print(esno_array)
print(fer)
print(ber)

fer_smooth = gaussian_filter1d(fer, sigma=1).tolist()

plt.plot(esno_array, fer, label="Original", alpha=0.5, linewidth=1)
plt.plot(esno_array, fer_smooth, label="Smoothed", linewidth=2)
plt.yscale("log")  # Логарифмическая шкала по Y
plt.xlabel("EsNo")
plt.ylabel("FER")
plt.legend()
plt.grid(True, which="both", linestyle="--")
plt.show()
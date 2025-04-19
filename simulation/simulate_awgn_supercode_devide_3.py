from BCJR import BCJRDecoder
from trellis_repo import get_trellis
from trellis4decoder import Trellis
from bpsk import bpsk_modulation, bpsk_demodulation
from awgn import awgn_llr
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from copy import deepcopy
from random import randrange
import numpy as np
import concurrent.futures
from trellis_repo import get_trellis

if __name__ == "__main__":

    ESNO_START = -10
    ESNO_END = 10
    ESNO_STEP = 0.2
    WRONG_DECODING_NUMBER = 120
    SUPERCODE_ITERATIONS = 4

    # trellis1 = Trellis("../matricies/BCH_MATRIX_N_15_K_7_PART_1_3.csv")
    # trellis1.build_trellis()
    # trellis2 = Trellis("../matricies/BCH_MATRIX_N_15_K_7_PART_2_3.csv")
    # trellis2.build_trellis()
    # trellis3 = Trellis("../matricies/BCH_MATRIX_N_15_K_7_PART_3_3.csv")
    # trellis3.build_trellis()
    trellis1 = get_trellis('/home/k111/BCJR_complex/simulation/trellis_binaries/BCH_MATRIX_N_31_K_16_part_1_3')
    trellis2 = get_trellis('/home/k111/BCJR_complex/simulation/trellis_binaries/BCH_MATRIX_N_31_K_16_part_2_3')
    trellis3 = get_trellis('/home/k111/BCJR_complex/simulation/trellis_binaries/BCH_MATRIX_N_31_K_16_part_3_3')


    N = len(trellis1.vex) - 1

    # Задаем нулевое кодовое слово
    # codeword_initial = [0] * N
    # codeword_initial = [0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0]  # BCH(15, 5)
    # codeword_initial = [1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0] # BCH(15, 7)
    codeword_initial = [1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0]   # BCH(31, 16)
    # codeword_initial = [0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0] # BCH(31, 26)

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

    decoder1 = BCJRDecoder(trellis1.edg)
    decoder2 = BCJRDecoder(trellis2.edg)
    decoder3 = BCJRDecoder(trellis3.edg)
    TITLE = f'Decoding SUPERCODE, WRONG_DECODING_NUMBER = {WRONG_DECODING_NUMBER}, ESNO_END = {ESNO_END}, iter: {SUPERCODE_ITERATIONS}, matrix: 3, BCH(31,16)'
    print('\n',TITLE,'\n')
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for (i, esno) in enumerate(esno_array):
            tests_passed, wrong_decoding, errors_at_all = 0, 0, 0
            print()
            print(esno_array)
            print(fer)
            print(ber)
            print(f"\n-------------------- EsNo = {esno} --------------------")

            # Создаём пул процессов ОДИН раз перед началом тестов
            while wrong_decoding < WRONG_DECODING_NUMBER:
                tests_passed += 1

                # Для заданного отношения сигнал-шум считаем llr
                llr, sigma2 = awgn_llr(codeword_modulated, esno)

                llr_initial = deepcopy(llr)
                llr_result = [0] * N

                llr_prev_super_1 = [0] * N
                llr_prev_super_2 = [0] * N
                llr_prev_super_3 = [0] * N

                for k in range(SUPERCODE_ITERATIONS):
                    for j in range(N):
                        llr[j] = llr_initial[j] + llr_result[j]

                    llr_to_decode_1 = [llr[j] - llr_prev_super_1[j] for j in range(N)]
                    llr_to_decode_2 = [llr[j] - llr_prev_super_2[j] for j in range(N)]
                    llr_to_decode_3 = [llr[j] - llr_prev_super_3[j] for j in range(N)]

                    # Запускаем два декодера в параллельных процессах
                    future1 = executor.submit(decoder1.decode, llr_to_decode_1, sigma2)
                    future2 = executor.submit(decoder2.decode, llr_to_decode_2, sigma2)
                    future3 = executor.submit(decoder3.decode, llr_to_decode_3, sigma2)

                    # Ожидание результатов
                    llr_out_1 = future1.result()
                    llr_out_2 = future2.result()
                    llr_out_3 = future3.result()

                    # Обработка результатов
                    for j in range(N):
                        llr_out_1[j] -= llr_to_decode_1[j]
                        llr_prev_super_1[j] = llr_out_1[j]

                        llr_out_2[j] -= llr_to_decode_2[j]
                        llr_prev_super_2[j] = llr_out_2[j]

                        llr_out_3[j] -= llr_to_decode_3[j]
                        llr_prev_super_3[j] = llr_out_3[j]

                        llr_result[j] = llr_out_1[j] + llr_out_2[j] + llr_out_3[j]

                # Декодированное кодовое слово в бинарном виде
                codeword_result = bpsk_demodulation(llr_result)


                # считаем кол-во ошибок
                errors = 0
                for j in range(N):
                    if codeword_result[j] != codeword_initial[j]:
                        errors += 1
                # print("error left", errors)

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

    fer_smooth = gaussian_filter1d(fer, sigma=1).tolist()

    plt.plot(esno_array, fer, label="Original", alpha=0.5, linewidth=1)
    plt.plot(esno_array, fer_smooth, label="Smoothed", linewidth=2)
    plt.yscale("log")  # Логарифмическая шкала по Y
    plt.xlabel("EsNo")
    plt.ylabel("FER")
    plt.legend()
    plt.grid(True, which="both", linestyle="--")
    # plt.show()
import numpy as np


def awgn_llr(codeword, snr_db):
    signal = np.array(codeword)
    P_signal = np.mean(signal**2)  # Мощность сигнала (BPSK = 1)
    P_noise = P_signal * 10 ** (-snr_db / 10)  # Мощность шума
    sigma2 = P_noise  # Дисперсия шума
    noise = np.random.normal(0, np.sqrt(sigma2 / 2), size=signal.shape)  # AWGN
    received_signal = signal + noise  # Сигнал после канала

    # Вычисление LLR (по формуле LLR = 2y / sigma^2)
    llr_values = 2 * received_signal / sigma2

    return llr_values.astype(list), sigma2


def awgn_llr_complex(codeword, snr_db):
    signal = np.array(codeword, dtype=complex)
    noiseVar = 10 ** (-snr_db / 10) # дисперсия sigma^2
    # sigma2 = noiseVar
    noise = (np.random.randn(*signal.shape) + 1j * np.random.randn(*signal.shape)) / np.sqrt(2)

    received_signal = signal + np.sqrt(noiseVar) * noise

    llr_values = 2 * received_signal.real / noiseVar

    return llr_values.astype(list), noiseVar

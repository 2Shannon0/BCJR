def bpsk_modulation(codeword: list) -> list:
    return [-2 * c + 1 for c in codeword]


def bpsk_demodulation(llr: list) -> list:
    return [int(lr < 0) for lr in llr]
import numpy as np

def generate_gallager_ldpc(K, N, w_c, w_r):
    assert N > K, "N должно быть больше K"
    M = N - K
    assert M * w_r == N * w_c, "Несоответствие плотностей: M*w_r должно быть равно N*w_c"

    J = w_c
    L = w_r

    rows_per_block = M // J
    H_blocks = []

    H0 = np.zeros((rows_per_block, N), dtype=int)
    cols_per_row = N // L

    for i in range(rows_per_block):
        for j in range(L):
            col = (i * L + j) % N
            H0[i, col] = 1
    H_blocks.append(H0)

    for _ in range(1, J):
        P = np.random.permutation(N)
        H_perm = H0[:, P]
        H_blocks.append(H_perm)

    H = np.vstack(H_blocks)

    return H


def gcd(first, second):
    assert first > second

    if second == 0:
        return first

    return gcd(second, first % second)


def find_possible_weights(N, K):
    M = N - K

    g_c_d = gcd(N, M)

    up = M // g_c_d
    down = N // g_c_d

    print("w_c", "w_r", "dencity", end="\n", sep="\t")
    for i in range(M):
        w_c = up * (i + 1)
        w_r = down * (i + 1)
        density = w_c * w_r / (N * M)
        if density > 1:
            break
        print(w_c, w_r, density, end="\n", sep="\t")


def print_matrix(m):
    for i in range(len(m)):
        for j in range(len(m[0])):
            print(m[i][j], sep="", end="\t")
        print()


K = 196
N = 420

find_possible_weights(N, K)

w_c = 8
w_r = 15

H = generate_gallager_ldpc(K, N, w_c, w_r)
print_matrix(H)
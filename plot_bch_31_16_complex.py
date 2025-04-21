import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

#-----------------------BCH(31,26) original----------------------------------------------------------------
esno_bch_31_16 = []
fer_bch_31_16 = []
ber_bch_31_16 = []
plt.plot(esno_bch_31_16, fer_bch_31_16, label="BCJR(31,16)", linewidth=2)
# plt.plot(esno_bch_31_16, ber_bch_31_16, label="BCJR(31,26)", linewidth=2)


#-----------------------BCH(31,26) supercode----------------------------------------------------------------
esno_bch_31_16_supercode_2_matrix_1_iter = []
fer_bch_31_16_supercode_2_matrix_1_iter = []
ber_bch_31_16_supercode_2_matrix_1_iter = []
plt.plot(esno_bch_31_16_supercode_2_matrix_1_iter, fer_bch_31_16_supercode_2_matrix_1_iter, label="BCJR(31,16) supercode(2,1)", linewidth=2)
# plt.plot(esno_bch_31_16_supercode_2_matrix_1_iter, ber_bch_31_16_supercode_2_matrix_1_iter, label="BCJR(31,16) supercode(2,1)", linewidth=2)


esno_bch_31_16_supercode_2_matrix_2_iter = []
fer_bch_31_16_supercode_2_matrix_2_iter = []
ber_bch_31_16_supercode_2_matrix_2_iter = []
plt.plot(esno_bch_31_16_supercode_2_matrix_2_iter, fer_bch_31_16_supercode_2_matrix_2_iter, label="BCJR(31,16) supercode(2,2)", linewidth=2)
# plt.plot(esno_bch_31_16_supercode_2_matrix_2_iter, ber_bch_31_16_supercode_2_matrix_2_iter, label="BCJR(31,16) supercode(2,2)", linewidth=2)


esno_bch_31_16_supercode_2_matrix_3_iter = []
fer_bch_31_16_supercode_2_matrix_3_iter = []
ber_bch_31_16_supercode_2_matrix_3_iter = []
plt.plot(esno_bch_31_16_supercode_2_matrix_3_iter, fer_bch_31_16_supercode_2_matrix_3_iter, label="BCJR(31,16) supercode(2,3)", linewidth=2)
# plt.plot(esno_bch_31_16_supercode_2_matrix_3_iter, ber_bch_31_16_supercode_2_matrix_3_iter, label="BCJR(31,16) supercode(2,3)", linewidth=2)

esno_bch_31_16_supercode_2_matrix_4_iter = []
fer_bch_31_16_supercode_2_matrix_4_iter = []
ber_bch_31_16_supercode_2_matrix_4_iter = []
plt.plot(esno_bch_31_16_supercode_2_matrix_4_iter, fer_bch_31_16_supercode_2_matrix_4_iter, label="BCJR(31,16) supercode(2,4)", linewidth=2)
# plt.plot(esno_bch_31_16_supercode_2_matrix_4_iter, ber_bch_31_16_supercode_2_matrix_4_iter, label="BCJR(31,16) supercode(2,4)", linewidth=2)


esno_bch_31_16_supercode_2_matrix_5_iter = []
fer_bch_31_16_supercode_2_matrix_5_iter = []
ber_bch_31_16_supercode_2_matrix_5_iter = []
plt.plot(esno_bch_31_16_supercode_2_matrix_5_iter, fer_bch_31_16_supercode_2_matrix_5_iter, label="BCJR(31,16) supercode(2,5)", linewidth=2)
# plt.plot(esno_bch_31_16_supercode_2_matrix_5_iter, ber_bch_31_16_supercode_2_matrix_5_iter, label="BCJR(31,16) supercode(2,5)", linewidth=2)


plt.yscale("log")  # Логарифмическая шкала по Y
plt.xlabel("EsNo [dB]")
plt.ylabel("FER")
plt.title("BCH(31, 26)")
plt.legend()
plt.grid(True, which="both", linestyle="--")
plt.savefig('modeling_results/BCH(31,26)_complex_awgn_FER.png', dpi=300, bbox_inches='tight')
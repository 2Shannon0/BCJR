import bcjr_decoder
from trellis4decoder import Trellis
from BCJR import BCJRDecoder




def main():
    print(dir(bcjr_decoder))
    t_h = Trellis('/Users/aleksejbandukov/Documents/python/BCJR_Project/matricies/file_hamming.csv')
    t_h.build_trellis()

    bcjr_python = BCJRDecoder(t_h.edg)

    llr_in = [1, 1.2, -0.5, 1, 1, 1, 0.5]
    sigma2 = 0.8

    # Вызываем функцию
    result_c = bcjr_decoder.decode(bcjr_python.edg, bcjr_python.edg_bpsk, llr_in, sigma2)

    result_pyt = bcjr_python.decode(llr_in, sigma2)

    print("Decoded Output with c:", result_c)
    print("Decoded Output with python:", result_pyt)


if __name__ == "__main__":
    main()

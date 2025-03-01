import pickle
from trellis4decoder import Trellis


FILE_PATH = "../matricies/BCH_MATRIX_N_31_K_16_DEFAULT.csv"
TRELLIS_NAME = "trellis_bch_31_16"


def save_trellis(file_path, trellis_name):
    t_h = Trellis(matrix_file=file_path)
    t_h.build_trellis()

    with(open(trellis_name, "wb")) as f:
        pickle.dump(t_h, f)


def get_trellis(trellis_name):
    with(open(trellis_name, "rb")) as f:
        trellis = pickle.load(f)
        return trellis


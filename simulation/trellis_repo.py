import pickle
from trellis4decoder import Trellis


FILE_PATH = "/home/k111/BCJR/matricies/BCH_MATRIX_N_31_K_16_half_2.csv"
TRELLIS_NAME = "BCH_MATRIX_N_31_K_16_half_2"


def save_trellis(file_path, trellis_name):
    t_h = Trellis(matrix_file=file_path)
    t_h.build_trellis()
    # t_h.plot_sections()

    with(open(f'trellis_binaries/{trellis_name}', "wb")) as f:
        pickle.dump(t_h, f)
    print(f'\nРешетка "{trellis_name}" сохранена!')


def get_trellis(trellis_name):
    with(open(trellis_name, "rb")) as f:
        trellis = pickle.load(f)
        return trellis

# save_trellis(FILE_PATH, TRELLIS_NAME)

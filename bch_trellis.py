from simulation.trellis4decoder import Trellis

t_h = Trellis('matricies/BCH_MATRIX_N_15_K_7_DEFAULT.csv')

t_h.build_trellis()
t_h.plot_sections()

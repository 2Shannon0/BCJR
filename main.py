from BCJR_visible import BCJRDecoder
from simulation.trellis4decoder import Trellis
import numpy as np


t_h = Trellis('matricies/file_hamming.csv')

t_h.build_trellis()
# t_h.plot_sections()
print('len(t_h.vex): ',len(t_h.vex),'\n')
print('t_h.vex[0]): ',t_h.vex[0],'\n')
print('t_h.vex[1]): ',t_h.vex[1],'\n')
print('t_h.vex[2]): ',t_h.vex[2],'\n')


# print(t_h.edg)

# t_h2 = Trellis(vex=t_h.vex, edg=t_h.edg)

# t_h2.plot_sections()

H = np.array([[1, 1, 0, 1, 1, 0, 0],
              [1, 0, 1, 1, 0, 1, 0],
              [0, 1, 1, 1, 0, 0, 1]])

decoder = BCJRDecoder(t_h.vex, t_h.edg, H, sigma=0.5)
# received = np.array([1.2, -0.8, 0.5, -1.1, 0.9, -0.7, 0.3])
decoder.decode([1, 1.2, -0.5, 1, 1, 1, 0.5], 0.8)
# print("Декодированное слово:", decoded)
# received = np.array([1.2, -0.8, 0.5, -1.1, 0.9, -0.7, 0.3])  # Пример принятых значений
# decoded = decoder.decode(received)
# print("Декодированное слово:", decoded)
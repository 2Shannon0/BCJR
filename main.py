from BCJR import BCJRDecoder
from trellis4decoder import Trellis
import numpy as np


t_h = Trellis('matricies/file_hamming.csv')

t_h.build_trellis()
# t_h.plot_sections()
# print(t_h.vex)
# print(t_h.edg)

H = np.array([[1, 1, 0, 1, 1, 0, 0],
              [1, 0, 1, 1, 0, 1, 0],
              [0, 1, 1, 1, 0, 0, 1]])
    
decoder = BCJRDecoder(t_h.vex, t_h.edg, H, sigma=0.5)
received = np.array([1.2, -0.8, 0.5, -1.1, 0.9, -0.7, 0.3])  # Пример принятых значений
decoded = decoder.decode(received)
print("Декодированное слово:", decoded)
from simulation.trellis4decoder import Trellis

t_h = Trellis('matricies/LDPC_16_10.csv')

t_h.build_trellis()

for i in range(len(t_h.edg)):
    for j in range(len(t_h.edg[i])):
        print(t_h.edg[i][j])
    print()

t_h.plot_sections()

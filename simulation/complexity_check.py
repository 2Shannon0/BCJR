from trellis_repo import get_trellis
from trellis4decoder import Trellis

bch_orig_trellis = Trellis("../matricies/BCH_MATRIX_N_15_K_7_DEFAULT.csv")
bch_orig_trellis.build_trellis()

edge_num_1 = 0

del bch_orig_trellis.edg[-1]

for layer in bch_orig_trellis.edg:
    for edge in layer:
        edge_num_1 += 1

bch_super_trellis = Trellis("../matricies/BCH_MATRIX_N_15_K_7_PART_1_4.csv")
bch_super_trellis.build_trellis()

edge_num_2 = 0

del bch_super_trellis.edg[-1]

for layer in bch_super_trellis.edg:
    for edge in layer:
        edge_num_2 += 1

print(edge_num_1, edge_num_2)

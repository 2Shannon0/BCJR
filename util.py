import csv
import numpy as np

# Находит вершины из списка vex, которые соединены с другими вершинами в edg
def vex_connected( vex, edg ):
    edgs = set()
    for e in edg:
        v0, a, v1 = e
        v0_str = gfn_array_to_str(v0)
        if v0_str in vex:
            edgs.add(v0_str)

    return edgs

def gfn_array_to_str(gfns) -> str:
    result = ""
    for gfn in gfns:
        result += str(gfn)
    return result

# Преобразует массив чисел в строку, добавляя индекс i, если он указан.
def v2str( in_arr, i=None ):
    name = ""
    if i is not None: name += str(i) + '_'
    for s in in_arr:    
        name += str(int(s))
    return name

# Проверяет, содержится ли target в vec_list.
def in_list( vec_list, target ):
    for vec in vec_list:
        if np.array_equal(vec,target): return True
    return False

# Загружает CSV-файл в numpy-массив.
def read_csv( filename ):
    with open( filename, newline='') as f:
        reader = csv.reader(f)
        symbols = []
        for row in reader:
            rows = []
            for symbol in row:
                rows.append(int(symbol))
            symbols.append(rows)
    return np.array(symbols)

# Загружает CSV-файл, транспонирует данные и адаптирует их под код GFn.
def read_mat( filename ):
    with open( filename, newline='') as f:
        reader = csv.reader(f)
        symbols = []
        for row in reader:
            rows = []
            for symbol in row:
                symbol_str = []
                for gf in symbol:
                    symbol_str.append(int(gf))
                rows.append(symbol_str)
            symbols.append(rows)
        rows_nparr = np.swapaxes( symbols, 0, 1 )
        size = int(int(rows_nparr.size)/len(symbols))
    return rows_nparr

# преобразует массив битов в целое число
def arr2int( v ):
    if isinstance(v, str):
        return int(v, 2)
    total = 0
    for e, digit in enumerate(v):
        total *= 2**digit.nbit
        total += int(digit)

    return total

# Возвращает примитивный полином для GF(2^n)
def prmt_ply( nbit ):
    if nbit==1: return np.array([0,1])
    if nbit==2: return np.array([1,1,1])
    if nbit==3: return np.array([1,1,0,1])
    if nbit==4: return np.array([1,1,0,0,1])
    if nbit==5: return np.array([1,0,1,0,0,1])
    if nbit==6: return np.array([1,1,0,0,0,0,1])
    if nbit==7: return np.array([1,0,0,1,0,0,0,1])
    err_msg = "No primitive polynomial for nbit = ", str(nbit)
    raise ValueError(err_msg)

# Выполняет деление с остатком в поле GF(2), важно для операций в GFn.
def gf2_remainder( a, b ):
    while 1:

        # Remainder is 0, return with padding or slicing
        if np.argwhere(a==1).size == 0:
            if a.size < b.size:
                return np.append( a, np.zeros(b.size-a.size-1, dtype=int) )
            else:
                return a[:b.size-1]

        msb = np.max(np.argwhere(a==1),axis=0)[0]

        # Degree of remainder is less than divisor
        if msb < len(b)-1:
            if a.size < b.size:
                return np.append( a, np.zeros(b.size-a.size-1, dtype=int) )
            else:
                return a[:b.size-1]

        # Do padding to align the MSB of remainder to the dividend
        # remainder = np.append( b, np.zeros( msb - b.size + 1 ) )
        # remainder = np.append( np.zeros(a.size - msb - 1), remainder )
        remainder = np.append( np.zeros( msb - b.size + 1 ), b )
        remainder = np.append( remainder, np.zeros(a.size - msb - 1)).astype(int)
        # print("remainder = ", remainder.flatten())        
        result = np.empty_like(a)
        for i, x in enumerate(result):
            result[i] = np.remainder( a[i]+remainder[i], 2)
        a = result

# Редуцирует a по модулю примитивного полинома, чтобы остаться в пределах GF(2^n). Используется в арифметических операциях GFn
def fit_gfn( a, nbit ):
    b = prmt_ply(nbit)
    a = np.remainder(a,2).astype(int)
    if a.size < b.size:
        return np.append( a, np.zeros(b.size-a.size-1, dtype=int) )
    else:
        return gf2_remainder(a,b)
    
def print_matrix(m):
    for i in range(m.shape[0]):
      for j in range(m.shape[1]):
          print(int(m[i][j]), end=" ")
      print()

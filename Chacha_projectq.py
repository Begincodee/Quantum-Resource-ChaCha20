from projectq import MainEngine
from projectq.ops import *
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator, Simulator
from projectq.meta import Loop, Compute, Uncompute, Control, LogicalQubitIDTag
from math import pi

def MAJ(eng, a, b, c):
    CNOT | (a, b)
    CNOT | (a, c)
    Toffoli | (c, b, a)
    return b
    
def UMA(eng, a, b, c):
    X | b
    CNOT | (c, b)
    Toffoli | (c, b, a)
    X | b
    CNOT | (a, c)
    CNOT | (a, b)
    return b

def Mod_add_QFT(eng, b, a, n):
    QFT | a
    for i in range(n):
        for j in range(n - i):
            with Control(eng, [b[j]]):
                if i >= j:
                    angle = 2 * pi / (2 ** (i - j + 1))
                    R(angle) | a[i]
    get_inverse(QFT) | a
    return a

def Mod_add(eng, a, b, n):
    c = eng.allocate_qubit()
    MAJ(eng, a[0], b[0], c)
    for i in range(1, n):
        MAJ(eng, a[i], b[i], a[i-1])
    for i in range(1, n):
        UMA(eng, a[n - i], b[n - i], a[n - i - 1])
    UMA(eng, a[0], b[0], c)
    return b


def cyclic_left_shift_16(eng, b):
    for i in range(16):
        Swap | (b[i], b[(i + 16)])
    return b

def cyclic_left_shift_12(eng, b):
    index = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3,
             4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
    new_b = []
    for i in range(32):
        new_b.append(b[index[i]])
    return new_b

def cyclic_left_shift_8(eng, b):
    index = [24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8,
             9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    new_b = []
    for i in range(32):
        new_b.append(b[index[i]])
    return new_b

def cyclic_left_shift_7(eng, b):
    index = [25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
             10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
    new_b = []
    for i in range(32):
        new_b.append(b[index[i]])
    return new_b

def quarter_round(eng, a, b, c, d):
    
    Mod_add(eng, b, a, 32)
    for i in range(32):
        CNOT | (a[i], d[i])
    d = cyclic_left_shift_16(eng, d)
    
    Mod_add(eng, d, c, 32)
    for i in range(32):
        CNOT | (c[i], b[i])
    b = cyclic_left_shift_12(eng, b)
    
    Mod_add(eng, b, a, 32)
    for i in range(32):
        CNOT | (a[i], d[i])
    d = cyclic_left_shift_8(eng, d)
    
    Mod_add(eng, d, c, 32)
    for i in range(32):
        CNOT | (c[i], b[i])
    b = cyclic_left_shift_7(eng, b)
    
    return [a, b, c, d]

def quarter_round_QFT(eng, a, b, c, d):
    
    Mod_add_QFT(eng, b, a, 32)
    for i in range(32):
        CNOT | (a[i], d[i])
    d = cyclic_left_shift_16(eng, d)
    
    Mod_add_QFT(eng, d, c, 32)
    for i in range(32):
        CNOT | (c[i], b[i])
    b = cyclic_left_shift_12(eng, b)
    
    Mod_add_QFT(eng, b, a, 32)
    for i in range(32):
        CNOT | (a[i], d[i])
    d = cyclic_left_shift_8(eng, d)
    
    Mod_add_QFT(eng, d, c, 32)
    for i in range(32):
        CNOT | (c[i], b[i])
    b = cyclic_left_shift_7(eng, b)
    
    return [a, b, c, d]

def ChaCha20(eng, key, const, count, nonce, out, n):
    
    for i in range(n * 16):
        if i < n * 4:
            CNOT | (const[i], out[i])
        elif  n * 4 <= i < n * 12:
            CNOT | (key[i - n * 4], out[i])
        elif n * 12 <= i < n * 13:
            CNOT | (count[i - n * 12], out[i])
        elif n * 13 <= i < n * 16:
            CNOT | (nonce[i - n * 13], out[i])
    
    for i in range(10):
        # Column rounds
        [out[0:n], out[4*n:5*n], out[8*n:9*n], out[12*n:13*n]] = quarter_round(eng, out[0:n], out[4*n:5*n], out[8*n:9*n], out[12*n:13*n])
        [out[n: 2*n], out[5*n:6*n], out[9*n:10*n], out[13*n:14*n]] = quarter_round(eng, out[n: 2*n], out[5*n:6*n], out[9*n:10*n], out[13*n:14*n])
        [out[2*n:3*n], out[6*n:7*n], out[10*n:11*n], out[14*n:15*n]] = quarter_round(eng, out[2*n:3*n], out[6*n:7*n], out[10*n:11*n], out[14*n:15*n])
        [out[3*n:4*n], out[7*n:8*n], out[11*n:12*n], out[15*n:16*n]] = quarter_round(eng, out[3*n:4*n], out[7*n:8*n], out[11*n:12*n], out[15*n:16*n])
        
        # Diagonal rounds
        [out[0:n], out[5*n:6*n], out[10*n:11*n], out[15*n:16*n]] = quarter_round(eng, out[0:n], out[5*n:6*n], out[10*n:11*n], out[15*n:16*n])
        [out[1*n:2*n], out[6*n:7*n], out[11*n:12*n], out[12*n:13*n]] = quarter_round(eng, out[1*n:2*n], out[6*n:7*n], out[11*n:12*n], out[12*n:13*n])
        [out[2*n:3*n], out[7*n:8*n], out[8*n:9*n], out[13*n:14*n]] = quarter_round(eng, out[2*n:3*n], out[7*n:8*n], out[8*n:9*n], out[13*n:14*n])
        [out[3*n:4*n], out[4*n:5*n], out[9*n:10*n], out[14*n:15*n]] = quarter_round(eng, out[3*n:4*n], out[4*n:5*n], out[9*n:10*n], out[14*n:15*n])
    
    # Mod_Add the initial state to the output
    Mod_add(eng,const, out[0: 4*n], 4*n)
    Mod_add(eng,key, out[4*n: 12*n], 8*n)
    Mod_add(eng,count, out[12*n: 13*n], n)
    Mod_add(eng,nonce, out[13*n: 16*n], 3*n)
    
    return out

def ChaCha20_QFT(eng, key, const, count, nonce, out, n):
    
    for i in range(n * 16):
        if i < n * 4:
            CNOT | (const[i], out[i])
        elif  n * 4 <= i < n * 12:
            CNOT | (key[i - n * 4], out[i])
        elif n * 12 <= i < n * 13:
            CNOT | (count[i - n * 12], out[i])
        elif n * 13 <= i < n * 16:
            CNOT | (nonce[i - n * 13], out[i])
    
    for i in range(10):
        # Column rounds
        [out[0:n], out[4*n:5*n], out[8*n:9*n], out[12*n:13*n]] = quarter_round_QFT(eng, out[0:n], out[4*n:5*n], out[8*n:9*n], out[12*n:13*n])
        [out[n: 2*n], out[5*n:6*n], out[9*n:10*n], out[13*n:14*n]] = quarter_round_QFT(eng, out[n: 2*n], out[5*n:6*n], out[9*n:10*n], out[13*n:14*n])
        [out[2*n:3*n], out[6*n:7*n], out[10*n:11*n], out[14*n:15*n]] = quarter_round_QFT(eng, out[2*n:3*n], out[6*n:7*n], out[10*n:11*n], out[14*n:15*n])
        [out[3*n:4*n], out[7*n:8*n], out[11*n:12*n], out[15*n:16*n]] = quarter_round_QFT(eng, out[3*n:4*n], out[7*n:8*n], out[11*n:12*n], out[15*n:16*n])
        
        # Diagonal rounds
        [out[0:n], out[5*n:6*n], out[10*n:11*n], out[15*n:16*n]] = quarter_round_QFT(eng, out[0:n], out[5*n:6*n], out[10*n:11*n], out[15*n:16*n])
        [out[1*n:2*n], out[6*n:7*n], out[11*n:12*n], out[12*n:13*n]] = quarter_round_QFT(eng, out[1*n:2*n], out[6*n:7*n], out[11*n:12*n], out[12*n:13*n])
        [out[2*n:3*n], out[7*n:8*n], out[8*n:9*n], out[13*n:14*n]] = quarter_round_QFT(eng, out[2*n:3*n], out[7*n:8*n], out[8*n:9*n], out[13*n:14*n])
        [out[3*n:4*n], out[4*n:5*n], out[9*n:10*n], out[14*n:15*n]] = quarter_round_QFT(eng, out[3*n:4*n], out[4*n:5*n], out[9*n:10*n], out[14*n:15*n])
    
    # Mod_Add the initial state to the output
    Mod_add_QFT(eng,const, out[0: 4*n], 4*n)
    Mod_add_QFT(eng,key, out[4*n: 12*n], 8*n)
    Mod_add_QFT(eng,count, out[12*n: 13*n], n)
    Mod_add_QFT(eng,nonce, out[13*n: 16*n], 3*n)
    
    return out

def print_hex(eng, qubits):
    for i in reversed(range(8)):
        temp = 0
        temp = temp+int(qubits[4*i+3])*8
        temp = temp+int(qubits[4*i+2])*4
        temp = temp+int(qubits[4*i+1])*2
        temp = temp+int(qubits[4*i])

        temp = hex(temp)
        y = temp.replace("0x", "")
        print(y, end='')

def print_state(eng, qubits):
    for i in reversed(range(128)):
        temp = 0
        temp = temp+int(qubits[4*i+3])*8
        temp = temp+int(qubits[4*i+2])*4
        temp = temp+int(qubits[4*i+1])*2
        temp = temp+int(qubits[4*i])

        temp = hex(temp)
        y = temp.replace("0x", "")
        print(y, end='')
        if(i%8 == 0 and i > 0):
            print(' 0x', end='')
        
        
def to_qreg(eng, a, data, bit):
    for i in range(bit):
        if (data >> i & 1):
            X | a[i]

def Run(eng):
    key = eng.allocate_qureg(256)
    const = eng.allocate_qureg(128)
    count = eng.allocate_qureg(32)
    nonce = eng.allocate_qureg(96)
    out = eng.allocate_qureg(512)
    
    if (not resource_check):
        to_qreg(eng, count, 0x1, 32)
        to_qreg(eng, const, 0x0, 128)
    
    out = ChaCha20(eng, key, const, count, nonce, out, 32)
    All(Measure) | out
    
    if (not resource_check):
        print('Ciphertext : 0x', end='')
        print_state(eng, out)

def get_circuit_depth(drawer):
    # Lấy danh sách các thao tác trên từng qubit
    qubit_lines = drawer._qubit_lines  # dict: qubit_id -> list of CircuitItem
    # Đếm số thao tác lớn nhất trên một qubit (tức là depth)
    return max(len(cmds) for cmds in qubit_lines.values()) if qubit_lines else 0

if __name__ == "__main__":
    
    global resource_check
    
    print("Running ChaCha20 on ProjectQ...\n")
    Simulator = ClassicalSimulator()
    eng = MainEngine(Simulator)
    resource_check = 0
    Run(eng)
    print('\n')

    print('Estimate cost...')
    Resource = ResourceCounter()
    eng = MainEngine(Resource)
    resource_check = 1
    Run(eng)
    print(Resource)
    print('\n')
    
    drawer = CircuitDrawer()
    eng = MainEngine(drawer)
    resource_check = 1
    Run(eng)
    eng.flush()
    depth = get_circuit_depth(drawer)
    print("Circuit depth (layers):", depth)
    print('\n')
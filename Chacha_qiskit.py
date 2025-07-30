from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit import QuantumCircuit, transpile
from qiskit.visualization import plot_histogram
from qiskit_aer import  Aer, AerSimulator
import matplotlib.pyplot as plt
from math import pi
from qiskit.circuit.library import QFT, PermutationGate

# backend = Aer.get_backend('statevector_simulator')

n_qubits = 32

def mod_add(circ, b, a):
    qft = QFT(num_qubits= n_qubits, do_swaps=False).to_gate()
    iqft = qft.inverse()
    circ.append(qft, a)
    # circ.barrier()
    # Controlled phase rotations: add a to b in Fourier basis
    for i in range(n_qubits):
        for j in range(n_qubits - i):
            if i >= j:
                angle = 2 * pi / (2 ** (i - j + 1))
                circ.cp(angle, b[j], a[i])
    # circ.barrier()
    circ.append(iqft, a)
    # circ.barrier()
    
def cyclic_left_shift_12(circ, qreg):
    for i in range(6):
        circ.swap(qreg[20+i], qreg[31-i])
    # circ.barrier()
    for i in range(10):
        circ.swap(qreg[0+i], qreg[19-i])
    # circ.barrier()
    for i in range(16):
        circ.swap(qreg[0+i], qreg[31-i])
    # circ.barrier()

def cyclic_left_shift_8(circ, qreg):
    # Bước 1: Đảo ngược 8 bit cao nhất (bit 24-31)
    for i in range(4):
        circ.swap(qreg[24+i], qreg[31-i])
    
    # Bước 2: Đảo ngược 24 bit thấp nhất (bit 0-23)
    for i in range(12):
        circ.swap(qreg[0+i], qreg[23-i])
    
    # Bước 3: Đảo ngược toàn bộ 32 bit
    for i in range(16):
        circ.swap(qreg[0+i], qreg[31-i])

def cyclic_left_shift_7(circ, qreg):
    # Bước 1: Đảo ngược 7 bit cao nhất (bit 25-31)
    for i in range(3):  # Chỉ cần 3 phép SWAP để đảo 7 bit
        circ.swap(qreg[25+i], qreg[31-i])
    
    # Bước 2: Đảo ngược 25 bit thấp nhất (bit 0-24)
    for i in range(12):
        circ.swap(qreg[0+i], qreg[24-i])
    
    # Bước 3: Đảo ngược toàn bộ 32 bit
    for i in range(16):
        circ.swap(qreg[0+i], qreg[31-i])

def quarter_round(circ, a, b, c, d):

    mod_add(circ, b, a)
    circ.barrier()
    for i in range(n_qubits):
        circ.cx(a[i], d[i])
    circ.barrier()
    for i in range(16):
        circ.swap(d[i], d[i + 16])
    circ.barrier()

    mod_add(circ, d, c)
    circ.barrier()
    for i in range(n_qubits):
        circ.cx(c[i], b[i])
    circ.barrier()
    cyclic_left_shift_12(circ, b)
    circ.barrier()
    
    mod_add(circ, b, a)
    circ.barrier()
    for i in range(n_qubits):
        circ.cx(a[i], d[i])
    circ.barrier()
    cyclic_left_shift_8(circ, d)
    circ.barrier()
    
    mod_add(circ, d, c)
    circ.barrier()
    for i in range(n_qubits):
        circ.cx(c[i], b[i])
    circ.barrier()
    cyclic_left_shift_7(circ, b)
    circ.barrier()

def ChaCha20(circ, key, const, count, nonce, out, cl_out):
    for i in range(n_qubits * 16):
        if i < n_qubits * 4:
            circ.cx(const[i], out[i])
        elif n_qubits * 4 <= i < n_qubits * 12:
            circ.cx(key[i - n_qubits * 4], out[i])
        elif n_qubits * 12 <= i < n_qubits * 13:
            circ.cx(count[i - n_qubits * 12], out[i])
        elif n_qubits * 13 <= i < n_qubits * 16:
            circ.cx(nonce[i - n_qubits * 13], out[i])
    circ.barrier()
    
    # Vòng 1/10
    quarter_round(circ, out[0:n_qubits], out[n_qubits:n_qubits * 2], out[n_qubits * 2:n_qubits * 3], out[n_qubits * 3:n_qubits * 4])
    quarter_round(circ, out[n_qubits * 4:n_qubits * 5], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 7:n_qubits * 8])
    quarter_round(circ, out[n_qubits * 8:n_qubits * 9], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 11:n_qubits * 12])
    quarter_round(circ, out[n_qubits * 12:n_qubits * 13], out[n_qubits * 13:n_qubits * 14], out[n_qubits * 14:n_qubits * 15], out[n_qubits * 15:n_qubits * 16])
    
    quarter_round(circ, out[0:n_qubits], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 15:n_qubits * 16])
    quarter_round(circ, out[n_qubits:n_qubits * 2], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 11:n_qubits * 12], out[n_qubits * 12:n_qubits * 13])
    quarter_round(circ, out[n_qubits * 2:n_qubits * 3], out[n_qubits * 7:n_qubits * 8], out[n_qubits * 8:n_qubits * 9], out[n_qubits * 13:n_qubits * 14])
    quarter_round(circ, out[n_qubits * 3:n_qubits * 4], out[n_qubits * 4:n_qubits * 5], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 14:n_qubits * 15])
    
    circ.barrier()
    
    # Vòng 2/10
    quarter_round(circ, out[0:n_qubits], out[n_qubits:n_qubits * 2], out[n_qubits * 2:n_qubits * 3], out[n_qubits * 3:n_qubits * 4])
    quarter_round(circ, out[n_qubits * 4:n_qubits * 5], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 7:n_qubits * 8])
    quarter_round(circ, out[n_qubits * 8:n_qubits * 9], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 11:n_qubits * 12])
    quarter_round(circ, out[n_qubits * 12:n_qubits * 13], out[n_qubits * 13:n_qubits * 14], out[n_qubits * 14:n_qubits * 15], out[n_qubits * 15:n_qubits * 16])
    
    quarter_round(circ, out[0:n_qubits], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 15:n_qubits * 16])
    quarter_round(circ, out[n_qubits:n_qubits * 2], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 11:n_qubits * 12], out[n_qubits * 12:n_qubits * 13])
    quarter_round(circ, out[n_qubits * 2:n_qubits * 3], out[n_qubits * 7:n_qubits * 8], out[n_qubits * 8:n_qubits * 9], out[n_qubits * 13:n_qubits * 14])
    quarter_round(circ, out[n_qubits * 3:n_qubits * 4], out[n_qubits * 4:n_qubits * 5], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 14:n_qubits * 15])
    
    circ.barrier()
    
    # Vòng 3/10
    quarter_round(circ, out[0:n_qubits], out[n_qubits:n_qubits * 2], out[n_qubits * 2:n_qubits * 3], out[n_qubits * 3:n_qubits * 4])
    quarter_round(circ, out[n_qubits * 4:n_qubits * 5], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 7:n_qubits * 8])
    quarter_round(circ, out[n_qubits * 8:n_qubits * 9], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 11:n_qubits * 12])
    quarter_round(circ, out[n_qubits * 12:n_qubits * 13], out[n_qubits * 13:n_qubits * 14], out[n_qubits * 14:n_qubits * 15], out[n_qubits * 15:n_qubits * 16])
    
    quarter_round(circ, out[0:n_qubits], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 15:n_qubits * 16])
    quarter_round(circ, out[n_qubits:n_qubits * 2], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 11:n_qubits * 12], out[n_qubits * 12:n_qubits * 13])
    quarter_round(circ, out[n_qubits * 2:n_qubits * 3], out[n_qubits * 7:n_qubits * 8], out[n_qubits * 8:n_qubits * 9], out[n_qubits * 13:n_qubits * 14])
    quarter_round(circ, out[n_qubits * 3:n_qubits * 4], out[n_qubits * 4:n_qubits * 5], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 14:n_qubits * 15])
    
    circ.barrier()
    
    # Vòng 4/10
    quarter_round(circ, out[0:n_qubits], out[n_qubits:n_qubits * 2], out[n_qubits * 2:n_qubits * 3], out[n_qubits * 3:n_qubits * 4])
    quarter_round(circ, out[n_qubits * 4:n_qubits * 5], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 7:n_qubits * 8])
    quarter_round(circ, out[n_qubits * 8:n_qubits * 9], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 11:n_qubits * 12])
    quarter_round(circ, out[n_qubits * 12:n_qubits * 13], out[n_qubits * 13:n_qubits * 14], out[n_qubits * 14:n_qubits * 15], out[n_qubits * 15:n_qubits * 16])
    
    quarter_round(circ, out[0:n_qubits], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 15:n_qubits * 16])
    quarter_round(circ, out[n_qubits:n_qubits * 2], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 11:n_qubits * 12], out[n_qubits * 12:n_qubits * 13])
    quarter_round(circ, out[n_qubits * 2:n_qubits * 3], out[n_qubits * 7:n_qubits * 8], out[n_qubits * 8:n_qubits * 9], out[n_qubits * 13:n_qubits * 14])
    quarter_round(circ, out[n_qubits * 3:n_qubits * 4], out[n_qubits * 4:n_qubits * 5], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 14:n_qubits * 15])
    
    circ.barrier()
    
    #vòng 5/10
    quarter_round(circ, out[0:n_qubits], out[n_qubits:n_qubits * 2], out[n_qubits * 2:n_qubits * 3], out[n_qubits * 3:n_qubits * 4])
    quarter_round(circ, out[n_qubits * 4:n_qubits * 5], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 7:n_qubits * 8])
    quarter_round(circ, out[n_qubits * 8:n_qubits * 9], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 11:n_qubits * 12])
    quarter_round(circ, out[n_qubits * 12:n_qubits * 13], out[n_qubits * 13:n_qubits * 14], out[n_qubits * 14:n_qubits * 15], out[n_qubits * 15:n_qubits * 16])
    
    quarter_round(circ, out[0:n_qubits], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 15:n_qubits * 16])
    quarter_round(circ, out[n_qubits:n_qubits * 2], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 11:n_qubits * 12], out[n_qubits * 12:n_qubits * 13])
    quarter_round(circ, out[n_qubits * 2:n_qubits * 3], out[n_qubits * 7:n_qubits * 8], out[n_qubits * 8:n_qubits * 9], out[n_qubits * 13:n_qubits * 14])
    quarter_round(circ, out[n_qubits * 3:n_qubits * 4], out[n_qubits * 4:n_qubits * 5], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 14:n_qubits * 15])
    
    circ.barrier()
    
    # Vòng 6/10
    quarter_round(circ, out[0:n_qubits], out[n_qubits:n_qubits * 2], out[n_qubits * 2:n_qubits * 3], out[n_qubits * 3:n_qubits * 4])
    quarter_round(circ, out[n_qubits * 4:n_qubits * 5], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 7:n_qubits * 8])
    quarter_round(circ, out[n_qubits * 8:n_qubits * 9], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 11:n_qubits * 12])
    quarter_round(circ, out[n_qubits * 12:n_qubits * 13], out[n_qubits * 13:n_qubits * 14], out[n_qubits * 14:n_qubits * 15], out[n_qubits * 15:n_qubits * 16])
    
    quarter_round(circ, out[0:n_qubits], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 15:n_qubits * 16])
    quarter_round(circ, out[n_qubits:n_qubits * 2], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 11:n_qubits * 12], out[n_qubits * 12:n_qubits * 13])
    quarter_round(circ, out[n_qubits * 2:n_qubits * 3], out[n_qubits * 7:n_qubits * 8], out[n_qubits * 8:n_qubits * 9], out[n_qubits * 13:n_qubits * 14])
    quarter_round(circ, out[n_qubits * 3:n_qubits * 4], out[n_qubits * 4:n_qubits * 5], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 14:n_qubits * 15])
    
    circ.barrier()
    
    # Vòng 7/10
    quarter_round(circ, out[0:n_qubits], out[n_qubits:n_qubits * 2], out[n_qubits * 2:n_qubits * 3], out[n_qubits * 3:n_qubits * 4])
    quarter_round(circ, out[n_qubits * 4:n_qubits * 5], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 7:n_qubits * 8])
    quarter_round(circ, out[n_qubits * 8:n_qubits * 9], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 11:n_qubits * 12])
    quarter_round(circ, out[n_qubits * 12:n_qubits * 13], out[n_qubits * 13:n_qubits * 14], out[n_qubits * 14:n_qubits * 15], out[n_qubits * 15:n_qubits * 16])
    
    quarter_round(circ, out[0:n_qubits], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 15:n_qubits * 16])
    quarter_round(circ, out[n_qubits:n_qubits * 2], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 11:n_qubits * 12], out[n_qubits * 12:n_qubits * 13])
    quarter_round(circ, out[n_qubits * 2:n_qubits * 3], out[n_qubits * 7:n_qubits * 8], out[n_qubits * 8:n_qubits * 9], out[n_qubits * 13:n_qubits * 14])
    quarter_round(circ, out[n_qubits * 3:n_qubits * 4], out[n_qubits * 4:n_qubits * 5], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 14:n_qubits * 15])
    
    circ.barrier()
    
    # Vòng 8/10
    quarter_round(circ, out[0:n_qubits], out[n_qubits:n_qubits * 2], out[n_qubits * 2:n_qubits * 3], out[n_qubits * 3:n_qubits * 4])
    quarter_round(circ, out[n_qubits * 4:n_qubits * 5], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 7:n_qubits * 8])
    quarter_round(circ, out[n_qubits * 8:n_qubits * 9], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 11:n_qubits * 12])
    quarter_round(circ, out[n_qubits * 12:n_qubits * 13], out[n_qubits * 13:n_qubits * 14], out[n_qubits * 14:n_qubits * 15], out[n_qubits * 15:n_qubits * 16])
    
    quarter_round(circ, out[0:n_qubits], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 15:n_qubits * 16])
    quarter_round(circ, out[n_qubits:n_qubits * 2], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 11:n_qubits * 12], out[n_qubits * 12:n_qubits * 13])
    quarter_round(circ, out[n_qubits * 2:n_qubits * 3], out[n_qubits * 7:n_qubits * 8], out[n_qubits * 8:n_qubits * 9], out[n_qubits * 13:n_qubits * 14])
    quarter_round(circ, out[n_qubits * 3:n_qubits * 4], out[n_qubits * 4:n_qubits * 5], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 14:n_qubits * 15])
    
    circ.barrier()
    
    # Vògn 9/10
    quarter_round(circ, out[0:n_qubits], out[n_qubits:n_qubits * 2], out[n_qubits * 2:n_qubits * 3], out[n_qubits * 3:n_qubits * 4])
    quarter_round(circ, out[n_qubits * 4:n_qubits * 5], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 7:n_qubits * 8])
    quarter_round(circ, out[n_qubits * 8:n_qubits * 9], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 11:n_qubits * 12])
    quarter_round(circ, out[n_qubits * 12:n_qubits * 13], out[n_qubits * 13:n_qubits * 14], out[n_qubits * 14:n_qubits * 15], out[n_qubits * 15:n_qubits * 16])
    
    quarter_round(circ, out[0:n_qubits], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 15:n_qubits * 16])
    quarter_round(circ, out[n_qubits:n_qubits * 2], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 11:n_qubits * 12], out[n_qubits * 12:n_qubits * 13])
    quarter_round(circ, out[n_qubits * 2:n_qubits * 3], out[n_qubits * 7:n_qubits * 8], out[n_qubits * 8:n_qubits * 9], out[n_qubits * 13:n_qubits * 14])
    quarter_round(circ, out[n_qubits * 3:n_qubits * 4], out[n_qubits * 4:n_qubits * 5], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 14:n_qubits * 15])
    
    circ.barrier()
    
    # Vòng 10/10
    quarter_round(circ, out[0:n_qubits], out[n_qubits:n_qubits * 2], out[n_qubits * 2:n_qubits * 3], out[n_qubits * 3:n_qubits * 4])
    quarter_round(circ, out[n_qubits * 4:n_qubits * 5], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 7:n_qubits * 8])
    quarter_round(circ, out[n_qubits * 8:n_qubits * 9], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 11:n_qubits * 12])
    quarter_round(circ, out[n_qubits * 12:n_qubits * 13], out[n_qubits * 13:n_qubits * 14], out[n_qubits * 14:n_qubits * 15], out[n_qubits * 15:n_qubits * 16])
    
    quarter_round(circ, out[0:n_qubits], out[n_qubits * 5:n_qubits * 6], out[n_qubits * 10:n_qubits * 11], out[n_qubits * 15:n_qubits * 16])
    quarter_round(circ, out[n_qubits:n_qubits * 2], out[n_qubits * 6:n_qubits * 7], out[n_qubits * 11:n_qubits * 12], out[n_qubits * 12:n_qubits * 13])
    quarter_round(circ, out[n_qubits * 2:n_qubits * 3], out[n_qubits * 7:n_qubits * 8], out[n_qubits * 8:n_qubits * 9], out[n_qubits * 13:n_qubits * 14])
    quarter_round(circ, out[n_qubits * 3:n_qubits * 4], out[n_qubits * 4:n_qubits * 5], out[n_qubits * 9:n_qubits * 10], out[n_qubits * 14:n_qubits * 15])
    
    circ.barrier()
    
    for i in range(n_qubits * 16):
        if i < n_qubits * 4:
            circ.cx(const[i], out[i])
        elif n_qubits * 4 <= i < n_qubits * 12:
            circ.cx(key[i - n_qubits * 4], out[i])
        elif n_qubits * 12 <= i < n_qubits * 13:
            circ.cx(count[i - n_qubits * 12], out[i])
        elif n_qubits * 13 <= i < n_qubits * 16:
            circ.cx(nonce[i - n_qubits * 13], out[i])
    
    circ.barrier()
    
    circ.measure(out, cl_out)

if __name__ == "__main__":
    key = QuantumRegister(n_qubits * 8, 'key')
    const = QuantumRegister(n_qubits * 4, 'const')
    count = QuantumRegister(n_qubits, 'count')
    nonce = QuantumRegister(n_qubits * 3, 'nonce')
    out = QuantumRegister(n_qubits * 16, 'out')
    b = QuantumRegister(n_qubits, 'b')
    
    cl_out = ClassicalRegister(n_qubits, 'cl_out')
    
    circ = QuantumCircuit(key, const, count, nonce, out, b, cl_out)
    # circ.x(count[0])    
    # mod_add(circ, count, b)
    # circ.measure(b, cl_out)
    
    # ChaCha20(circ, key, const, count, nonce, out, cl_out)
    # for i in range(n_qubits * 16):
    #     if i < n_qubits * 4:
    #         circ.cx(const[i], out[i])
    #     elif n_qubits * 4 <= i < n_qubits * 12:
    #         circ.cx(key[i - n_qubits * 4], out[i])
    #     elif n_qubits * 12 <= i < n_qubits * 13:
    #         circ.cx(count[i - n_qubits * 12], out[i])
    #     elif n_qubits * 13 <= i < n_qubits * 16:
    #         circ.cx(nonce[i - n_qubits * 13], out[i])
    # circ.barrier()
    
    # Vòng 1/10
    quarter_round(circ, out[0:n_qubits], out[n_qubits:n_qubits * 2], out[n_qubits * 2:n_qubits * 3], out[n_qubits * 3:n_qubits * 4])
    # circ_compiled = transpile(circ, backend=backend, optimization_level=3)
    # print(circ.draw())
    circ.draw('mpl')
    plt.show()
    # job = backend.run(circ_compiled)
    # result = job.result()
    # counts = result.get_counts(circ_compiled)
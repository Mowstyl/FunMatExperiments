import numpy as np
import scipy.sparse as sp
import doki
import re
import faulthandler

from qsimov import SimpleGate, Funmatrix
from mpi4py import MPI


comm = MPI.COMM_WORLD
faulthandler.enable()
_h_regex = re.compile(r"^H[0-9]*$")
gates = {}

def root_print(*args, **kwargs):
    if comm.rank == 0:
        print(*args, **kwargs)

def get_gate(gate):
    gate_mat = None
    
    if gate in gates:
        return gates[gate]
    if not _h_regex.match(gate):
        raw_sgate = SimpleGate(gate)
        shape = raw_sgate.matrix.shape
        matrix = [[complex(raw_sgate.matrix[i, j]) for j in range(shape[1])] for i in range(shape[0])]
        gate_mat = doki.funmatrix_create(matrix, False)
    else:
        nq = 1
        if len(gate[1:]) > 0:
            nq = int(gate[1:])
        if nq < 1:
            raise ValueError("Hadamard gate for less than one qubit is impossible")
        gate_mat = doki.funmatrix_hadamard(nq, False)
    gates[gate] = gate_mat
    return gate_mat


def get_identity(nq):
    aux = doki.funmatrix_identity(nq, False)
    return aux


def get_system(nq):
    rho = doki.funmatrix_densityzero(nq, False)
    return rho


def get_projector(nq, qubit, value):
    proj = None
    if not bool(value):
        proj = get_system(nq)
    else:
        proj = doki.funmatrix_create([[0, 0], [0, 1]], False)
        left, right = 0, 0
        if qubit > 0:
            left = qubit
        if qubit < nq-1:
            right = nq-qubit-1
        if left != 0 or right != 0:
            proj = doki.funmatrix_eyekron(proj, left, right, False)
    return proj


def apply_gate(nq, sys, gate, qubit, num_controls=0, verbose=False):
    newgate = get_gate(gate)
    nqg = np.log2(doki.funmatrix_shape(newgate, False)[0])
    if nqg % 1 != 0:
        raise ValueError(f"Wrong gate size: 2^{nqg}")
    nqg = int(nqg)
    while num_controls > 0:
        newgate = doki.funmatrix_addcontrol(newgate, verbose)
        num_controls -= 1
        nqg += 1
    left, right = 0, 0
    if qubit > 0:
        left = qubit
    if qubit < nq-nqg:
        right = nq-qubit-nqg
    if left != 0 or right != 0:
        newgate = doki.funmatrix_eyekron(newgate, left, right, False)
    gatetrans = doki.funmatrix_dagger(newgate, verbose)
    res = doki.funmatrix_matmul(newgate, sys, verbose)
    res = doki.funmatrix_matmul(res, gatetrans, verbose)
    return res

def trace(mat, size, verbose=False):
    my_trace = np.zeros(1)
    for i in range(comm.rank, size, comm.size):
        my_trace += doki.funmatrix_get(mat, i, i, verbose).real
    trace = np.zeros(1)
    comm.Reduce(
        my_trace,
        trace,
        op = MPI.SUM,
        root = 0
    )
    comm.Bcast(trace, root=0)
    return trace[0]

def measure(nq, sys, qubit, verbose=False):
    P = get_projector(nq, qubit, 1)
    P_rho = doki.funmatrix_matmul(P, sys, verbose)
    size = 2**nq
    p = trace(P_rho, size, verbose)
    # p = doki.funmatrix_trace(P_rho, verbose).real
    roll = np.empty(1)
    if comm.rank == 0:
        roll = np.random.rand(1)
    comm.Bcast(roll, root=0)
    res = roll[0] < p
    if not res:
        P = get_projector(nq, qubit, 0)
        P_rho = doki.funmatrix_matmul(P, sys, verbose)
        p = trace(P_rho, size, verbose)
        #p = doki.funmatrix_trace(P_rho, verbose).real
    numerator = doki.funmatrix_matmul(P_rho, P, verbose)
    aux = doki.funmatrix_scalar_div(numerator, p, verbose)
    return (res, aux)


def main():
    root_print("Starting!")
    sys = get_system(2)
    root_print("rho0 =", Funmatrix(sys, "rho0")[:])
    sys2 = apply_gate(2, sys, "H", 0)
    root_print("rho1 =", Funmatrix(sys2, "rho1")[:])
    del sys
    sys3 = apply_gate(2, sys2, "X", 0, num_controls=1)
    root_print("rho2 =", Funmatrix(sys3, "rho2")[:])
    del sys2
    res, sys4 = measure(2, sys3, 0)
    root_print("res =", res)
    root_print("rho3 =", Funmatrix(sys4, "rho3")[:])
    del sys3
    del sys4

    sys = get_system(1)
    root_print("rho0 =", Funmatrix(sys, "rho0")[:])
    sys2 = apply_gate(1, sys, "H", 0)
    root_print("rho1 =", Funmatrix(sys2, "rho1")[:])
    del sys
    res, sys3 = measure(1, sys2, 0)
    del sys2
    root_print("res =", res)
    root_print("rho2 =", Funmatrix(sys3, "rho2")[:])
    #del sys3


if __name__ == "__main__":
    main()

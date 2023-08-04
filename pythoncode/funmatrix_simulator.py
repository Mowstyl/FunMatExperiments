import doki
import re
import numpy as np
#import faulthandler

from qsimov import SimpleGate, Funmatrix
from mpi4py import MPI


#faulthandler.enable()  # So that segfaults and the likes are shown

_h_regex = re.compile(r"^H[0-9]*$")
gates = {}
comm = MPI.COMM_WORLD


def get_gate(gate):
    gate_mat = None
    
    if gate in gates:
        return gates[gate]
    if not _h_regex.match(gate):
        raw_sgate = SimpleGate(gate)
        shape = raw_sgate.matrix.shape
        matrix = [[complex(raw_sgate.matrix[i, j]) for j in range(shape[1])] for i in range(shape[0])]
        gate_mat = doki.funmatrix_create(matrix, False)
        del matrix
    else:
        nq = 1
        if len(gate[1:]) > 0:
            nq = int(gate[1:])
        if nq < 1:
            raise ValueError("Hadamard gate for less than one qubit is impossible")
        gate_mat = doki.funmatrix_hadamard(nq, False)
    gates[gate] = gate_mat
    return gate_mat


def get_system(nq):
    rho = doki.funmatrix_statezero(nq, False)
    return rho, nq


def apply_gate(sys, gate, targets, controls=set(), anticontrols=set(), verbose=False):
    sys, nq = sys
    newgate = get_gate(gate)
    if type(targets) != list:
        targets = [targets]
    if type(controls) != set:
        controls = set([controls])
    if type(anticontrols) != set:
        anticontrols = set([anticontrols])
    res = doki.funmatrix_apply(sys, newgate, targets, controls, anticontrols, verbose)
    return res, nq

def sqsum(mat, size, qubit, value=1, verbose=False):
    my_sum = np.zeros(1, dtype=float)
    powqb = 1 << qubit
    if value == 1:
        i = powqb
    count = 0
    for i in range(value * powqb, size, powqb << 1):
        for j in range(i, i + powqb):
            if count % comm.size == comm.rank:
                my_sum[0] += (doki.funmatrix_get(mat, j, 0, verbose)**2).real
            count += 1
    tot_sum = np.empty(1, dtype=float)
    comm.Reduce(
        my_sum,
        tot_sum,
        op = MPI.SUM,
        root = 0
    )
    comm.Bcast(tot_sum, root=0)
    return tot_sum[0]


def measure(sys, qubit, verbose=False):
    sys, nq = sys
    size = 2**nq
    p = sqsum(sys, size, qubit, value=1)
    roll = np.empty(1)
    if comm.rank == 0:
        roll = np.random.rand(1)
    comm.Bcast(roll, root=0)
    res = roll[0] < p
    del roll
    value = int(res)
    if not res:
        p = sqsum(sys, size, qubit, value=0)
    P_sys = doki.funmatrix_projection(sys, qubit, value, verbose)
    aux = doki.funmatrix_scalar_div(P_sys, np.sqrt(p), verbose)
    return res, (aux, nq)



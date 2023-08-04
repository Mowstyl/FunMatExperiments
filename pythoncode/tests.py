import tracemalloc
import time as t
import pandas as pd
import mpmath as mp
import numpy as np
import sympy as sp
#import faulthandler
import os
import random as rnd
import time as t
import pandas as pd

from funmatrix_simulator import apply_gate, measure, get_system, gates
from mpi4py import MPI


#faulthandler.enable()  # So that segfaults and the likes are shown
comm = MPI.COMM_WORLD
filters = (
    tracemalloc.Filter(inclusive=False, filename_pattern="<frozen importlib._bootstrap>"),
    tracemalloc.Filter(inclusive=False, filename_pattern="<frozen importlib._bootstrap_external>"),
    tracemalloc.Filter(inclusive=False, filename_pattern="<frozen _collections_abc>"),
    tracemalloc.Filter(inclusive=False, filename_pattern="<unknown>"),
    tracemalloc.Filter(inclusive=False, filename_pattern=tracemalloc.__file__),
    tracemalloc.Filter(inclusive=False, filename_pattern=os.path.join(os.path.dirname(sp.__file__), "*")),
    tracemalloc.Filter(inclusive=False, filename_pattern=os.path.join(os.path.dirname(np.__file__), "*")),
    tracemalloc.Filter(inclusive=False, filename_pattern=os.path.join(os.path.dirname(mp.__file__), "*"))
    )


def root_print(*args, **kwargs):
    if comm.rank == 0:
        print(*args, **kwargs, flush=True)


def grover2(sys):
    sys = apply_gate(sys, "X", targets=1, anticontrols=0)
    return sys


def grover3XOR(sys):
    sys = apply_gate(sys, "X", targets=2, controls=0)
    sys = apply_gate(sys, "X", targets=2, controls=1)
    return sys


def grover3OR(sys):
    sys = apply_gate(sys, "X", targets=2, controls=0)
    sys = apply_gate(sys, "X", targets=2, anticontrols=1)
    sys = apply_gate(sys, "X", targets=2, controls=0, anticontrols=1)
    return sys


def grover4(sys):
    sys = apply_gate(sys, "X", targets=3, controls=0, anticontrols=1)
    sys = apply_gate(sys, "X", targets=3, controls=2, anticontrols=0)
    return sys


def grover5(sys):
    sys = apply_gate(sys, "X", targets=4, controls=0, anticontrols=1)
    sys = apply_gate(sys, "X", targets=3, controls=2, anticontrols=4)
    sys = apply_gate(sys, "X", targets=4, controls=0, anticontrols=1)
    sys = apply_gate(sys, "X", targets=3)
    return sys


def grover6(sys):
    sys = apply_gate(sys, "X", targets=5, controls={0, 2}, anticontrols=1)
    sys = apply_gate(sys, "X", targets=4, anticontrols={3, 5})
    sys = apply_gate(sys, "X", targets=5, controls={0, 2}, anticontrols=1)
    sys = apply_gate(sys, "X", targets=4)
    return sys


def grover7(sys):
    sys = apply_gate(sys, "X", targets=4, controls=0, anticontrols=1)
    sys = apply_gate(sys, "X", targets=5, controls=0, anticontrols=2)
    sys = apply_gate(sys, "X", targets=6, controls=1, anticontrols=2)
    sys = apply_gate(sys, "X", targets=3, anticontrols={4, 5, 6})
    sys = apply_gate(sys, "X", targets=6, controls=1, anticontrols=2)
    sys = apply_gate(sys, "X", targets=5, controls=0, anticontrols=2)
    sys = apply_gate(sys, "X", targets=4, controls=0, anticontrols=1)
    sys = apply_gate(sys, "X", targets=3)
    return sys


def grover8(sys):
    sys = apply_gate(sys, "X", targets=5, controls=0, anticontrols=2)
    sys = apply_gate(sys, "X", targets=6, controls=3, anticontrols=2)
    sys = apply_gate(sys, "X", targets=7, controls=0, anticontrols=1)
    sys = apply_gate(sys, "X", targets=4, anticontrols={5, 6, 7})
    sys = apply_gate(sys, "X", targets=7, controls=0, anticontrols=1)
    sys = apply_gate(sys, "X", targets=6, controls=3, anticontrols=2)
    sys = apply_gate(sys, "X", targets=5, controls=0, anticontrols=2)
    sys = apply_gate(sys, "X", targets=4)
    return sys


def grover9(sys):
    sys = apply_gate(sys, "X", targets=5, controls=2, anticontrols=0)
    sys = apply_gate(sys, "X", targets=6, controls=2, anticontrols=3)
    sys = apply_gate(sys, "X", targets=7, anticontrols={5, 6})
    sys = apply_gate(sys, "X", targets=8, controls={0, 1})
    sys = apply_gate(sys, "X", targets=4, anticontrols={7, 8})
    sys = apply_gate(sys, "X", targets=8, controls={0, 1})
    sys = apply_gate(sys, "X", targets=7, anticontrols={5, 6})
    sys = apply_gate(sys, "X", targets=6, controls=2, anticontrols=3)
    sys = apply_gate(sys, "X", targets=5, controls=2, anticontrols=0)
    sys = apply_gate(sys, "X", targets=4)
    return sys


def grover10(sys):
    sys = apply_gate(sys, "X", targets=6, controls={0, 2}, anticontrols=1)
    sys = apply_gate(sys, "X", targets=7, controls=1)
    sys = apply_gate(sys, "X", targets=7, controls=2)
    sys = apply_gate(sys, "X", targets=8, anticontrols={1, 4})
    sys = apply_gate(sys, "X", targets=9, controls={0, 3, 4, 7}, anticontrols={6, 8})
    sys = apply_gate(sys, "X", targets=5, controls=2, anticontrols=9)
    sys = apply_gate(sys, "X", targets=9, controls={0, 3, 4, 7}, anticontrols={6, 8})
    sys = apply_gate(sys, "X", targets=8, anticontrols={1, 4})
    sys = apply_gate(sys, "X", targets=7, controls=2)
    sys = apply_gate(sys, "X", targets=7, controls=1)
    sys = apply_gate(sys, "X", targets=6, controls={0, 2}, anticontrols=1)
    return sys


def IAM(sys, nq, last=False):
    sys = apply_gate(sys, f"H{nq-1}", targets=[i for i in range(nq-1)])
    last_id = nq - 1
    sys = apply_gate(sys, "X", targets=last_id)
    sys = apply_gate(sys, "Z", targets=last_id, anticontrols={i for i in range(last_id)})
    if not last:
        sys = apply_gate(sys, "X", targets=last_id)
    sys = apply_gate(sys, f"H{nq-1}", targets=[i for i in range(nq-1)])

    return sys


oracles = [grover2, grover3OR, grover4, grover5, grover6, grover7, grover8, grover9, grover10]
no_ancilla = [2, 3, 4, 4, 5, 4, 5, 5, 6]


def grover(nq):
    if nq < 2 or nq > 10:
        raise ValueError("That test has not been designed")
    UwU_function = oracles[nq-2]
    not_ancilla = no_ancilla[nq-2]
    start_snap = tracemalloc.take_snapshot()
    start_time = t.time()
    sys = get_system(nq + 1)
    sys = apply_gate(sys, f"H{not_ancilla}", targets=[i for i in range(not_ancilla)])
    sys = UwU_function(sys)
    sys = IAM(sys, not_ancilla, last=True)
    results = []
    for i in range(not_ancilla - 1):
        res, sys = measure(sys, i)
        results.append(res)
    stop_time = t.time()
    stop_snap = tracemalloc.take_snapshot()
    start_stats = start_snap.filter_traces(filters).statistics("lineno")
    stop_stats = stop_snap.filter_traces(filters).statistics("lineno")
    start_size = sum(stat.size for stat in start_stats)
    stop_size = sum(stat.size for stat in stop_stats)
    return (results, stop_time - start_time, stop_size - start_size)

def bell_pair():
    sys = get_system(2)
    sys2 = apply_gate(sys, "H", targets=0)
    del sys
    sys3 = apply_gate(sys2, "X", targets=0, controls=1)
    del sys2
    res, sys4 = measure(sys3, 0)
    return (res, sys4)


def deutsch_jozsa(nq, is_constant):
    start_snap = tracemalloc.take_snapshot()
    start_time = t.time()
    sys = get_system(nq + 1)
    sys = apply_gate(sys, "X", targets=nq)
    sys = apply_gate(sys, f"H{nq+1}", targets=[i for i in range(nq+1)])
    if is_constant: # f(x1, ..., xn) = True
        sys = apply_gate(sys, "X", targets=nq)
    else:           # f(x1, ..., xn) = xn
        sys = apply_gate(sys, "X", targets=nq, controls=nq-1)
    sys = apply_gate(sys, f"H{nq}", targets=[i for i in range(nq)])
    results = []
    for i in range(nq):
        res, sys = measure(sys, i)
        results.append(res)
    stop_time = t.time()
    stop_snap = tracemalloc.take_snapshot()
    start_stats = start_snap.filter_traces(filters).statistics("lineno")
    stop_stats = stop_snap.filter_traces(filters).statistics("lineno")
    start_size = sum(stat.size for stat in start_stats)
    stop_size = sum(stat.size for stat in stop_stats)
    return (results, stop_time - start_time, stop_size - start_size)


def bernstein_vazirani(bit_str):
    nq = len(bit_str)
    start_snap = tracemalloc.take_snapshot()
    start_time = t.time()
    sys = get_system(nq + 1)
    sys = apply_gate(sys, "X", targets=nq)
    sys = apply_gate(sys, f"H{nq+1}", targets=[i for i in range(nq+1)])
    for i in range(nq):
        if bit_str[i]:
            sys = apply_gate(sys, "X", targets=nq, controls=i)
    sys = apply_gate(sys, f"H{nq}", targets=[i for i in range(nq)])
    results = []
    for i in range(nq):
        res, sys = measure(sys, i)
        results.append(res)
    stop_time = t.time()
    stop_snap = tracemalloc.take_snapshot()
    start_stats = start_snap.filter_traces(filters).statistics("lineno")
    stop_stats = stop_snap.filter_traces(filters).statistics("lineno")
    start_size = sum(stat.size for stat in start_stats)
    stop_size = sum(stat.size for stat in stop_stats)
    return (results, stop_time - start_time, stop_size - start_size)


def main():
    tracemalloc.start()
    root_print("Warmup...")
    for i in range(10):
        _, _ = bell_pair()
        _, _, _ = deutsch_jozsa(1, True)
        _, _, _ = bernstein_vazirani([False])
        _, _, _ = grover(2)
    global gates
    gates.clear()
    root_print("Warm!")

    min_qb = 2
    max_qb = 10
    data = {"num_qubits": [], "DJ_times": [], "DJ_sizes": [], "BV_times": [], "BV_sizes": [], "grover_times": [], "grover_sizes": []}
    roll = None
    while True:
        for n in range(min_qb, max_qb + 1):
            root_print(f"{n} qubit tests")
            num_iters = 1
            '''
            if n < 8:
                num_iters += int(255/2**n)
            root_print(f"{num_iters} iterations...")
            dec = num_iters // 100
            if dec <= 1:
                dec = num_iters // 10
            dec = max(1, dec)
            point = 1
            root_print(f"ACK each {dec} iterations")
            '''
            for it in range(num_iters):
                if comm.rank == 0:
                    roll = np.random.randint(2, size=1)
                else:
                    roll = np.empty(1, dtype=int)
                comm.Bcast(roll, root=0)
                is_constant = bool(roll[0])
                res, time, size = deutsch_jozsa(n-1, is_constant)
                if any(res) == is_constant:
                    raise ValueError(f"Error when executing DJ with is_constant = {is_constant}, n = {n}")
                data["DJ_times"].append(time)
                data["DJ_sizes"].append(size)
                if comm.rank == 0:
                    roll = np.random.randint(2, size=n-1)
                else:
                    roll = np.empty(n-1, dtype=int)
                comm.Bcast(roll, root=0)
                bit_str = [bool(bit) for bit in roll]
                res, time, size = bernstein_vazirani(bit_str)
                if bit_str != res:
                    raise ValueError(f"Error when executing BV with n = {n}, bit_str = {bit_str}")
                data["BV_times"].append(time)
                data["BV_sizes"].append(size)
                if n > 1:
                    res, time, size = grover(n)
                    data["grover_times"].append(time)
                    data["grover_sizes"].append(size)
                else:
                    data["grover_times"].append(np.nan)
                    data["grover_sizes"].append(np.nan)
                data["num_qubits"].append(n)
                '''
                if point * dec == it:
                    root_print(f"{(100 * (it + 1)) // num_iters}%")
                    point += 1
                '''
        if comm.rank == 0:
            pd.DataFrame(data).to_csv(f"fmtests/raw/data_{int(t.time())}_{min_qb}_{max_qb}.csv", sep=';', index=False)


if __name__ == "__main__":
    main()

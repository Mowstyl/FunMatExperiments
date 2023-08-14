import time as t
#import faulthandler
import os
import random as rnd
import time as t
import numpy as np
import pandas as pd
import qsimov as qj
import doki


grover2 = qj.QGate(2, 0, "Uw2")
grover2.add_operation("X", targets=1, anticontrols=0)

grover3XOR = qj.QGate(3, 0, "UwXOR")
grover3XOR.add_operation("X", targets=2, controls=0)
grover3XOR.add_operation("X", targets=2, controls=1)

grover3OR = qj.QGate(3, 0, "UwOR")
grover3OR.add_operation("X", targets=2, controls=0)
grover3OR.add_operation("X", targets=2, anticontrols=1)
grover3OR.add_operation("X", targets=2, controls=0, anticontrols=1)


grover4 = qj.QGate(4, 0, "Uw4")
grover4.add_operation("X", targets=3, controls=0, anticontrols=1)
grover4.add_operation("X", targets=3, controls=2, anticontrols=0)


grover5 = qj.QGate(5, 0, "Uw5")
grover5.add_operation("X", targets=4, controls=0, anticontrols=1)
grover5.add_operation("X", targets=3, controls=2, anticontrols=4)
grover5.add_operation("X", targets=4, controls=0, anticontrols=1)
grover5.add_operation("X", targets=3)


grover6 = qj.QGate(6, 0, "Uw6")
grover6.add_operation("X", targets=5, controls={0, 2}, anticontrols=1)
grover6.add_operation("X", targets=4, anticontrols={3, 5})
grover6.add_operation("X", targets=5, controls={0, 2}, anticontrols=1)
grover6.add_operation("X", targets=4)


grover7 = qj.QGate(7, 0, "Uw7")
grover7.add_operation("X", targets=4, controls=0, anticontrols=1)
grover7.add_operation("X", targets=5, controls=0, anticontrols=2)
grover7.add_operation("X", targets=6, controls=1, anticontrols=2)
grover7.add_operation("X", targets=3, anticontrols={4, 5, 6})
grover7.add_operation("X", targets=6, controls=1, anticontrols=2)
grover7.add_operation("X", targets=5, controls=0, anticontrols=2)
grover7.add_operation("X", targets=4, controls=0, anticontrols=1)
grover7.add_operation("X", targets=3)


grover8 = qj.QGate(8, 0, "Uw8")
grover8.add_operation("X", targets=5, controls=0, anticontrols=2)
grover8.add_operation("X", targets=6, controls=3, anticontrols=2)
grover8.add_operation("X", targets=7, controls=0, anticontrols=1)
grover8.add_operation("X", targets=4, anticontrols={5, 6, 7})
grover8.add_operation("X", targets=7, controls=0, anticontrols=1)
grover8.add_operation("X", targets=6, controls=3, anticontrols=2)
grover8.add_operation("X", targets=5, controls=0, anticontrols=2)
grover8.add_operation("X", targets=4)


grover9 = qj.QGate(9, 0, "Uw9")
grover9.add_operation("X", targets=5, controls=2, anticontrols=0)
grover9.add_operation("X", targets=6, controls=2, anticontrols=3)
grover9.add_operation("X", targets=7, anticontrols={5, 6})
grover9.add_operation("X", targets=8, controls={0, 1})
grover9.add_operation("X", targets=4, anticontrols={7, 8})
grover9.add_operation("X", targets=8, controls={0, 1})
grover9.add_operation("X", targets=7, anticontrols={5, 6})
grover9.add_operation("X", targets=6, controls=2, anticontrols=3)
grover9.add_operation("X", targets=5, controls=2, anticontrols=0)
grover9.add_operation("X", targets=4)


grover10 = qj.QGate(10, 0, "Uw10")
grover10.add_operation("X", targets=6, controls={0, 2}, anticontrols=1)
grover10.add_operation("X", targets=7, controls=1)
grover10.add_operation("X", targets=7, controls=2)
grover10.add_operation("X", targets=8, anticontrols={1, 4})
grover10.add_operation("X", targets=9, controls={0, 3, 4, 7}, anticontrols={6, 8})
grover10.add_operation("X", targets=5, controls=2, anticontrols=9)
grover10.add_operation("X", targets=9, controls={0, 3, 4, 7}, anticontrols={6, 8})
grover10.add_operation("X", targets=8, anticontrols={1, 4})
grover10.add_operation("X", targets=7, controls=2)
grover10.add_operation("X", targets=7, controls=1)
grover10.add_operation("X", targets=6, controls={0, 2}, anticontrols=1)


def IAM(nq, last=False):
    gate = qj.QGate(nq, 0, f"IAM{nq}")
    for i in range(nq-1):
        gate.add_operation("H", targets=i)
    last_id = nq - 1
    gate.add_operation("X", targets=last_id)
    gate.add_operation("Z", targets=last_id, anticontrols={i for i in range(last_id)})
    if not last:
        gate.add_operation("X", targets=last_id)
    for i in range(nq-1):
        gate.add_operation("H", targets=i)

    return gate


oracles = [grover2, grover3OR, grover4, grover5, grover6, grover7, grover8, grover9, grover10]
no_ancilla = [2, 3, 4, 4, 5, 4, 5, 5, 6]


def bell_pair():
    sys = qj.QRegistry(2)
    sys = sys.apply_gate("H", targets=0)
    sys = sys.apply_gate("X", targets=1, controls=0)
    res, sys = sys.measure([0, 1])
    return (res, sys)


def deutsch_jozsa(nq, is_constant):
    start_time = t.time()
    
    circ = qj.QRegistry(nq + 1)
    circ = circ.apply_gate("X", targets=nq)
    for i in range(nq+1):
        circ = circ.apply_gate("H", targets=i)
    if is_constant: # f(x1, ..., xn) = True
        circ = circ.apply_gate("X", targets=nq)
    else:           # f(x1, ..., xn) = xn
        circ = circ.apply_gate("X", targets=nq, controls=nq-1)
    for i in range(nq):
        circ = circ.apply_gate("H", targets=i)
    circ2, results = circ.measure([i for i in range(nq)])
    
    stop_time = t.time()
    
    return (results[:nq], stop_time - start_time, 2 * doki.registry_mem(circ.reg, False))


def bernstein_vazirani(bit_str):
    nq = len(bit_str)
    start_time = t.time()
    
    circ = qj.QRegistry(nq + 1)
    circ = circ.apply_gate("X", targets=nq)
    for i in range(nq+1):
        circ = circ.apply_gate("H", targets=i)
    for i in range(nq):
        if bit_str[i]:
            circ = circ.apply_gate("X", targets=nq, controls=i)
    for i in range(nq):
        circ = circ.apply_gate("H", targets=i)
    circ2, results = circ.measure([i for i in range(nq)])
    
    stop_time = t.time()
    
    return (results[:nq], stop_time - start_time, 2 * doki.registry_mem(circ.reg, False))


def grover(nq):
    if nq < 2 or nq > 10:
        raise ValueError("That test has not been designed")
    UwU_function = oracles[nq-2]
    not_ancilla = no_ancilla[nq-2]
    start_time = t.time()
    
    circ = qj.QRegistry(nq)
    for i in range(not_ancilla):
        circ = circ.apply_gate("H", targets=i)
    circ = circ.apply_gate(UwU_function, targets=[i for i in range(nq)])
    circ = circ.apply_gate(IAM(not_ancilla, last=True), targets=[i for i in range(not_ancilla)])
    circ2, results = circ.measure([i for i in range(not_ancilla-1)])
    
    stop_time = t.time()
    
    return (results[:not_ancilla-1], stop_time - start_time, 2 * doki.registry_mem(circ.reg, False))


def main():
    print("Warmup...")
    for i in range(10):
        _, _ = bell_pair()
        _, _, _ = deutsch_jozsa(1, True)
        _, _, _ = bernstein_vazirani([False])
        _, _, _ = grover(2)
    print("Warm!")

    min_qb = 2
    max_qb = 10
    data = {"num_qubits": [], "DJ_times": [], "DJ_sizes": [], "BV_times": [], "BV_sizes": [], "grover_times": [], "grover_sizes": []}
    roll = None
    while True:
        for n in range(min_qb, max_qb + 1):
            print(f"{n} qubit tests")
            num_iters = 1
            for it in range(num_iters):
                roll = np.random.randint(2)
                is_constant = bool(roll)
                res, time, size = deutsch_jozsa(n-1, is_constant)
                if any(res) == is_constant:
                    print(is_constant)
                    print(res)
                    raise ValueError(f"Error when executing DJ with is_constant = {is_constant}, n = {n}")
                data["DJ_times"].append(time)
                data["DJ_sizes"].append(size)
                roll = np.random.randint(2, size=n-1)
                bit_str = [bool(bit) for bit in roll]
                res, time, size = bernstein_vazirani(bit_str)
                if bit_str != res:
                    print(bit_str)
                    print(res)
                    raise ValueError(f"Error when executing BV with n = {n}, bit_str = {bit_str}")
                data["BV_times"].append(time)
                data["BV_sizes"].append(size)
                res, time, size = grover(n)
                data["grover_times"].append(time)
                data["grover_sizes"].append(size)
                data["num_qubits"].append(n)
                '''
                if point * dec == it:
                    root_print(f"{(100 * (it + 1)) // num_iters}%")
                    point += 1
                '''
        pd.DataFrame(data).to_csv(f"fmtests/raw/dataQJ_{int(t.time())}_{min_qb}_{max_qb}.csv", sep=';', index=False)


if __name__ == "__main__":
    main()


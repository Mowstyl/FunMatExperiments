"""
Module usagetests_r_mpi - Provides functions to compare matrix performance.

Compares time and memory of dense, sparse and functional matrices
"""
import gates as g
import numpy as np
import scipy.sparse as sparse
import sys
import time as t

from funmatrix import Funmatrix
from mpi4py import MPI
from usagetests_h_mpi import deleteMem, deleteMemS, deleteMemF, \
                             accessTest, saveCSV


def getGate(code, nq):
    """Generate the functional matrix associated with the specified code.

    0 -> H
    1 -> Rx(Pi/nq)
    2 -> Ry(Pi/nq)
    3 -> Rz(Pi/nq)
    4 -> C-X
    5 -> C-Y
    6 -> C-Z
    """
    gate = None
    if code == 0:
        gate = g.H(1)
    elif code == 1:
        gate = g.Rx(np.pi/nq)
    elif code == 2:
        gate = g.Ry(np.pi/nq)
    elif code == 3:
        gate = g.Rz(np.pi/nq)
    elif code == 4:
        gate = g.CU(g.X())
    elif code == 5:
        gate = g.CU(g.Y())
    elif code == 6:
        gate = g.CU(g.Z())
    return gate


def gateMat(nq, gate_names, is_sparse=False, is_functional=False):
    """Generate a quantum gate and returns it with the creation time.

    Generates the matrix (dense, sparse or functional) associated with a
    quantum gate composed by the gates listed in gate_names, for the
    specified number of qubits.
    Return -> (matrix, time needed to calculate the kronecker product)
    """
    if not is_sparse and not is_functional:
        def gate_getter(gate_name):
            return getGate(gate_name, nq)[:]
        kron = np.kron
    elif not is_sparse and is_functional:
        def gate_getter(gate_name):
            return getGate(gate_name, nq)
        kron = Funmatrix.__pow__
    elif is_sparse and not is_functional:
        def gate_getter(gate_name):
            return sparse.dok_matrix(getGate(gate_name, nq)[:])

        def kron(gate1, gate2):
            return sparse.kron(gate1, gate2, "dok")
    else:
        raise ValueError("Functional sparse matrices do not exist!")
    gates = [gate_getter(gate) for gate in gate_names]
    t1 = t.time()
    mat = gates[0]
    for gate in gates[1:]:
        mat = kron(mat, gate)
    t2 = t.time()
    del gates
    t1 = int(t1 * 1000000000)  # To nanoseconds
    t2 = int(t2 * 1000000000)  # To nanoseconds

    return (mat, t2 - t1)


def runTests(minq, maxq, gateq):
    """Run all the tests and return the raw results."""
    res = []

    for nqubits in range(minq, maxq+1):
        dmat, tcd = gateMat(nqubits, gateq[nqubits])
        smat, tcs = gateMat(nqubits, gateq[nqubits], is_sparse=True)
        fmat, tcf = gateMat(nqubits, gateq[nqubits], is_functional=True)
        td = accessTest(dmat)
        sd = deleteMem(dmat)
        ts = accessTest(smat)
        ss = deleteMemS(smat)
        tf = accessTest(fmat)
        sf = deleteMemF(fmat)
        res.append((nqubits, sd, ss, sf, tcd, tcs, tcf, td, ts, tf))

    return res


def main():
    """Execute everything that has to be executed."""
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        def printRoot(thing):
            print(thing, flush=True)
    else:
        def printRoot(thing):
            pass

    argv = sys.argv[1:]
    if len(argv) < 3 or 5 < len(argv):
        printRoot("Usage: " + sys.argv[0] +
                  " <iterations> <minQubits> <maxQubits> <verbose> <seed>")
        return

    iterations = int(argv[0])
    if iterations <= 0:
        printRoot("<iterations> must be an integer greater than zero")
        return
    minq = int(argv[1])
    if minq <= 0:
        printRoot("<minq> must be an integer greater than zero")
        return
    maxq = int(argv[2])
    if maxq < minq:
        printRoot("<maxq> must be an integer greater than or equal to <minq>")
        return

    if len(argv) >= 4:
        verboseStr = argv[3].lower()
        if verboseStr != "true" and verboseStr != "false":
            printRoot("<verbose> must be true or false")
            return
    else:
        verboseStr = "false"
    verbose = verboseStr == "true"
    if verbose:
        def printVerbose(x):
            print(x, flush=True)
    else:
        def printVerbose(x):
            pass

    if len(argv) == 5:
        seed = int(argv[4])
        np.random.seed(seed=seed)
    else:
        seed = None

    baseIters = iterations // size
    remainder = iterations % size
    myIters = baseIters
    if rank < remainder:
        myIters += 1

    ntests = maxq - minq + 1
    resu = []
    st = t.time()

    gates1 = [0, 1, 2, 3]  # H, Rx, Ry, Rz
    gates2 = [4, 5, 6]  # C-X, C-Y, C-Z
    gates = {1: gates1, 2: gates2}
    gateq = None

    if rank == 0:
        gateq = {}
        for nqubits in range(minq, maxq+1):
            gateq[nqubits] = []
            tot = 0
            for i in range(nqubits):
                maxsize = nqubits - tot
                asizes = [size for size in gates.keys() if size <= maxsize]
                size = np.random.choice(asizes)
                gate = np.random.choice(gates[size])
                gateq[nqubits].append(gate)
                tot += size
                if tot == nqubits:
                    break
    printRoot(gateq)
    gateq = comm.bcast(gateq, root=0)

    for i in range(0, myIters + 1):
        if i == 0:
            printVerbose("Job[" + str(rank) + "]: Warming up...")
        res = runTests(minq, maxq, gateq)
        if i > 0:
            resu.append(res)
        nt = t.time()
        printVerbose("Job[" + str(rank) + "]: " +
                     str((i * 100) / myIters) + " % done...")
        if i < myIters:
            printVerbose("Job[" + str(rank) + "]: Elapsed time: " +
                         t.strftime('%H:%M:%S', t.gmtime(nt - st)))
            ait = (nt - st) / (i + 1)  # Average iteration time
            rem = (myIters - i) * ait
            printVerbose("Job[" + str(rank) + "]: Remaining time: " +
                         t.strftime('%H:%M:%S', t.gmtime(rem)))
    rarr = np.array(resu)
    allres = comm.gather(rarr, root=0)
    if rank == 0:
        resu = np.concatenate(allres).tolist()
        printVerbose("Done in " + t.strftime('%H:%M:%S', t.gmtime(nt - st)))

        fres = [[resu[0][i][0],
                 np.mean([resu[j][i][1] for j in range(iterations)]),
                 np.var([resu[j][i][1] for j in range(iterations)]),
                 np.std([resu[j][i][1] for j in range(iterations)]),
                 np.mean([resu[j][i][2] for j in range(iterations)]),
                 np.var([resu[j][i][2] for j in range(iterations)]),
                 np.std([resu[j][i][2] for j in range(iterations)]),
                 np.mean([resu[j][i][3] for j in range(iterations)]),
                 np.var([resu[j][i][3] for j in range(iterations)]),
                 np.std([resu[j][i][3] for j in range(iterations)]),
                 np.mean([resu[j][i][4] for j in range(iterations)]),
                 np.var([resu[j][i][4] for j in range(iterations)]),
                 np.std([resu[j][i][4] for j in range(iterations)]),
                 np.mean([resu[j][i][5] for j in range(iterations)]),
                 np.var([resu[j][i][5] for j in range(iterations)]),
                 np.std([resu[j][i][5] for j in range(iterations)]),
                 np.mean([resu[j][i][6] for j in range(iterations)]),
                 np.var([resu[j][i][6] for j in range(iterations)]),
                 np.std([resu[j][i][6] for j in range(iterations)]),
                 np.mean([resu[j][i][7] for j in range(iterations)]),
                 np.var([resu[j][i][7] for j in range(iterations)]),
                 np.std([resu[j][i][7] for j in range(iterations)]),
                 np.mean([resu[j][i][8] for j in range(iterations)]),
                 np.var([resu[j][i][8] for j in range(iterations)]),
                 np.std([resu[j][i][8] for j in range(iterations)]),
                 np.mean([resu[j][i][9] for j in range(iterations)]),
                 np.var([resu[j][i][9] for j in range(iterations)]),
                 np.std([resu[j][i][9] for j in range(iterations)])]
                for i in range(ntests)]
        infostr = str(minq) + "-" + str(maxq) + "_" + \
            t.strftime("%Y%m%d%H%M%S")
        saveCSV(resu, ("NumQubits", "SizeOfDMat", "SizeOfSMat", "SizeOfFMat",
                       "CreationTimeDMat", "CreationTimeSMat",
                       "CreationTimeFMat", "AccessTimeDMat", "AccessTimeSMat",
                       "AccessTimeFMat"),
                "dataRandom_" + infostr + ".csv", nested=True)
        saveCSV(fres, ("NumQubits",
                       "SizeOfDMat", "VarSizeOfDMat", "StdSizeOfDMat",
                       "SizeOfSMat", "VarSizeOfSMat", "StdSizeOfSMat",
                       "SizeOfFMat", "VarSizeOfFMat", "StdSizeOfFMat",
                       "AverageCreationTimeDMat", "VarCreationTimeDMat",
                       "StdCreationTimeDMat", "AverageCreationTimeSMat",
                       "VarCreationTimeSMat", "StdCreationTimeSMat",
                       "AverageCreationTimeFMat", "VarCreationTimeFMat",
                       "StdCreationTimeFMat", "AverageAccessTimeDMat",
                       "VarAccessTimeDMat", "StdAccessTimeDMat",
                       "AverageAccessTimeSMat", "VarAccessTimeSMat",
                       "StdAccessTimeSMat", "AverageAccessTimeFMat",
                       "VarAccessTimeFMat", "StdAccessTimeFMat"),
                "resultRandom_" + infostr + ".csv", verbose=verbose)


if __name__ == "__main__":
    main()

import gates as g
import csv
import numpy as np
import os
import psutil as psu
import sys
import time as t

from mpi4py import MPI


def getGate(code, nq):
    # [0, 1, 2, 3, 4, 5, 6] -> [H, Rx, Ry, Rz, C-X, C-Y, C-Z]
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


def gateMat(nq, gates, pinfo):
    t1 = 0.0
    t2 = 0.0
    preMat = pinfo.memory_info().rss
    gates = [getGate(gate, nq)[:] for gate in gates]
    t1 = t.time()
    mat = gates[0]
    for gate in gates[1:]:
        mat = np.kron(mat, gate)
    t2 = t.time()
    postMat = pinfo.memory_info().rss
    del gates
    t1 = int(t1 * 1000000000)  # To nanoseconds
    t2 = int(t2 * 1000000000)  # To nanoseconds

    return (mat, postMat - preMat, t2 - t1)


def gateFunc(nq, gates, pinfo):
    t1 = 0.0
    t2 = 0.0
    fmats = [getGate(gate, nq) for gate in gates]
    preMat = pinfo.memory_info().rss
    t1 = t.time()
    fmat = fmats[0]
    for gate in fmats[1:]:
        fmat = fmat ** gate
    t2 = t.time()
    postMat = pinfo.memory_info().rss
    t1 = int(t1 * 1000000000)  # To nanoseconds
    t2 = int(t2 * 1000000000)  # To nanoseconds

    return (fmat, postMat - preMat, t2 - t1)


def deleteMem(mat, pinfo):
    # postMat = pinfo.memory_info().rss
    # aux = mat[0,0]
    sizeof = sys.getsizeof(mat)
    del mat
    # preMat = pinfo.memory_info().rss

    return sizeof


def deleteMemF(fmat, pinfo):
    # postMat = pinfo.memory_info().rss
    # aux = mat[0,0]
    sizeof = fmat.getsizeof()
    del fmat
    # preMat = pinfo.memory_info().rss

    return sizeof


def accessTest(mat):
    aux = None
    t1 = t.time()
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            aux = mat[i, j]
    t2 = t.time()
    t1 = int(t1 * 1000000000)  # To nanoseconds
    t2 = int(t2 * 1000000000)  # To nanoseconds
    return int((t2 - t1)/(mat.shape[0]*mat.shape[1]))


def runTests(minq, maxq, gateq):
    proc = psu.Process(os.getpid())

    res = []

    for nqubits in range(minq, maxq+1):
        mat, rm, tcm = gateMat(nqubits, gateq[nqubits], proc)
        tm = accessTest(mat)
        dm = deleteMem(mat, proc)
        fmat, rf, tcf = gateFunc(nqubits, gateq[nqubits], proc)
        tf = accessTest(fmat)
        df = deleteMemF(fmat, proc)
        res.append((nqubits, rm, rf, dm, df, tcm, tcf, tm, tf))

    return res


def saveAll(results, filename):
    with open(filename, 'w', newline='') as result_file:
        csv_writer = csv.writer(result_file, delimiter=';')
        csv_writer.writerow(("NumQubits", "RAMUsageMat", "RAMUsageFMat",
                             "SizeOfMat", "SizeOfFMat", "MeanCreationTimeMat",
                             "MeanCreationTimeFMat", "MeanAccessTimeMat",
                             "MeanAccessTimeFMat"))
        for result in results:
            for tup in result:
                tup[0] = int(tup[0])
                csv_writer.writerow(tup)


def saveResults(result, filename, verbose):
    with open(filename, 'w', newline='') as result_file:
        csv_writer = csv.writer(result_file, delimiter=';')
        csv_writer.writerow(("NumQubits",
                             "SizeOfMat",
                             "VarSizeOfMat",
                             "StdSizeOfMat",
                             "SizeOfFMat",
                             "VarSizeOfFMat",
                             "StdSizeOfFMat",
                             "MeanCreationTimeMat",
                             "VarMeanCreationTimeMat",
                             "StdMeanCreationTimeMat",
                             "MeanCreationTimeFMat",
                             "VarMeanCreationTimeFMat",
                             "StdMeanCreationTimeFMat",
                             "MeanAccessTimeMat",
                             "VarMeanAccessTimeMat",
                             "StdMeanAccessTimeMat",
                             "MeanAccessTimeFMat",
                             "VarMeanAccessTimeFMat",
                             "StdMeanAccessTimeFMat"))
        for tup in result:
            tup[0] = int(tup[0])
            if (verbose):
                print(tup, flush=True)
            csv_writer.writerow(tup)


def main():
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    def printRoot(thing):
        if rank == 0:
            print(thing, flush=True)

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
            print(x)
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
                 np.std([resu[j][i][8] for j in range(iterations)])]
                for i in range(ntests)]
        infostr = str(minq) + "-" + str(maxq) + "_" + \
            t.strftime("%Y%m%d%H%M%S")
        saveAll(resu, "dataRandom_" + infostr + ".csv")
        saveResults(fres, "resultRandom_" + infostr + ".csv", verbose)


if __name__ == "__main__":
    main()

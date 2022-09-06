"""
Module usagetests_h_mpi - Provides functions to compare matrix performance.

Compares time and memory of dense, sparse and functional matrices
"""
import csv
import gates as g
import numpy as np
import sys
import scipy.sparse as sparse
import time as t

from mpi4py import MPI


def hadamardMat(nq):
    """Create dense and sparse matrix of H gate, returns metrics.

    Creates both dense and sparse matrix that represents the H quantum
    gate for the given number of qubits. Also measures access time and
    memory requirements for the specified functional matrix.

    Return -> (dense matrix,
               sparse matrix,
               functional matrix,
               average access time of functional matrix)
    """
    fmat = g.H(nq)
    t1 = 0.0
    t2 = 0.0
    t1 = t.time()
    mat = g.H(nq)[:]
    t2 = t.time()
    smat = sparse.dok_matrix(mat)
    t1 = int(t1 * 1000000000)  # To nanoseconds
    t2 = int(t2 * 1000000000)  # To nanoseconds

    return (mat, smat, fmat, (t2 - t1))


def deleteMem(mat):
    """Delete dense matrix and returns its size."""
    sizeof = sys.getsizeof(mat)
    del mat
    return sizeof


def deleteMemS(smat):
    """Delete sparse matrix and returns its size."""
    sizeof = sys.getsizeof(smat)
    del smat
    return sizeof


def deleteMemF(fmat):
    """Delete functional matrix and returns its size."""
    sizeof = fmat.getsizeof()
    del fmat
    return sizeof


def accessTest(mat):
    """Access all the elements of a matrix. Return the time needed (ns)."""
    t1 = t.time()
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            _ = mat[i, j]
    t2 = t.time()
    t1 = int(t1 * 1000000000)  # To nanoseconds
    t2 = int(t2 * 1000000000)  # To nanoseconds
    return t2 - t1


def runTests(minq, maxq):
    """Run all the tests and return the raw results."""
    res = []

    for nqubits in range(minq, maxq+1):
        mat, smat, fmat, tf = hadamardMat(nqubits)
        tm = accessTest(mat)
        ts = accessTest(smat)
        dm = deleteMem(mat)
        ds = deleteMem(smat)
        df = deleteMemF(fmat)
        res.append((nqubits, dm, ds, df, tm, ts, tf))

    return res


def saveCSV(results, column_names, filename, nested=False, verbose=False):
    """Save the results to the specified file in csv format."""
    def writerow(data):
        data[0] = int(data[0])
        if (verbose):
            print(data, flush=True)
        csv_writer.writerow(data)

    with open(filename, 'w', newline='') as result_file:
        csv_writer = csv.writer(result_file, delimiter=';')
        csv_writer.writerow(column_names)
        for result in results:
            if nested:
                for tup in result:
                    writerow(tup)
            else:
                writerow(result)


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
    if len(argv) != 4:
        printRoot("Usage: " + sys.argv[0] +
                  " <iterations> <minQubits> <maxQubits> <verbose>")
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
    verboseStr = argv[3].lower()
    if verboseStr != "true" and verboseStr != "false":
        printRoot("<verbose> must be true or false")
        return
    verbose = verboseStr == "true"

    if verbose:
        def printVerbose(x):
            print(x, flush=True)
    else:
        def printVerbose(x):
            pass

    baseIters = iterations // size
    remainder = iterations % size
    myIters = baseIters
    if rank < remainder:
        myIters += 1

    ntests = maxq - minq + 1
    resu = []
    st = t.time()
    for i in range(0, myIters + 1):
        if i == 0:
            printVerbose("Job[" + str(rank) + "]: Warming up...")
        res = runTests(minq, maxq)
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
                 np.std([resu[j][i][6] for j in range(iterations)])]
                for i in range(ntests)]
        infostr = str(minq) + "-" + str(maxq) + "_" + \
            t.strftime("%Y%m%d%H%M%S")
        saveCSV(resu, ("NumQubits", "SizeOfDMat", "SizeOfSMat", "SizeOfFMat",
                       "AccessTimeDMat", "AccessTimeSMat", "AccessTimeFMat"),
                "dataHadamard_" + infostr + ".csv", nested=True)
        saveCSV(fres, ("NumQubits",
                       "SizeOfDMat", "VarSizeOfDMat", "StdSizeOfDMat",
                       "SizeOfSMat", "VarSizeOfSMat", "StdSizeOfSMat",
                       "SizeOfFMat", "VarSizeOfFMat", "StdSizeOfFMat",
                       "AverageAccessTimeDMat", "VarAccessTimeDMat",
                       "StdAccessTimeDMat", "AverageAccessTimeSMat",
                       "VarAccessTimeSMat", "StdAccessTimeSMat",
                       "AverageAccessTimeFMat", "VarAccessTimeFMat",
                       "StdMeanAccessTimeFMat"),
                "resultHadamard_" + infostr + ".csv", verbose=verbose)


if __name__ == "__main__":
    main()

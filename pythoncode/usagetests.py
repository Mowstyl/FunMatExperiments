import csv
import gates as g
import numpy as np
import os
import psutil as psu
import sys
import time as t


def hadamardMat(nq, pinfo):
    fmat = g.H(nq)
    t1 = 0.0
    t2 = 0.0
    preMat = pinfo.memory_info().rss
    t1 = t.time()
    mat = g.H(nq)[:]
    t2 = t.time()
    postMat = pinfo.memory_info().rss
    del fmat
    t1 = int(t1 * 10000000) # To nanoseconds
    t2 = int(t2 * 10000000) # To nanoseconds

    return (mat, postMat - preMat, (t2 - t1)/(mat.shape[0]*mat.shape[1]))

def hadamardFunc(nq, pinfo):
    preMat = pinfo.memory_info().rss
    fmat = g.H(nq)
    postMat = pinfo.memory_info().rss

    return (fmat, postMat - preMat)

def deleteMem(mat, pinfo):
    #postMat = pinfo.memory_info().rss
    #aux = mat[0,0]
    sizeof = sys.getsizeof(mat)
    del mat
    #preMat = pinfo.memory_info().rss

    return sizeof

def deleteMemF(fmat, pinfo):
    #postMat = pinfo.memory_info().rss
    #aux = mat[0,0]
    sizeof = fmat.getsizeof()
    del fmat
    #preMat = pinfo.memory_info().rss

    return sizeof

def accessTest(mat):
    aux = None
    t1 = t.time()
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            aux = mat[i, j]
    t2 = t.time()
    t1 = int(t1 * 10000000) # To nanoseconds
    t2 = int(t2 * 10000000) # To nanoseconds
    return (t2 - t1)/(mat.shape[0]*mat.shape[1])

def runTests(minq, maxq):
    proc = psu.Process(os.getpid())

    res = []

    for nqubits in range(minq, maxq+1):
        mat, rm, tf = hadamardMat(nqubits, proc)
        tm = accessTest(mat)
        dm = deleteMem(mat, proc)
        fmat, rf = hadamardFunc(nqubits, proc)
        df = deleteMemF(fmat, proc)
        res.append((nqubits, rm, rf, dm, df, tm, tf))

    return res

def saveAll(results, filename):
    with open(filename, 'w', newline='') as result_file:
        csv_writer = csv.writer(result_file, delimiter=';')
        csv_writer.writerow(("NumQubits", "RAMUsageMat", "RAMUsageFMat", "SizeOfMat", "SizeOfFMat", "MeanAccesTimeMat", "MeanAccesTimeFMat"))
        for result in results:
            for tup in result:
                csv_writer.writerow(tup)

def saveResults(result, filename):
    with open(filename, 'w', newline='') as result_file:
        csv_writer = csv.writer(result_file, delimiter=';')
        csv_writer.writerow(("NumQubits",
                             "SizeOfMat",
                             "VarSizeOfMat",
                             "StdSizeOfMat",
                             "SizeOfFMat",
                             "VarSizeOfFMat",
                             "StdSizeOfFMat",
                             "MeanAccesTimeMat",
                             "VarMeanAccesTimeMat",
                             "StdMeanAccesTimeMat",
                             "MeanAccesTimeFMat",
                             "VarMeanAccesTimeFMat",
                             "StdMeanAccesTimeFMat"))
        for tup in result:
            print(tup)
            csv_writer.writerow(tup)

def main():
    iterations = 100
    minq = 1
    maxq = 5

    ntests = maxq - minq + 1
    resu = []
    print("Starting!")
    st = t.time()
    for i in range(0, iterations + 1):
        if i == 0:
            print("Warming up...")
        res = runTests(minq, maxq)
        if i > 0:
            resu.append(res)
        nt = t.time()
        print(str((i * 100) / iterations) + " % done...")
        if i < iterations:
            print("Elapsed time: " + t.strftime('%H:%M:%S', t.gmtime(nt - st)))
            #print("Iterations: " + str(i + 1))
            #print("Remaining: " + str(iterations - i))
            ait = (nt - st) / (i + 1) # Average iteration time
            #print("Average: " + str(ait))
            rem = (iterations - i) * ait
            print("Remaining time: " + t.strftime('%H:%M:%S', t.gmtime(rem)))
    print("Done in " + t.strftime('%H:%M:%S', t.gmtime(nt - st)))

    fres = [(resu[0][i][0],
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
             np.std([resu[j][i][6] for j in range(iterations)])) for i in range(ntests)]
    saveAll(resu, "data.csv")
    saveResults(fres, "result.csv")

if __name__ == "__main__":
    main()

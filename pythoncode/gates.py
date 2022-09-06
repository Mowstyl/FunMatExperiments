import ctypes as ct
import funmatrix as fm
import numpy as np
import platform as plat

# DLL Load
if plat.system() == "Windows":
    extension = ".dll"
else:
    extension = ".so"

__qgates__ = ct.CDLL("./libqgates" + extension)

c_double_p = ct.POINTER(ct.c_double)

__cIdentity__ = __qgates__.Identity
__cIdentity__.argtypes = [ct.c_int]
__cIdentity__.restype = ct.c_void_p


def I(n):
    mat = __cIdentity__(ct.c_int(n))

    return fm.Funmatrix(mat, "I(" + str(n) + ")")


__cWalsh__ = __qgates__.Walsh
__cWalsh__.argtypes = [ct.c_int]
__cWalsh__.restype = ct.c_void_p


def Walsh(n):
    mat = __cWalsh__(ct.c_int(n))

    return fm.Funmatrix(mat, "Walsh(" + str(n) + ")")


__cHadamard__ = __qgates__.H
__cHadamard__.argtypes = [ct.c_int]
__cHadamard__.restype = ct.c_void_p


def H(n):
    mat = __cHadamard__(ct.c_int(n))

    return fm.Funmatrix(mat, "H(" + str(n) + ")")


__cQFT__ = __qgates__.QFT
__cQFT__.argtypes = [ct.c_int]
__cQFT__.restype = ct.c_void_p


def QFT(n):
    mat = __cQFT__(ct.c_int(n))

    return fm.Funmatrix(mat, "QFT(" + str(n) + ")")


__cX__ = __qgates__.X
__cX__.argtypes = []
__cX__.restype = ct.c_void_p


def X():
    mat = __cX__()

    return fm.Funmatrix(mat, "X")


__cY__ = __qgates__.Y
__cY__.argtypes = []
__cY__.restype = ct.c_void_p


def Y():
    mat = __cY__()

    return fm.Funmatrix(mat, "Y")


__cZ__ = __qgates__.Z
__cZ__.argtypes = []
__cZ__.restype = ct.c_void_p


def Z():
    mat = __cZ__()

    return fm.Funmatrix(mat, "Z")


__cRX__ = __qgates__.RX
__cRX__.argtypes = [ct.c_double]
__cRX__.restype = ct.c_void_p


def Rx(angle):
    mat = __cRX__(ct.c_double(angle))

    return fm.Funmatrix(mat, "Rx(" + str(angle) + ")")


__cRY__ = __qgates__.RY
__cRY__.argtypes = [ct.c_double]
__cRY__.restype = ct.c_void_p


def Ry(angle):
    mat = __cRY__(ct.c_double(angle))

    return fm.Funmatrix(mat, "Ry(" + str(angle) + ")")


__cRZ__ = __qgates__.RZ
__cRZ__.argtypes = [ct.c_double]
__cRZ__.restype = ct.c_void_p


def Rz(angle):
    mat = __cRZ__(ct.c_double(angle))

    return fm.Funmatrix(mat, "Rz(" + str(angle) + ")")


__cSqrtX__ = __qgates__.SqrtX
__cSqrtX__.argtypes = []
__cSqrtX__.restype = ct.c_void_p


def SqrtX():
    mat = __cSqrtX__()

    return fm.Funmatrix(mat, "SqrtX")


__cCNOT__ = __qgates__.CNOT
__cCNOT__.argtypes = [ct.c_int, ct.c_int]
__cCNOT__.restype = ct.c_void_p

def CNOT():
    mat = __cCNOT__(ct.c_int(1), ct.c_int(0))

    return fm.Funmatrix(mat, "CNOT")

__cSWAP__ = __qgates__.SWAP
__cSWAP__.argtypes = []
__cSWAP__.restype = ct.c_void_p

def SWAP():
    mat = __cSWAP__()

    return fm.Funmatrix(mat, "SWAP")

__cISWAP__ = __qgates__.ISWAP
__cISWAP__.argtypes = []
__cISWAP__.restype = ct.c_void_p

def iSWAP():
    mat = __cISWAP__()

    return fm.Funmatrix(mat, "iSWAP")

__cU__ = __qgates__.U
__cU__.argtypes = [ct.c_double, ct.c_double, ct.c_double]
__cU__.restype = ct.c_void_p

def U(theta, phi, lamb):
    mat = __cU__(ct.c_double(theta), ct.c_double(phi), ct.c_double(lamb))

    return fm.Funmatrix(mat, "U(" + str(theta) + "," + str(phi) + "," + str(lamb) + ")")

__cU2__ = __qgates__.U2
__cU2__.argtypes = [ct.c_double, ct.c_double]
__cU2__.restype = ct.c_void_p

def U2(theta, phi):
    mat = __cU2__(ct.c_double(theta), ct.c_double(phi))

    return fm.Funmatrix(mat, "U2(" + str(theta) + "," + str(phi) + ")")

__cU1__ = __qgates__.U1
__cU1__.argtypes = [ct.c_double]
__cU1__.restype = ct.c_void_p

def U1(theta):
    mat = __cU1__(ct.c_double(theta))

    return fm.Funmatrix(mat, "U1(" + str(theta) + ")")

__cIAA__ = __qgates__.IAA
__cIAA__.argtypes = [ct.c_int]
__cIAA__.restype = ct.c_void_p

def IAA(n):
    mat = __cIAA__(ct.c_int(n))

    return fm.Funmatrix(mat, "IAA(" + str(n) + ")")

__cPyCustomGate__ = __qgates__.PyCustomGate
__cPyCustomGate__.argtypes = [c_double_p, c_double_p, ct.c_uint, ct.c_uint]
__cPyCustomGate__.restype = ct.c_void_p

def CustomGate(matrix, name="UNNAMED"):
    nm = np.array(matrix)
    r, c = nm.shape
    size = r * c
    flatreal = nm.reshape(size)
    flatimag = None
    if flatreal.dtype.type is np.dtype(np.complex).type:
        flatimag = np.array([num.imag for num in flatreal])
        flatreal = np.array([num.real for num in flatreal])
    else:
        flatimag = np.array([0.0 for num in flatreal])

    double_array = ct.c_double * size

    cfr = double_array(*flatreal)
    cfi = double_array(*flatimag)

    mat = __cPyCustomGate__(cfr, cfi, ct.c_uint(r), ct.c_uint(size))

    return fm.Funmatrix(mat, name)

__cCU__ = __qgates__.CU
__cCU__.argtypes = [ct.c_void_p]
__cCU__.restype = ct.c_void_p

def CU(gate):
    mat = __cCU__(ct.c_void_p(gate.m))

    return fm.Funmatrix(mat, "C-" + gate.name)

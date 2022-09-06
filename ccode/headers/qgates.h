#ifndef __QGATES_H
#define __QGATES_H

#include <complex.h>
#include "precision.h"
#include "funmatrix.h"


__attribute__ ((const))
unsigned int      intLog2           (unsigned int      N);

__attribute__ ((pure))
int               getGateQubits     (FunctionalMatrix *U);

__attribute__ ((const))
REAL _Complex     _IdentityFunction (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *unused3 __attribute__((unused)));

FunctionalMatrix* Identity          (int               n);

__attribute__ ((const))
REAL _Complex     _WalshFunction    (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          size,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     void             *unused2 __attribute__((unused)));

FunctionalMatrix* Walsh             (int               n);

FunctionalMatrix* H                 (int               n);

__attribute__ ((const))
REAL _Complex     _QFTFunction      (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          size,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     void             *unused2 __attribute__((unused)));

FunctionalMatrix* QFT               (int               n);

__attribute__ ((const))
REAL _Complex     _XFunction        (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *unused3 __attribute__((unused)));

FunctionalMatrix* X                 ();

__attribute__ ((const))
REAL _Complex     _YFunction        (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *unused3 __attribute__((unused)));

FunctionalMatrix* Y                 ();

__attribute__ ((const))
REAL _Complex     _ZFunction        (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *unused3 __attribute__((unused)));

FunctionalMatrix* Z                 ();

__attribute__ ((pure))
REAL _Complex     _RXFunction       (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *angle);

FunctionalMatrix* RX                (REAL              phi);

__attribute__ ((pure))
REAL _Complex     _RYFunction       (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *angle);

FunctionalMatrix* RY                (REAL              phi);

__attribute__ ((pure))
REAL _Complex     _RZFunction       (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *angle);

FunctionalMatrix* RZ                (REAL              phi);


__attribute__ ((const))
REAL _Complex     _SXFunction       (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *unused3 __attribute__((unused)));

FunctionalMatrix* SqrtX             ();

__attribute__ ((pure))
REAL _Complex     _CNOTFunction     (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *dinv);

FunctionalMatrix* CNOT              (int               distance,
                                     int               inversion);

__attribute__ ((const))
REAL _Complex     _SWAPFunction     (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *unused3 __attribute__((unused)));

FunctionalMatrix* SWAP              ();

__attribute__ ((const))
REAL _Complex     _ISWAPFunction    (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *unused3 __attribute__((unused)));

FunctionalMatrix* ISWAP             ();

__attribute__ ((pure))
REAL _Complex     _UFunction        (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *angles);

FunctionalMatrix* U                 (REAL              theta,
                                     REAL              phi,
                                     REAL              lambda);

__attribute__ ((pure))
REAL _Complex     _U2Function       (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *angles);

FunctionalMatrix* U2                (REAL              phi,
                                     REAL              lambda);

__attribute__ ((pure))
REAL _Complex     _U1Function       (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *rawangle);

FunctionalMatrix* U1                (REAL              angle);

__attribute__ ((const))
REAL _Complex     _InvFunction      (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          nrows,
                                     UINTEGER          unused __attribute__((unused)),
                                     void             *unused2 __attribute__((unused)));

FunctionalMatrix* IAA               (int               nqubits); // Inversion about the average

__attribute__ ((pure))
REAL _Complex     _CustomGate       (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          nrows,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     void             *matrix_2d);

FunctionalMatrix* CustomGate        (REAL _Complex    *matrix_2d,
                                     UINTEGER          nrows);

FunctionalMatrix* PyCustomGate      (REAL             *re_2d,
                                     REAL             *im_2d,
                                     UINTEGER          nrows,
                                     UINTEGER          size);

__attribute__ ((pure))
REAL _Complex     _CUFunction       (UINTEGER          i,
                                     UINTEGER          j,
                                     UINTEGER          unused1 __attribute__((unused)),
                                     UINTEGER          unused2 __attribute__((unused)),
                                     void             *RawU);

FunctionalMatrix* CU                (FunctionalMatrix *U); // Gives the Controlled-U gate

#endif

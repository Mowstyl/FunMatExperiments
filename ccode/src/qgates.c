#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "precision.h"
#include "qgates.h"
#include "funmatrix.h"

REAL _Complex
_IdentityFunction (UINTEGER  i,
                   UINTEGER  j,
                   UINTEGER  unused1 __attribute__((unused)),
                   UINTEGER  unused2 __attribute__((unused)),
                   void     *unused3 __attribute__((unused)))
{
  short number;

  number = 0;
  if (i == j)
    number = 1;

  return number;
}

FunctionalMatrix*
Identity (int n)
{
  FunctionalMatrix* pFM;
  UINTEGER size;

  size = 2 << (n - 1); // 2^n
  pFM = new_FunctionalMatrix(size, size, &_IdentityFunction, NULL);

  return pFM;
}

REAL _Complex
_WalshFunction (UINTEGER  i,
                UINTEGER  j,
                UINTEGER  size,
                UINTEGER  unused1 __attribute__((unused)),
                void     *unused2 __attribute__((unused)))
{
  short number;
  UINTEGER mid;

  mid = size / 2;
  number = 1;
  if (size == 2)
    {
      if (i == 1 && j == 1)
        number = -1;
    }
  else
    {
      if (i >= mid && j >= mid)
        {
          number = -_WalshFunction(i - mid, j - mid, mid, 0, NULL);
        }
      else
        {
          if (i >= mid)
            i = i - mid;

          if (j >= mid)
            j = j - mid;

          number = _WalshFunction(i, j, mid, 0, NULL);
        }
    }

  return number;
}

FunctionalMatrix*
Walsh (int n)
{
  FunctionalMatrix* pFM;
  UINTEGER size;

  size = 2 << (n - 1); // 2^n
  pFM = new_FunctionalMatrix(size, size, &_WalshFunction, NULL);

  return pFM;
}

FunctionalMatrix*
H (int n)
{
  FunctionalMatrix *pFM, *aux;
  UINTEGER size;

  size = 2 << (n - 1); // 2^n
  aux = new_FunctionalMatrix(size, size, &_WalshFunction, NULL);
  pFM = mprod(1/sqrt(size), aux);
  free(aux);
  aux = NULL;

  return pFM;
}

REAL _Complex
_QFTFunction (UINTEGER  i,
              UINTEGER  j,
              UINTEGER  size,
              UINTEGER  unused1 __attribute__((unused)),
              void     *unused2 __attribute__((unused)))
{
  REAL c = i * j * Q_TAU / size;
  return COS(c) + SIN(c) * I;
}

FunctionalMatrix*
QFT (int n)
{
  FunctionalMatrix *pFM, *aux;
  UINTEGER size;

  size = 2 << (n - 1); // 2^n

  aux = new_FunctionalMatrix(size, size, &_QFTFunction, NULL);
  pFM = mprod(1/sqrt(size), aux);
  free(aux);
  aux = NULL;

  return pFM;
}

REAL _Complex
_XFunction (UINTEGER  i,
            UINTEGER  j,
            UINTEGER  unused1 __attribute__((unused)),
            UINTEGER  unused2 __attribute__((unused)),
            void     *unused3 __attribute__((unused)))
{
  short val;

  val = i - j;
  if (val < 0)
    val = -val;

  return val;
}

FunctionalMatrix*
X ()
{
  FunctionalMatrix *pFM;

  pFM = new_FunctionalMatrix(2, 2, &_XFunction, NULL);

  return pFM;
}

REAL _Complex
_YFunction (UINTEGER  i,
            UINTEGER  j,
            UINTEGER  unused1 __attribute__((unused)),
            UINTEGER  unused2 __attribute__((unused)),
            void     *unused3 __attribute__((unused)))
{
  return (short)(i - j) * I;
}

FunctionalMatrix*
Y ()
{
  FunctionalMatrix *pFM;

  pFM = new_FunctionalMatrix(2, 2, &_YFunction, NULL);

  return pFM;
}

REAL _Complex
_ZFunction (UINTEGER  i,
            UINTEGER  j,
            UINTEGER  unused1 __attribute__((unused)),
            UINTEGER  unused2 __attribute__((unused)),
            void     *unused3 __attribute__((unused)))
{
  return 1 - (short)(i + j);
}

FunctionalMatrix*
Z ()
{
  FunctionalMatrix *pFM;

  pFM = new_FunctionalMatrix(2, 2, &_ZFunction, NULL);

  return pFM;
}

REAL _Complex
_RXFunction (UINTEGER  i,
             UINTEGER  j,
             UINTEGER  unused1 __attribute__((unused)),
             UINTEGER  unused2 __attribute__((unused)),
             void     *angle)
{ // here phi must be a pointer to a real number, the angle
  REAL phi = *(REAL*) angle;

  return (i == j) * COS(phi / 2) - (i != j) * SIN(phi / 2) * I;
}

FunctionalMatrix*
RX (REAL phi)
{
  FunctionalMatrix *pFM;
  REAL *value = (REAL*) malloc (sizeof(REAL));
  *value = phi;

  pFM = new_FunctionalMatrix(2, 2, &_RXFunction, value);

  return pFM;
}

REAL _Complex
_RYFunction (UINTEGER  i,
             UINTEGER  j,
             UINTEGER  unused1 __attribute__((unused)),
             UINTEGER  unused2 __attribute__((unused)),
             void     *angle)
{ // here phi must be a pointer to a real number, the angle
  REAL phi = *(REAL*) angle;

  return (i == j) * COS(phi/2) + (i != j) * SIN(phi/2) * (short) (i - j) + 0 * I;
}

FunctionalMatrix*
RY (REAL phi)
{
  FunctionalMatrix *pFM;
  REAL *value = (REAL*) malloc (sizeof(REAL));
  *value = phi;

  pFM = new_FunctionalMatrix(2, 2, &_RYFunction, value);

  return pFM;
}

REAL _Complex
_RZFunction (UINTEGER  i,
             UINTEGER  j,
             UINTEGER  unused1 __attribute__((unused)),
             UINTEGER  unused2 __attribute__((unused)),
             void     *angle)
{ // here angle must be a pointer to a real number, the angle¡
  REAL phi2 = *(REAL*) angle;

  phi2 = phi2 / 2;

  return (i == j) * (COS (i * phi2 - !i * phi2) + SIN (i * phi2 - !i * phi2) * I);
}

FunctionalMatrix*
RZ (REAL phi)
{
  FunctionalMatrix *pFM;
  REAL *value = (REAL*) malloc (sizeof(REAL));
  *value = phi;

  pFM = new_FunctionalMatrix(2, 2, &_RZFunction, value);

  return pFM;
}

REAL _Complex
_SXFunction (UINTEGER  i,
             UINTEGER  j,
             UINTEGER  unused1 __attribute__((unused)),
             UINTEGER  unused2 __attribute__((unused)),
             void     *unused3 __attribute__((unused)))
{
  REAL _Complex val;

  if (i == j)
    val = 0.5 + 0.5 * I;
  else
    val = 0.5 - 0.5 * I;

  return val;
}

FunctionalMatrix*
SqrtX ()
{
  FunctionalMatrix *pFM;

  pFM = new_FunctionalMatrix(2, 2, &_SXFunction, NULL);

  return pFM;
}

REAL _Complex
_CNOTFunction (UINTEGER  i,
               UINTEGER  j,
               UINTEGER  unused1 __attribute__((unused)),
               UINTEGER  unused2 __attribute__((unused)),
               void     *dinv)
{
  int *distinv = (int*) dinv,
       bitsi[2],
       bitsj[2],
       aux = 1,
       k;

  bitsi[!distinv[1]] = i & 1;
  bitsj[!distinv[1]] = j & 1;
  i = i >> 1;
  j = j >> 1;
  for (k = 1; k < distinv[0]; k++)
    {
      aux = aux && ((i & 1) == (j & 1));
      i = i >> 1;
      j = j >> 1;
    }
  bitsi[distinv[1]] = i;
  bitsj[distinv[1]] = j;

  return aux && bitsi[0] == bitsj[0] && ((!bitsi[0] && bitsi[1] == bitsj[1]) || (bitsi[0] && bitsi[1] != bitsj[1]));
}

FunctionalMatrix*
CNOT (int distance,
      int inversion)
{
  FunctionalMatrix *pFM;
  int *distinv = (int*) malloc (2 * sizeof(int));
  UINTEGER size = 2 << distance;
  distinv[0] = distance;
  distinv[1] = inversion;

  pFM = new_FunctionalMatrix(size, size, &_CNOTFunction, distinv);

  return pFM;
}

REAL _Complex
_SWAPFunction (UINTEGER  i,
               UINTEGER  j,
               UINTEGER  unused1 __attribute__((unused)),
               UINTEGER  unused2 __attribute__((unused)),
               void     *unused3 __attribute__((unused)))
{
  return ((i == j && (i == 0 || j == 3)) || ((i == 1 && j == 2) || (i == 2 && j == 1))) + 0 * I;
}

FunctionalMatrix*
SWAP ()
{
  FunctionalMatrix *pFM;

  pFM = new_FunctionalMatrix(4, 4, &_SWAPFunction, NULL);

  return pFM;
}

REAL _Complex
_ISWAPFunction (UINTEGER  i,
                UINTEGER  j,
                UINTEGER  unused1 __attribute__((unused)),
                UINTEGER  unused2 __attribute__((unused)),
                void     *unused3 __attribute__((unused)))
{
  return (i == j && (i == 0 || j == 3)) + ((i == 1 && j == 2) || (i == 2 && j == 1)) * I;
}

FunctionalMatrix*
ISWAP ()
{
  FunctionalMatrix *pFM;

  pFM = new_FunctionalMatrix(4, 4, &_ISWAPFunction, NULL);

  return pFM;
}

REAL _Complex
_UFunction (UINTEGER  i,
            UINTEGER  j,
            UINTEGER  unused1 __attribute__((unused)),
            UINTEGER  unused2 __attribute__((unused)),
            void     *angles)
{ // here angle must be a pointer to 3 real numbers, the angles
  REAL *rawangles = (REAL*) angles;
  REAL _Complex val;

  if (i == 0 && j == 1)
		val = (-COS(rawangles[2]) - SIN(rawangles[2]) * I) * SIN(rawangles[0]/2);
	else if (i == 1 && j == 0)
		val = (COS(rawangles[1]) + SIN(rawangles[1]) * I) * SIN(rawangles[0]/2);
	else if (i == 1 && j == 1)
		val = (COS(rawangles[1]+rawangles[2]) + SIN(rawangles[1]+rawangles[2]) * I) * COS(rawangles[0]/2);
	else
		val = COS(rawangles[0]/2);

  return val;
}

FunctionalMatrix*
U (REAL theta,
   REAL phi,
   REAL lambda)
{
  FunctionalMatrix *pFM;
  REAL *value = (REAL*) malloc (3 * sizeof(REAL));
  value[0] = theta;
  value[1] = phi;
  value[2] = lambda;

  pFM = new_FunctionalMatrix (2, 2, &_UFunction, value);

  return pFM;
}

REAL _Complex
_U2Function (UINTEGER  i,
             UINTEGER  j,
             UINTEGER  unused1 __attribute__((unused)),
             UINTEGER  unused2 __attribute__((unused)),
             void     *angles)
{ // here angle must be a pointer to 2 real numbers, the angles
  REAL *rawangles = (REAL*) angles;
  REAL _Complex val;

  if (i == 0 && j == 1)
		val = (-COS (rawangles[1]) - SIN (rawangles[1]) * I) / Q_SQRT2;
	else if (i == 1 && j == 0)
		val = (COS (rawangles[0]) + SIN (rawangles[0]) * I) / Q_SQRT2;
	else if (i == 1 && j == 1)
		val = (COS (rawangles[0] + rawangles[1]) + SIN (rawangles[0] + rawangles[1]) * I) / Q_SQRT2;
	else
		val = Q_SQRT1_2;

  return val;
}

FunctionalMatrix*
U2 (REAL phi,
    REAL lambda)
{
  FunctionalMatrix *pFM;
  REAL *value = (REAL*) malloc (2 * sizeof(REAL));
  value[0] = phi;
  value[1] = lambda;

  pFM = new_FunctionalMatrix (2, 2, &_UFunction, value);

  return pFM;
}

REAL _Complex
_U1Function (UINTEGER  i,
             UINTEGER  j,
             UINTEGER  unused1 __attribute__((unused)),
             UINTEGER  unused2 __attribute__((unused)),
             void     *rawangle)
{ // here angle must be a pointer to 1 real number, the angle
  REAL angle = *((REAL*) rawangle);
  REAL _Complex val;

  if (i == 1 && j == 1)
    val = COS (angle) + SIN (angle) * I;
  else
    val = !abs ((INTEGER) (i-j));

  return val;
}

FunctionalMatrix*
U1 (REAL angle)
{
  FunctionalMatrix *pFM;
  REAL *value = (REAL*) malloc (sizeof(REAL));
  *value = angle;

  pFM = new_FunctionalMatrix (2, 2, &_UFunction, value);

  return pFM;
}

REAL _Complex
_InvFunction (UINTEGER  i,
              UINTEGER  j,
              UINTEGER  nrows,
              UINTEGER  unused __attribute__((unused)),
              void     *unused2 __attribute__((unused)))
{
  return -(i == j) + 2.0/nrows + 0 * I; // if i != j -> 2/N, else -> -1 + 2/N
}

// Inversion about the average
FunctionalMatrix*
IAA (int nqubits)
{
  UINTEGER size = 2 << (nqubits - 1);

  return new_FunctionalMatrix (size, size, &_InvFunction, NULL);
}

REAL _Complex
_CustomGate (UINTEGER  i,
             UINTEGER  j,
             UINTEGER  nrows,
             UINTEGER  unused1 __attribute__((unused)),
             void     *matrix_2d)
{
  REAL _Complex *custom_matrix;

  custom_matrix = (REAL _Complex*) matrix_2d;

  return custom_matrix[i * nrows + j];
}

FunctionalMatrix*
CustomGate (REAL _Complex *matrix_2d,
            UINTEGER       nrows)
{
  return new_FunctionalMatrix (nrows, nrows, &_CustomGate, matrix_2d);
}

FunctionalMatrix*
PyCustomGate (REAL     *re_2d,
              REAL     *im_2d,
              UINTEGER  nrows,
              UINTEGER  size)
{
  UINTEGER i;
  REAL _Complex *matrix_2d = (REAL _Complex*) malloc (size * sizeof (REAL _Complex));

  for (i = 0; i < size; i++)
    matrix_2d[i] = re_2d[i] + im_2d[i] * I;

  return CustomGate (matrix_2d, nrows);
}

REAL _Complex
_CUFunction (UINTEGER  i,
             UINTEGER  j,
             UINTEGER  unused1 __attribute__((unused)),
             UINTEGER  unused2 __attribute__((unused)),
             void     *RawU)
{
  REAL _Complex val;
  int result = 1;
  FunctionalMatrix *U = (FunctionalMatrix*) RawU;

  if (i < (UINTEGER) rows (U) || j < (UINTEGER) columns (U))
    val = i == j;
  else
    result = getitem (U, i - rows (U), j - columns (U), &val);

  if (!result)
    printf ("Error getting element ("UINTEGER_STRING_FORMAT", "UINTEGER_STRING_FORMAT") from U gate\n", i - rows (U), j - columns (U));

  return val;
}

FunctionalMatrix*
CU (FunctionalMatrix *U)
{
  return new_FunctionalMatrix (rows (U) * 2, columns (U) * 2, &_CUFunction, U);
}

unsigned int
intLog2 (unsigned int N)
{
   N |= (N >> 1);
   N |= (N >> 2);
   N |= (N >> 4);
   N |= (N >> 8);
   N |= (N >> 16);
   return (__builtin_popcount(N) - 1);
}

int
getGateQubits (FunctionalMatrix *U)
{
  int nQubits = -1;

  if (U != NULL && rows (U) == columns (U)) {
    nQubits = (int) intLog2((unsigned int) rows (U));
  }

  return nQubits;
}

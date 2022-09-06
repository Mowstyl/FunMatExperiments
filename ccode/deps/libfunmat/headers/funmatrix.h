#ifndef __FUNMATRIX_H
#define __FUNMATRIX_H

#include <complex.h>
#include "precision.h"

struct FMatrix
{
  /* Number of rows */
  UINTEGER          r;
  /* Number of columns */
  UINTEGER          c;
  /* Function that, given (i, j, nrows, ncolumns, *argv)
    returns the value of the element (i, j) of
    the matrix */
  REAL _Complex   (*f)     (UINTEGER,
                            UINTEGER,
                            UINTEGER,
                            UINTEGER,
                            void*);
  /* Scalar number s that will be multiplied by the result of f(i, j) or multiplied by A op B */
  REAL _Complex     s;
  /* Pointer to matrix A in case an operation is going to be performed A op B */
  struct FMatrix   *A;
  /* Pointer to matrix B in case an operation is going to be performed A op B */
  struct FMatrix   *B;
  /* Operation to apply between the matrices.
      0 -> Matrix addition                A + B
      1 -> Matrix subtraction             A - B
      2 -> Matrix multiplication          A * B
      3 -> Entity-wise multiplication     A.* B
      4 -> Kronecker product              A ⊗ B
  */
  short             op;
  /* Whether the matrix has to be transposed or not */
  short             transpose;
  /* Whether the matrix has to be complex conjugated or not */
  short             conjugate;
  /* Whether the matrix is simple or you have to perform an operation */
  short             simple;
  /* Extra arguments to pass to the function f */
  void             *argv;
};

typedef struct FMatrix FunctionalMatrix;

struct DMatrixForTrace
{
  /* Density Matrix */
  FunctionalMatrix *m;
  /* Element to trace out */
  int               e;
};

typedef struct DMatrixForTrace _MatrixElem;

/* Constructor */
FunctionalMatrix* new_FunctionalMatrix (UINTEGER           n_rows,
                                        UINTEGER           n_columns,
                                        REAL _Complex    (*fun)      (UINTEGER,
                                                                      UINTEGER,
                                                                      UINTEGER,
                                                                      UINTEGER,
                                                                      void*),
                                        void              *argv);

/*
 * Get the element (i, j) from the matrix a, and return the result in
 * the address pointed by sol. If a 0 is returned, something went wrong.
 */
 __attribute__ ((pure))
int               getitem              (FunctionalMatrix  *a,
                                        UINTEGER           i,
                                        UINTEGER           j,
                                        REAL _Complex     *sol);

__attribute__ ((pure))
int               getitemaux           (FunctionalMatrix  *a,
                                        UINTEGER           i,
                                        UINTEGER           j,
                                        REAL              *re,
                                        REAL              *im);

/* Addition */
FunctionalMatrix* madd                 (FunctionalMatrix  *a,
                                        FunctionalMatrix  *b);

/* Subtraction */
FunctionalMatrix* msub                 (FunctionalMatrix  *a,
                                        FunctionalMatrix  *b);

/* Scalar product */
FunctionalMatrix* mprod                (REAL _Complex      r,
                                        FunctionalMatrix  *a);


FunctionalMatrix* mprodaux             (REAL               re,
                                        REAL               im,
                                        FunctionalMatrix  *a);

/* Scalar division */
FunctionalMatrix* mdiv                 (REAL _Complex      r,
                                        FunctionalMatrix  *a);

/* Scalar product */
FunctionalMatrix* mdivaux              (REAL               re,
                                        REAL               im,
                                        FunctionalMatrix  *a);

/* Matrix multiplication */
FunctionalMatrix* matmul               (FunctionalMatrix  *a,
                                        FunctionalMatrix  *b);

/* Entity-wise multiplication */
FunctionalMatrix* ewmul                (FunctionalMatrix  *a,
                                        FunctionalMatrix  *b);

/* Kronecker product */
FunctionalMatrix* kron                 (FunctionalMatrix  *a,
                                        FunctionalMatrix  *b);

/* Transpose */
FunctionalMatrix* transpose            (FunctionalMatrix  *m);

/* Hermitian transpose */
FunctionalMatrix* dagger               (FunctionalMatrix  *m);

int               rows                 (FunctionalMatrix  *m);
int               columns              (FunctionalMatrix  *m);

__attribute__ ((const))
UINTEGER          _GetElemIndex        (int                value,
                                        UINTEGER           position,
                                        int                bit);

__attribute__ ((pure))
REAL _Complex     _PartialTFunct       (UINTEGER           i,
                                        UINTEGER           j,
                                        UINTEGER           unused1 __attribute__((unused)),
                                        UINTEGER           unused2 __attribute__((unused)),
                                        void              *items);

/* Partial trace */
FunctionalMatrix* partial_trace        (FunctionalMatrix  *m,
                                        int                elem);

/*
 * Calculates the number of bytes added to a string
 * using the result of the sprintf function.
 */
__attribute__ ((const))
int               _bytes_added         (int                sprintfRe);

/* Gets the size in memory */
__attribute__ ((pure))
size_t            getMemory            (FunctionalMatrix *fm);

/* Print matrix */
__attribute__ ((pure))
char*             FM_toString          (FunctionalMatrix  *a);

#endif

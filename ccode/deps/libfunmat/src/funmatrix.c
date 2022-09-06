#ifdef __MINGW32__
#define __USE_MINGW_ANSI_STDIO 1
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include "funmatrix.h"

/* Constructor */
FunctionalMatrix*
new_FunctionalMatrix (UINTEGER        n_rows,
                      UINTEGER        n_columns,
                      REAL _Complex (*fun)      (UINTEGER,
                                                 UINTEGER,
                                                 UINTEGER,
                                                 UINTEGER,
                                                 void*),
                      void           *argv)
{
  FunctionalMatrix* pFM = (FunctionalMatrix*) malloc (sizeof (FunctionalMatrix));

  pFM->r = n_rows;
  pFM->c = n_columns;
  pFM->f = fun;
  pFM->s = 1;
  pFM->transpose = 0;
  pFM->conjugate = 0;
  pFM->simple = 1;
  pFM->argv = argv;

  return pFM;
}

/* Get the element (i, j) from the matrix a */
int
getitem (FunctionalMatrix *a,
         UINTEGER          i,
         UINTEGER          j,
         REAL _Complex    *sol)
{
  unsigned int k;
  UINTEGER aux;
  int result = 1;
  REAL _Complex aux1 = 0,
                aux2 = 0;

  if (i < a->r && j < a->c)
    {
      if (a->transpose)
        {
          aux = i;
          i = j;
          j = aux;
        }

      if (a->simple)
        {
          *sol = a->f (i, j, a->r, a->c, a->argv);
        }
      else
        {
          switch (a->op)
            {
            case 0: /* Matrix addition */
              if (getitem (a->A, i, j, &aux1) && getitem (a->B, i, j, &aux2))
                {
                  *sol = aux1 + aux2;
                }
              else
                {
                  printf ("Error while operating!\n");
                  result = 0;
                  sol = NULL;
                }
              break;

            case 1: /* Matrix subtraction */
              if (getitem (a->A, i, j, &aux1) && getitem (a->B, i, j, &aux2))
                {
                  *sol = aux1 - aux2;
                }
              else
                {
                  printf ("Error while operating!\n");
                  result = 0;
                  sol = NULL;
                }
              break;

            case 2: /* Matrix multiplication  */
              *sol = 0;
              for (k = 0; k < a->A->c; k++)
              {
                if (getitem (a->A, i, k, &aux1) && getitem (a->B, k, j, &aux2))
                  {
                    *sol = *sol + aux1 * aux2;
                  }
                else
                  {
                    printf ("Error while operating!\n");
                    result = 0;
                    sol = NULL;
                  }
              }
              break;

            case 3: /* Entity-wise multiplication */
              if (getitem (a->A, i, j, &aux1) && getitem (a->B, i, j, &aux2))
                {
                  *sol = aux1 * aux2;
                }
              else
                {
                  printf ("Error while operating!\n");
                  result = 0;
                  sol = NULL;
                }
              break;

            case 4: /* Kronecker product */
              if (getitem (a->A, i/a->B->r, j/a->B->c, &aux1) && getitem (a->B, i%a->B->r, j%a->B->c, &aux2))
                {
                  *sol = aux1 * aux2;
                }
              else
                {
                  printf ("Error while operating!\n");
                  result = 0;
                  sol = NULL;
                }
              break;

            default:
              printf ("Unknown option: %d\n", a->op);
              result = 0;
              sol = NULL;
            }
        }

        if (result && a->conjugate)
          *sol = CONJ (*sol);
    }
  else
    {
      printf ("(" UINTEGER_STRING_FORMAT ", " UINTEGER_STRING_FORMAT ") is out of bounds!\n Matrix dimensions: (" UINTEGER_STRING_FORMAT ", " UINTEGER_STRING_FORMAT ")\n", i, j, a->r, a->c);
      result = 0;
      sol = NULL;
    }

  if (result)
    *sol = *sol * a->s;

  return result;
}

int
getitemaux (FunctionalMatrix *a,
            UINTEGER          i,
            UINTEGER          j,
            REAL              *re,
            REAL              *im)
{
  int result;
  REAL _Complex aux = 0;

  result = getitem (a, i, j, &aux);
  if (result) {
    *re = RE (aux);
    *im = IM (aux);
  }

  return result;
}

/* Addition */
FunctionalMatrix*
madd (FunctionalMatrix *a,
      FunctionalMatrix *b)
{
  FunctionalMatrix* pFM = NULL;

  if (a->r == b->r && a->c == b->c)
    { /* if the dimensions allign (nxm .* nxm)*/
      pFM = (FunctionalMatrix*) malloc (sizeof (FunctionalMatrix));
      pFM->r = a->r;
      pFM->c = a->c;
      pFM->s = 1;
      pFM->A = a;
      pFM->B = b;
      pFM->op = 0;
      pFM->simple = 0;
    }

  return pFM;
}

/* Subtraction */
FunctionalMatrix*
msub (FunctionalMatrix *a,
      FunctionalMatrix *b)
{
  FunctionalMatrix* pFM = NULL;

  if (a->r == b->r && a->c == b->c)
    { /* if the dimensions allign (nxm .* nxm)*/
      pFM = (FunctionalMatrix*) malloc (sizeof (FunctionalMatrix));
      pFM->r = a->r;
      pFM->c = a->c;
      pFM->s = 1;
      pFM->A = a;
      pFM->B = b;
      pFM->op = 0;
      pFM->simple = 0;
    }

  return pFM;
}

/* Scalar product */
FunctionalMatrix*
mprod (REAL _Complex     r,
       FunctionalMatrix *a)
{
  FunctionalMatrix* pFM = (FunctionalMatrix*) malloc (sizeof (FunctionalMatrix));

  pFM->r = a->r;
  pFM->c = a->c;
  pFM->f = a->f;
  pFM->s = a->s * r;
  pFM->A = a->A;
  pFM->B = a->B;
  pFM->op = a->op;
  pFM->simple = a->simple;
  pFM->argv = a->argv;

  return pFM;
}

FunctionalMatrix*
mprodaux (REAL              re,
          REAL              im,
          FunctionalMatrix *a)
{
  REAL _Complex r = re + im * I;

  return mprod (r, a);
}

/* Scalar division */
FunctionalMatrix*
mdiv(REAL _Complex     r,
     FunctionalMatrix *a)
{
  FunctionalMatrix* pFM = (FunctionalMatrix*) malloc (sizeof (FunctionalMatrix));

  pFM->r = a->r;
  pFM->c = a->c;
  pFM->f = a->f;
  pFM->s = a->s / r;
  pFM->A = a->A;
  pFM->B = a->B;
  pFM->op = a->op;
  pFM->simple = a->simple;
  pFM->argv = a->argv;

  return pFM;
}

FunctionalMatrix*
mdivaux (REAL              re,
         REAL              im,
         FunctionalMatrix *a)
{
  REAL _Complex r = re + im * I;

  return mdiv (r, a);
}

/* Matrix multiplication */
FunctionalMatrix*
matmul (FunctionalMatrix *a,
        FunctionalMatrix *b)
{
  FunctionalMatrix* pFM = NULL;

  if (a->c == b->r)
    { /* if the dimensions allign (uxv * vxw) */
      pFM = (FunctionalMatrix*) malloc (sizeof (FunctionalMatrix));
      pFM->r = a->r;
      pFM->c = b->c;
      pFM->s = 1;
      pFM->A = a;
      pFM->B = b;
      pFM->op = 2;
      pFM->simple = 0;
    }

  return pFM;
}

/* Entity-wise multiplication */
FunctionalMatrix*
ewmul (FunctionalMatrix *a,
       FunctionalMatrix *b)
{
  FunctionalMatrix* pFM = NULL;

  if (a->r == b->r && a->c == b->c)
    { /* if the dimensions allign (nxm .* nxm)*/
      pFM = (FunctionalMatrix*) malloc (sizeof (FunctionalMatrix));
      pFM->r = a->r;
      pFM->c = a->c;
      pFM->s = 1;
      pFM->A = a;
      pFM->B = b;
      pFM->op = 3;
      pFM->simple = 0;
    }
  else if (a->r == 1 && b->c == 1)
    { /* row .* column */
      pFM = matmul(b, a);
    }
  else if (b->r == 1 && a->c == 1)
    {
      pFM = matmul(a, b);
    }

  return pFM;
}

/* Kronecker product */
FunctionalMatrix*
kron(FunctionalMatrix *a,
     FunctionalMatrix *b)
{
  FunctionalMatrix* pFM = (FunctionalMatrix*) malloc (sizeof (FunctionalMatrix));

  pFM->r = a->r * b->r;
  pFM->c = a->c * b->c;
  pFM->s = 1;
  pFM->A = a;
  pFM->B = b;
  pFM->op = 4;
  pFM->simple = 0;

  return pFM;
}

/* Transpose */
FunctionalMatrix*
transpose (FunctionalMatrix  *m)
{
  FunctionalMatrix* pFM = (FunctionalMatrix*) malloc (sizeof (FunctionalMatrix));

  pFM->r = m->r;
  pFM->c = m->c;
  pFM->f = m->f;
  pFM->s = m->s;
  pFM->A = m->A;
  pFM->B = m->B;
  pFM->op = m->op;
  pFM->transpose = !m->transpose;
  pFM->conjugate = m->conjugate;
  pFM->simple = 0;

  return pFM;
}

/* Hermitian transpose */
FunctionalMatrix*
dagger (FunctionalMatrix  *m)
{
  FunctionalMatrix* pFM = (FunctionalMatrix*) malloc (sizeof (FunctionalMatrix));

  pFM->r = m->r;
  pFM->c = m->c;
  pFM->f = m->f;
  pFM->s = m->s;
  pFM->A = m->A;
  pFM->B = m->B;
  pFM->op = m->op;
  pFM->transpose = !m->transpose;
  pFM->conjugate = !m->conjugate;
  pFM->simple = 0;

  return pFM;
}

UINTEGER
_GetElemIndex (int      value,
               UINTEGER position,
               int      bit)
{
  UINTEGER index = 0,
           aux = 1;

  if ((value == 0 || value == 1) && bit >= 0)
    {
      if (bit != 0)
        aux = (UINTEGER) (2 << (bit - 1));
      index = position % aux + (position / aux) * (aux << 1) + value * aux;
    }

  return index;
}

REAL _Complex
_PartialTFunct (UINTEGER  i,
                UINTEGER  j,
                UINTEGER  unused1 __attribute__((unused)),
                UINTEGER  unused2 __attribute__((unused)),
                void     *items)
{
  REAL _Complex sol = 0,
                aux = 0;
  _MatrixElem *me;

  if (items != NULL)
    {
      me = (_MatrixElem*) items;

      if (getitem (me->m, _GetElemIndex (0, i, me->e), _GetElemIndex (0, j, me->e), &sol) &&
          getitem (me->m, _GetElemIndex (1, i, me->e), _GetElemIndex (1, j, me->e), &aux))
        {
          sol = sol + aux;
        }
    }

  return sol;
}

/* Partial trace */
FunctionalMatrix*
partial_trace (FunctionalMatrix  *m, int elem)
{
  FunctionalMatrix *pt = NULL;
  _MatrixElem *me = NULL;

  if (m != NULL && m->r == m->c && elem >= 0)
    {
      me = (_MatrixElem*) malloc (sizeof (_MatrixElem));
      me->m = m;
      me->e = elem;
      pt = new_FunctionalMatrix (m->r >> 1, m->c >> 1, _PartialTFunct, me);
    }

  return pt;
}


int
rows (FunctionalMatrix  *m)
{
  return m->r;
}

int
columns (FunctionalMatrix  *m)
{
  return m->c;
}

int
_bytes_added (int sprintfRe)
{
    return (sprintfRe > 0) ? sprintfRe : 0;
}

/* Gets the size in memory */
size_t
getMemory (FunctionalMatrix *m)
{
  size_t total;

  total = sizeof(*m);
  if (!m->simple)
    {
      total += getMemory(m->A);
      total += getMemory(m->B);
    }

  return total;
}

/* Print matrix */
char*
FM_toString (FunctionalMatrix *a)
{
  char *text;
  REAL _Complex it;
  UINTEGER i, j;
  int length = 0;
  const UINTEGER MAX_BUF = a->r * a->c * (2 * (DECIMAL_PLACES + 7) + 2) + 2;
  // numero de elementos (r * c) multiplicado por numero de cifras significativas establecidas para cada numero
  // por 2 (son complejos) mas 7 (1 del signo, otro del . y 5 del exponente e-001) mas 2, uno de la i y otro del
  // espacio/;/] que hay despues de cada numero. Al final se suman 2, uno para el corchete inicial y otro para \0.

  text = (char*) malloc (MAX_BUF);
  length += _bytes_added (snprintf (text + length, MAX_BUF - length, "["));
  for (i = 0; i < a->r; i++)
    {
      for (j = 0; j < a->c; j++)
        {
          if (getitem (a, i, j, &it))
            {
              if (IM (it) >= 0)
                {
                  length += _bytes_added (snprintf (text + length,
                                                    MAX_BUF - length,
                                                    REAL_STRING_FORMAT "+" REAL_STRING_FORMAT "i",
                                                    RE (it),
                                                    IM (it)));
                }
              else
                {
                  length += _bytes_added (snprintf (text + length,
                                                    MAX_BUF - length,
                                                    REAL_STRING_FORMAT REAL_STRING_FORMAT "i",
                                                    RE (it),
                                                    IM (it)));
                }
            }
          else
          {
            length += _bytes_added (snprintf (text + length, MAX_BUF - length, "ERR"));
          }
          if (j < a->c - 1)
            length += _bytes_added (snprintf (text + length, MAX_BUF - length, " "));
        }

      if (i < a->r - 1)
        length += _bytes_added (snprintf (text + length, MAX_BUF - length, ";"));
    }
  length += _bytes_added (snprintf (text + length, MAX_BUF - length, "]"));
  *(text+length) = '\0';

  return text;
}

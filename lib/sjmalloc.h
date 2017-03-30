// Author: Hou, Sian - sianhou1987@outlook.com

#ifndef SJI_SJMALLOC_H
#define SJI_SJMALLOC_H

#include "sjinc.h"

//! function for allocate the 1 demension array and initialize with zero
void *sjalloc1d(int n1, int size);

//! function for allocate the 2 demension array and initialize with zero
void **sjalloc2d(int n2, int n1, int size);

//! function for allocate the 3 demension array and initialize with zero
void ***sjalloc3d(int n3, int n2, int n1, int size);

//! function for free the 1 demension array
void sjfree1d(void *p);

void sjcheckfree1d(void *p);

//! function for free the 2 demension array
void sjfree2d(void **p);

void sjcheckfree2d(void **p);

//! function for free the 3 demension array
void sjfree3d(void ***p);

void sjcheckfree3d(void ***p);

//! Fast allocate the 1 demension using macro
#define sjmilloc1d(n1) (int *)sjalloc1d(n1, sizeof(int))
#define sjmflloc1d(n1) (float *)sjalloc1d(n1, sizeof(float))
#define sjmdlloc1d(n1) (double *)sjalloc1d(n1, sizeof(double))
#define sjmclloc1d(n1) (fcomplex *)sjalloc1d(n1, sizeof(fcomplex))
#define sjmzlloc1d(n1) (dcomplex *)sjalloc1d(n1, sizeof(dcomplex))

//! Fast allocate the 2 demension using macro
#define sjmilloc2d(n2, n1) (int **)sjalloc2d(n2, n1,  sizeof(int))
#define sjmflloc2d(n2, n1) (float **)sjalloc2d(n2, n1, sizeof(float))
#define sjmdlloc2d(n2, n1) (double **)sjalloc2d(n2, n1, sizeof(double))
#define sjmclloc2d(n2, n1) (fcomplex **)sjalloc2d(n2, n1, sizeof(fcomplex))
#define sjmzlloc2d(n2, n1) (dcomplex **)sjalloc2d(n2, n1, sizeof(dcomplex))

//! Fast allocate the 3 demension using macro
#define sjmilloc3d(n3, n2, n1) (int ***)sjalloc3d(n3, n2, n1, sizeof(int))
#define sjmflloc3d(n3, n2, n1) (float ***)sjalloc3d(n3, n2, n1, sizeof(float))
#define sjmdlloc3d(n3, n2, n1) (double ***)sjalloc3d(n3, n2, n1, sizeof(double))
#define sjmclloc3d(n3, n2, n1) (fcomplex ***)sjalloc3d(n3, n2, n1, sizeof(fcomplex))
#define sjmzlloc3d(n3, n2, n1) (dcomplex ***)sjalloc3d(n3, n2, n1, sizeof(dcomplex))

//! Fast free the 1 demension using macro
#define sjmfree1d(p) sjfree1d((void *) p)
#define sjmfree2d(p) sjfree2d((void **) p)
#define sjmfree3d(p) sjfree3d((void ***) p)

#endif //SJI_SJMALLOC_H

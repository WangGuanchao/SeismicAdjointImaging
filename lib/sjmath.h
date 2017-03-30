// Author: Hou, Sian - sianhou1987@outlook.com

#ifndef GJI_SJMATH_H
#define GJI_SJMATH_H

#include "sjinc.h"

#define SJMMAX(a, b) (a > b ? a : b)
#define SJMMIN(a, b) (a < b ? a : b)

int sjfindabsmaxf(float *a, int n);

void sjguasssmoothf2d(float **z, int n2, int n1, float alpha, int length, float **x);

void sjfilter2d(float **z, int n2, int n1, float **x, char *mode);

float sjvecdotf(int n, float a, float *x, float *y);

void sjvecaddf(float *z, int n, float a, float *x, float b, float *y); //! c = a + b

void sjvecsubf(float *z, int n, float a, float *x, float b, float *y); //! c = a - b

//! z[] = a * x[] * y[]
void sjvecmulf(float *z, int n, float a, float *y, float *x);

void sjvecdivf(float *z, int n, float a, float *x, float *y, float ep); //! c = a / (b + ep)

//! Z[] = 0.0f
void sjveczerof(float *z, int n);

void sjcgdirection(float *d, int n, float *g1, float *g0, int iter);

float sjcglength(int n, float *s, float *x, float err, int iter);

float sjcgbeta(int n, float *cg, float *g1, float *g0, int iter);

void sjcgsolver(float *z, int n, float *cg, float *g1, float *g0, int iter);

#endif //GJI_SJMATH_H

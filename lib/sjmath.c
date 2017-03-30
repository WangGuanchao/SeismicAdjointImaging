// Author: Hou, Sian - sianhou1987@outlook.com

#include "sjmath.h"

int sjmaxabscorf(float *pobs, float *pcal, int n, int maxshift) {

}

void sjguasssmoothf2d(float **z, int n2, int n1, float alpha, int length, float **x) {
    int ii, jj, ix, iz;

    float **g = (float **) sjalloc2d(2 * length + 1, 2 * length + 1, sizeof(float));
    for (ix = -length; ix <= length; ++ix)
        for (iz = -length; iz <= length; iz++)
            g[ix + length][iz + length] = alpha * expf(-alpha * (ix * ix + iz * iz)) / pi;
    float **p = (float **) sjalloc2d(n2 + 2 * length, n1 + 2 * length, sizeof(float));
    sjextend2d(p, n2, n1, length, length, length, length, x);
    for (ii = 0; ii < n2; ++ii) {
        for (jj = 0; jj < n1; ++jj) {
            z[ii][jj] = 0.0f;
            for (ix = -length; ix <= length; ++ix)
                for (iz = -length; iz <= length; iz++)
                    z[ii][jj] += g[ix + length][iz + length] * p[ii + length + ix][jj + length + iz];
        }
    }
    sjmfree2d(g);
    sjmfree2d(p);
}

void sjfilter2d(float **z, int n2, int n1, float **x, char *mode) {
    //! z = filter2(x)
    int ix, iz;
    float **p = (float **) sjalloc2d(n2, n1, sizeof(float));
    memcpy(p[0], x[0], n2 * n1 * sizeof(float));
    if (strcmp(mode, "laplace") == 0) {
        for (ix = 2; ix < n2 - 2; ++ix)
            for (iz = 2; iz < n1 - 2; ++iz)
                z[ix][iz] = -4.0f * p[ix][iz] + p[ix - 1][iz] + p[ix + 1][iz] + p[ix][iz - 1] + p[ix][iz + 1];
    }
    if (strcmp(mode, "guass5") == 0) {
        float a0 = 41.0f / 273.0f;
        float a1 = 26.0f / 273.0f;
        float a2 = 16.0f / 273.0f;
        float a3 = 7.0f / 273.0f;
        float a4 = 4.0f / 273.0f;
        float a5 = 1.0f / 273.0f;
        for (ix = 2; ix < n2 - 2; ++ix)
            for (iz = 2; iz < n1 - 2; ++iz)
                z[ix][iz] = a0 * (p[ix][iz])
                            + a1 * (p[ix + 1][iz] + p[ix - 1][iz] + p[ix][iz + 1] + p[ix][iz - 1])
                            + a2 * (p[ix + 1][iz + 1] + p[ix - 1][iz + 1] + p[ix + 1][iz - 1] + p[ix - 1][iz - 1])
                            + a3 * (p[ix + 2][iz] + p[ix - 2][iz] + p[ix][iz + 2] + p[ix][iz - 2])
                            + a4 * (p[ix + 2][iz + 1] + p[ix + 1][iz + 2] + p[ix - 2][iz + 1] + p[ix - 1][iz + 2] +
                                    p[ix + 2][iz - 1] + p[ix + 1][iz - 2] + p[ix + 2][iz - 1] + p[ix + 1][iz - 2])
                            + a5 * (p[ix + 2][iz + 2] + p[ix - 2][iz + 2] + p[ix + 2][iz - 2] + p[ix - 2][iz - 2]);
    }
    sjmfree2d(p);
}

int sjfindabsmaxf(float *a, int n) {
    //! Find the index of a vector' absolute maxmium
    int ii, index = 0;
    for (ii = 0; ii < n; ++ii)
        if (fabs(a[ii]) > fabs(a[index])) index = ii;
    return index;
}

float sjvecdotf(int n, float a, float *x, float *y) {
    //! z = sum( a * x[] * b[] )
    int ii;
    float z = 0.0f;
    for (ii = 0; ii < n; ++ii)
        z += a * x[ii] * y[ii];
    return z;
}

void sjvecaddf(float *z, int n, float a, float *x, float b, float *y) {
    //! z[] = a*x[] + b*y[]
    int ii;
    for (ii = 0; ii < n; ++ii)
        z[ii] = a * x[ii] + b * y[ii];
}

void sjvecsubf(float *z, int n, float a, float *x, float b, float *y) {
    //! z[] = a*x[] - b*y[]
    int ii;
    for (ii = 0; ii < n; ++ii)
        z[ii] = a * x[ii] - b * y[ii];
}

void sjvecmulf(float *z, int n, float a, float *y, float *x) {
    //! z[] = a * x[] * y[]
    int ii;
    for (ii = 0; ii < n; ++ii)
        z[ii] = a * x[ii] * y[ii];
}

void sjvecdivf(float *z, int n, float a, float *x, float *y, float ep) {
    //! z[] = a * x[] / (b[] + ep)
    int ii;
    for (ii = 0; ii < n; ++ii)
        z[ii] = a * x[ii] / (y[ii] + ep);
}

void sjveczerof(float *z, int n) {
    memset(z, 0, n * sizeof(float));
}

void sjcgdirection(float *d, int n, float *g1, float *g0, int iter) {
    //! by Guanchao wang
    int ii;
    float beta;
    if (iter == 0) {
        for (ii = 0; ii < n; ii++) {
            d[ii] = -1.0f * g1[ii];
        }
    } else {
        beta = (sjvecdotf(n, 1.0, g1, g1) - sjvecdotf(n, 1.0, g0, g1)) / sjvecdotf(n, 1.0, g0, g0);
        sjvecaddf(d, n, -1.0f, g1, beta, d);
    }
    for (ii = 0; ii < n; ++ii) {
        g0[ii] = g1[ii];
    }
}

float sjcglength(int n, float *s, float *x, float err, int iter) {
    //! by Guanchao wang
    int index_s, index_x;
    if (iter == 0) {
        return 0.0001f;
    } else {
        index_s = sjfindabsmaxf(s, n);
        index_x = sjfindabsmaxf(x, n);
        return err * fabsf(x[index_x]) / fabsf(s[index_s]);
    }
}

float sjcgbeta(int n, float *cg, float *g1, float *g0, int iter) {
    //! by Pengliang Yang
    if (iter == 0) {
        return 0.0f;
    } else {
        int ii;
        float a, b, c;
        a = b = c = 0;
        for (ii = 0; ii < n; ++ii) {
            a += g1[ii] * (g1[ii] - g0[ii]);
            b += cg[ii] * (g1[ii] - g0[ii]);
            c += g1[ii] * g1[ii];
        }

        float beta_HS = (fabsf(b) > 0.0f) ? (a / b) : 0.0f;
        float beta_DY = (fabsf(b) > 0.0f) ? (c / b) : 0.0f;
        return SJMMAX(0.0f, SJMMIN(beta_HS, beta_DY));
    }
}

void sjcgsolver(float *z, int n, float *cg, float *g1, float *g0, int iter) {

    float alpha, beta;

    //! Calculate CG direction
    if(iter==0) {
        beta = 0.0f;
        sjvecaddf(cg, n, -1.0f, g1, beta, cg);
    } else {
        beta = sjcgbeta(n, cg, g1, g0, iter);
        sjvecaddf(cg, n, -1.0f, g1, beta, cg);
    }
    memcpy(g0, g1, n * sizeof(float));

    //! Calculate CG step size
    alpha = sjcglength(n, cg, z, 0.1, iter);

    //! Update model
    sjvecaddf(z, n, 1.0f, z, alpha, cg);
}

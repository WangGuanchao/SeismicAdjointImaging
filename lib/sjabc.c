// Author: Hou, Sian - sianhou1987@outlook.com

#include "sjabc.h"

//! Initialize two order time hybrid abc for 1 demension simulation
int sjinitthabc1d(float *vp, float *vs, float ds, float dt, int nxb, float *gxl, float *gxr) {
    float v, b = 0.5, beta = 1.0;
    float qxz, rxz, qt, rt, qxzt, rxzt;

    //! Left bounadry condition
    v = vp[0] * dt / ds;
    qxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
    qt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
    qxzt = b / (b - 1.0f);

    v = vs[0] * dt / ds;
    rxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
    rt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
    rxzt = b / (b - 1.0f);

    gxl[0] = -(qxz + rxz);
    gxl[1] = -(qxz * rxz);
    gxl[2] = -(qt + rt);
    gxl[3] = -(qxz * rt + qt * rxz + qxzt + rxzt);
    gxl[4] = -(qxz * rxzt + rxz * qxzt);
    gxl[5] = -(qt * rt);
    gxl[6] = -(qt * rxzt + rt * qxzt);
    gxl[7] = -(qxzt * rxzt);

    //! Right bounadry condition
    v = vp[nxb - 1] * dt / ds;
    qxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
    qt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
    qxzt = b / (b - 1.0f);

    v = vs[nxb - 1] * dt / ds;
    rxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
    rt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
    rxzt = b / (b - 1.0f);

    gxr[0] = -(qxz + rxz);
    gxr[1] = -(qxz * rxz);
    gxr[2] = -(qt + rt);
    gxr[3] = -(qxz * rt + qt * rxz + qxzt + rxzt);
    gxr[4] = -(qxz * rxzt + rxz * qxzt);
    gxr[5] = -(qt * rt);
    gxr[6] = -(qt * rxzt + rt * qxzt);
    gxr[7] = -(qxzt * rxzt);

    return 1;
}

//! Apply two order time hybrid abc for 1 demension simulation
int sjapplythabc1d(float *fp, float *cp, float *pp, float *gxl, float *gxr, int nxb, int nb, int marg) {
    int ii, ix;
    float temp;

    for (ii = 0; ii <= nb; ii++) {
        //! Left & right bounadry condition
        ix = marg + nb - ii;
        temp = gxl[0] * fp[ix + 1] + gxl[1] * fp[ix + 2]
               + gxl[2] * cp[ix + 0] + gxl[3] * cp[ix + 1] + gxl[4] * cp[ix + 2]
               + gxl[5] * pp[ix + 0] + gxl[6] * pp[ix + 1] + gxl[7] * pp[ix + 2];
        fp[ix] = fp[ix] * ((float) (nb - 1 - ii)) / ((float) nb) + temp * ((float) (1 + ii)) / ((float) nb);

        ix = nxb - marg - nb + ii - 1;
        temp = gxr[0] * fp[ix - 1] + gxr[1] * fp[ix - 2]
               + gxr[2] * cp[ix + 0] + gxr[3] * cp[ix - 1] + gxr[4] * cp[ix - 2]
               + gxr[5] * pp[ix + 0] + gxr[6] * pp[ix - 1] + gxr[7] * pp[ix - 2];
        fp[ix] = fp[ix] * ((float) (nb - 1 - ii)) / ((float) nb) + temp * ((float) (1 + ii)) / ((float) nb);
    }

    return 1;
}

//! Initialize one order time hybrid abc for 2 demension simulation
int sjinitohabc2d(float **vp, float **vs, float ds, float dt, int nxb, int nzb,
                  float **gxl, float **gxr, float **gzu, float **gzb) {
    int ix, iz;
    float r, b = 0.5, beta1 = 0, beta2 = 0, beta = 1.0;

    //! Left bounadry condition*/
    for (iz = 0; iz < nzb; iz++) {
        beta1 = 1.0;
        beta2 = vp[0][iz] / vs[0][iz];
        beta = (beta1 + beta2) / 2.0f;
        r = vp[0][iz] * dt / ds;

        gxl[iz][0] = -(b * (beta + r) - r) / (beta + r) / (1.0f - b);
        gxl[iz][1] = -(b * (beta + r) - beta) / (beta + r) / (1.0f - b);
        gxl[iz][2] = -b / (b - 1.0f);
    }

    //! Right bounadry condition
    for (iz = 0; iz < nzb; iz++) {
        beta1 = 1.0;
        beta2 = vp[nxb - 1][iz] / vs[nxb - 1][iz];
        beta = (beta1 + beta2) / 2.0f;
        r = vp[nxb - 1][iz] * dt / ds;

        gxr[iz][0] = -(b * (beta + r) - r) / (beta + r) / (1.0f - b);
        gxr[iz][1] = -(b * (beta + r) - beta) / (beta + r) / (1.0f - b);
        gxr[iz][2] = -b / (b - 1.0f);
    }

    //! Upper bounadry condition
    for (ix = 0; ix < nxb; ix++) {
        beta1 = 1.0;
        beta2 = vp[ix][0] / vs[ix][0];
        beta = (beta1 + beta2) / 2.0f;
        r = vp[ix][0] * dt / ds;

        gzu[ix][0] = -(b * (beta + r) - r) / (beta + r) / (1.0f - b);
        gzu[ix][1] = -(b * (beta + r) - beta) / (beta + r) / (1.0f - b);
        gzu[ix][2] = -b / (b - 1.0f);
    }

    //! Below bounadry condition
    for (ix = 0; ix < nxb; ix++) {
        beta1 = 1.0;
        beta2 = vp[ix][nzb - 1] / vs[ix][nzb - 1];
        beta = (beta1 + beta2) / 2.0f;
        r = vp[ix][nzb - 1] * dt / ds;

        gzb[ix][0] = -(b * (beta + r) - r) / (beta + r) / (1.0f - b);
        gzb[ix][1] = -(b * (beta + r) - beta) / (beta + r) / (1.0f - b);
        gzb[ix][2] = -b / (b - 1.0f);
    }

    return 1;
}

//! Apply one order time hybrid abc for 2 demension simulation
int sjapplyohabc2d(float **fp, float **cp, float **gxl, float **gxr, float **gzu, float **gzb,
                   int nxb, int nzb, int nb, int marg) {
    int ii, ix, iz;
    float temp;

    for (ii = 0; ii <= nb; ii++) {
        //! Left & right bounadry condition
        for (iz = marg + nb - ii + 1; iz < nzb - marg - nb + ii - 1; iz++) {
            ix = marg + nb - ii;
            temp = gxl[iz][0] * fp[ix + 1][iz] + gxl[iz][1] * cp[ix + 0][iz] + gxl[iz][2] * cp[ix + 1][iz];
            fp[ix][iz] = fp[ix][iz] * ((float) (nb - ii)) / ((float) nb) + temp * ((float) (ii)) / ((float) nb);

            ix = nxb - marg - nb + ii - 1;
            temp = gxr[iz][0] * fp[ix - 1][iz] + gxr[iz][1] * cp[ix - 0][iz] + gxr[iz][2] * cp[ix - 1][iz];
            fp[ix][iz] = fp[ix][iz] * ((float) (nb - ii)) / ((float) nb) + temp * ((float) (ii)) / ((float) nb);
        }

        //! Upper & below bounadry condition
        for (ix = marg + nb - ii + 1; ix < nxb - marg - nb + ii - 1; ix++) {
            iz = marg + nb - ii;
            temp = gzu[ix][0] * fp[ix][iz + 1] + gzu[ix][1] * cp[ix][iz + 0] + gzu[ix][2] * cp[ix][iz + 1];
            fp[ix][iz] = fp[ix][iz] * ((float) (nb - ii)) / ((float) nb) + temp * ((float) (ii)) / ((float) nb);

            iz = nzb - marg - nb + ii - 1;
            temp = gzb[ix][0] * fp[ix][iz - 1] + gzb[ix][1] * cp[ix][iz - 0] + gzb[ix][2] * cp[ix][iz - 1];
            fp[ix][iz] = fp[ix][iz] * ((float) (nb - ii)) / ((float) nb) + temp * ((float) (ii)) / ((float) nb);
        }

        //! Corner
        ix = marg + nb - ii;
        iz = marg + nb - ii;
        temp = 0.5f * (gxl[iz][0] * fp[ix + 1][iz] + gxl[iz][1] * cp[ix + 0][iz] + gxl[iz][2] * cp[ix + 1][iz]) +
               0.5f * (gzu[ix][0] * fp[ix][iz + 1] + gzu[ix][1] * cp[ix][iz + 0] + gzu[ix][2] * cp[ix][iz + 1]);
        fp[ix][iz] = fp[ix][iz] * ((float) (nb - ii)) / ((float) nb) + temp * ((float) (ii)) / ((float) nb);

        ix = nxb - marg - nb + ii - 1;
        iz = marg + nb - ii;
        temp = 0.5f * (gxr[iz][0] * fp[ix - 1][iz] + gxr[iz][1] * cp[ix - 0][iz] + gxr[iz][2] * cp[ix - 1][iz]) +
               0.5f * (gzu[ix][0] * fp[ix][iz + 1] + gzu[ix][1] * cp[ix][iz + 0] + gzu[ix][2] * cp[ix][iz + 1]);
        fp[ix][iz] = fp[ix][iz] * ((float) (nb - ii)) / ((float) nb) + temp * ((float) (ii)) / ((float) nb);

        ix = marg + nb - ii;
        iz = nzb - marg - nb + ii - 1;
        temp = 0.5f * (gxl[iz][0] * fp[ix + 1][iz] + gxl[iz][1] * cp[ix + 0][iz] + gxl[iz][2] * cp[ix + 1][iz]) +
               0.5f * (gzb[ix][0] * fp[ix][iz - 1] + gzb[ix][1] * cp[ix][iz - 0] + gzb[ix][2] * cp[ix][iz - 1]);
        fp[ix][iz] = fp[ix][iz] * ((float) (nb - ii)) / ((float) nb) + temp * ((float) (ii)) / ((float) nb);

        ix = nxb - marg - nb + ii - 1;
        iz = nzb - marg - nb + ii - 1;
        temp = 0.5f * (gxr[iz][0] * fp[ix - 1][iz] + gxr[iz][1] * cp[ix - 0][iz] + gxr[iz][2] * cp[ix - 1][iz]) +
               0.5f * (gzb[ix][0] * fp[ix][iz - 1] + gzb[ix][1] * cp[ix][iz - 0] + gzb[ix][2] * cp[ix][iz - 1]);
        fp[ix][iz] = fp[ix][iz] * ((float) (nb - ii)) / ((float) nb) + temp * ((float) (ii)) / ((float) nb);
    }

    return 1;
}

//! Initialize two order time hybrid abc for 2 demension simulation
int sjinitthabc2d(float **vp, float **vs, float ds, float dt, int nxb, int nzb,
                  float **gxl, float **gxr, float **gzu, float **gzb) {
    int ix, iz;
    float v, b = 0.5, beta = 1.0;
    float qxz, rxz, qt, rt, qxzt, rxzt;

    //! Left bounadry condition
    for (iz = 0; iz < nzb; iz++) {
        v = vp[0][iz] * dt / ds;
        qxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
        qt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
        qxzt = b / (b - 1.0f);

        v = vs[0][iz] * dt / ds;
        rxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
        rt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
        rxzt = b / (b - 1.0f);

        gxl[iz][0] = -(qxz + rxz);
        gxl[iz][1] = -(qxz * rxz);
        gxl[iz][2] = -(qt + rt);
        gxl[iz][3] = -(qxz * rt + qt * rxz + qxzt + rxzt);
        gxl[iz][4] = -(qxz * rxzt + rxz * qxzt);
        gxl[iz][5] = -(qt * rt);
        gxl[iz][6] = -(qt * rxzt + rt * qxzt);
        gxl[iz][7] = -(qxzt * rxzt);
    }

    //! Right bounadry condition
    for (iz = 0; iz < nzb; iz++) {
        v = vp[nxb - 1][iz] * dt / ds;
        qxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
        qt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
        qxzt = b / (b - 1.0f);

        v = vs[nxb - 1][iz] * dt / ds;
        rxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
        rt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
        rxzt = b / (b - 1.0f);

        gxr[iz][0] = -(qxz + rxz);
        gxr[iz][1] = -(qxz * rxz);
        gxr[iz][2] = -(qt + rt);
        gxr[iz][3] = -(qxz * rt + qt * rxz + qxzt + rxzt);
        gxr[iz][4] = -(qxz * rxzt + rxz * qxzt);
        gxr[iz][5] = -(qt * rt);
        gxr[iz][6] = -(qt * rxzt + rt * qxzt);
        gxr[iz][7] = -(qxzt * rxzt);
    }

    //! Upper bounadry condition
    for (ix = 0; ix < nxb; ix++) {
        v = vp[ix][0] * dt / ds;
        qxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
        qt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
        qxzt = b / (b - 1.0f);

        v = vs[ix][0] * dt / ds;
        rxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
        rt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
        rxzt = b / (b - 1.0f);

        gzu[ix][0] = -(qxz + rxz);
        gzu[ix][1] = -(qxz * rxz);
        gzu[ix][2] = -(qt + rt);
        gzu[ix][3] = -(qxz * rt + qt * rxz + qxzt + rxzt);
        gzu[ix][4] = -(qxz * rxzt + rxz * qxzt);
        gzu[ix][5] = -(qt * rt);
        gzu[ix][6] = -(qt * rxzt + rt * qxzt);
        gzu[ix][7] = -(qxzt * rxzt);
    }

    //! Below bounadry condition
    for (ix = 0; ix < nxb; ix++) {
        v = vp[ix][nzb - 1] * dt / ds;
        qxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
        qt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
        qxzt = b / (b - 1.0f);

        v = vs[ix][nzb - 1] * dt / ds;
        rxz = (b * (beta + v) - v) / (beta + v) / (1.0f - b);
        rt = (b * (beta + v) - beta) / (beta + v) / (1.0f - b);
        rxzt = b / (b - 1.0f);

        gzb[ix][0] = -(qxz + rxz);
        gzb[ix][1] = -(qxz * rxz);
        gzb[ix][2] = -(qt + rt);
        gzb[ix][3] = -(qxz * rt + qt * rxz + qxzt + rxzt);
        gzb[ix][4] = -(qxz * rxzt + rxz * qxzt);
        gzb[ix][5] = -(qt * rt);
        gzb[ix][6] = -(qt * rxzt + rt * qxzt);
        gzb[ix][7] = -(qxzt * rxzt);
    }

    return 1;
}

//! Apply two order time hybrid abc for 2 demension simulation
int sjapplythabc2d(float **fp, float **cp, float **pp, float **gxl, float **gxr, float **gzu, float **gzb,
                   int nxb, int nzb, int nb, int marg) {
    int ii, ix, iz;
    float temp;

    for (ii = 0; ii <= nb; ii++) {
        //! Left & right bounadry condition
        for (iz = marg + nb - ii + 1; iz < nzb - marg - nb + ii - 1; iz++) {
            ix = marg + nb - ii;
            temp = gxl[iz][0] * fp[ix + 1][iz] + gxl[iz][1] * fp[ix + 2][iz]
                   + gxl[iz][2] * cp[ix + 0][iz] + gxl[iz][3] * cp[ix + 1][iz] + gxl[iz][4] * cp[ix + 2][iz]
                   + gxl[iz][5] * pp[ix + 0][iz] + gxl[iz][6] * pp[ix + 1][iz] + gxl[iz][7] * pp[ix + 2][iz];
            fp[ix][iz] = fp[ix][iz] * ((float) (nb - 1 - ii)) / ((float) nb) + temp * ((float) (1 + ii)) / ((float) nb);

            ix = nxb - marg - nb + ii - 1;
            temp = gxr[iz][0] * fp[ix - 1][iz] + gxr[iz][1] * fp[ix - 2][iz]
                   + gxr[iz][2] * cp[ix + 0][iz] + gxr[iz][3] * cp[ix - 1][iz] + gxr[iz][4] * cp[ix - 2][iz]
                   + gxr[iz][5] * pp[ix + 0][iz] + gxr[iz][6] * pp[ix - 1][iz] + gxr[iz][7] * pp[ix - 2][iz];
            fp[ix][iz] = fp[ix][iz] * ((float) (nb - 1 - ii)) / ((float) nb) + temp * ((float) (1 + ii)) / ((float) nb);
        }

        //! Upper & below bounadry condition
        for (ix = marg + nb - ii + 1; ix < nxb - marg - nb + ii - 1; ix++) {
            iz = marg + nb - ii;
            temp = gzu[ix][0] * fp[ix][iz + 1] + gzu[ix][1] * fp[ix][iz + 2]
                   + gzu[ix][2] * cp[ix][iz + 0] + gzu[ix][3] * cp[ix][iz + 1] + gzu[ix][4] * cp[ix][iz + 2]
                   + gzu[ix][5] * pp[ix][iz + 0] + gzu[ix][6] * pp[ix][iz + 1] + gzu[ix][7] * pp[ix][iz + 2];
            fp[ix][iz] = fp[ix][iz] * ((float) (nb - 1 - ii)) / ((float) nb) + temp * ((float) (1 + ii)) / ((float) nb);

            iz = nzb - marg - nb + ii - 1;
            temp = gzb[ix][0] * fp[ix][iz - 1] + gzb[ix][1] * fp[ix][iz - 2]
                   + gzb[ix][2] * cp[ix][iz - 0] + gzb[ix][3] * cp[ix][iz - 1] + gzb[ix][4] * cp[ix][iz - 2]
                   + gzb[ix][5] * pp[ix][iz - 0] + gzb[ix][6] * pp[ix][iz - 1] + gzb[ix][7] * pp[ix][iz - 2];
            fp[ix][iz] = fp[ix][iz] * ((float) (nb - 1 - ii)) / ((float) nb) + temp * ((float) (1 + ii)) / ((float) nb);
        }

        //! Corner
        ix = marg + nb - ii;
        iz = marg + nb - ii;
        temp = 0.5f * (gxl[iz][0] * fp[ix + 1][iz] + gxl[iz][1] * fp[ix + 2][iz]
                       + gxl[iz][2] * cp[ix + 0][iz] + gxl[iz][3] * cp[ix + 1][iz] + gxl[iz][4] * cp[ix + 2][iz]
                       + gxl[iz][5] * pp[ix + 0][iz] + gxl[iz][6] * pp[ix + 1][iz] + gxl[iz][7] * pp[ix + 2][iz]) +
               0.5f * (gzu[ix][0] * fp[ix][iz + 1] + gzu[ix][1] * fp[ix][iz + 2]
                       + gzu[ix][2] * cp[ix][iz + 0] + gzu[ix][3] * cp[ix][iz + 1] + gzu[ix][4] * cp[ix][iz + 2]
                       + gzu[ix][5] * pp[ix][iz + 0] + gzu[ix][6] * pp[ix][iz + 1] + gzu[ix][7] * pp[ix][iz + 2]);
        fp[ix][iz] = fp[ix][iz] * ((float) (nb - 1 - ii)) / ((float) nb) + temp * ((float) (1 + ii)) / ((float) nb);

        ix = nxb - marg - nb + ii - 1;
        iz = marg + nb - ii;
        temp = 0.5f * (gxr[iz][0] * fp[ix - 1][iz] + gxr[iz][1] * fp[ix - 2][iz]
                       + gxr[iz][2] * cp[ix - 0][iz] + gxr[iz][3] * cp[ix - 1][iz] + gxr[iz][4] * cp[ix - 2][iz]
                       + gxr[iz][5] * pp[ix - 0][iz] + gxr[iz][6] * pp[ix - 1][iz] + gxr[iz][7] * pp[ix - 2][iz]) +
               0.5f * (gzu[ix][0] * fp[ix][iz + 1] + gzu[ix][1] * fp[ix][iz + 2]
                       + gzu[ix][2] * cp[ix][iz + 0] + gzu[ix][3] * cp[ix][iz + 1] + gzu[ix][4] * cp[ix][iz + 2]
                       + gzu[ix][5] * pp[ix][iz + 0] + gzu[ix][6] * pp[ix][iz + 1] + gzu[ix][7] * pp[ix][iz + 2]);
        fp[ix][iz] = fp[ix][iz] * ((float) (nb - 1 - ii)) / ((float) nb) + temp * ((float) (1 + ii)) / ((float) nb);

        ix = marg + nb - ii;
        iz = nzb - marg - nb + ii - 1;
        temp = 0.5f * (gxl[iz][0] * fp[ix + 1][iz] + gxl[iz][1] * fp[ix + 2][iz]
                       + gxl[iz][2] * cp[ix + 0][iz] + gxl[iz][3] * cp[ix + 1][iz] + gxl[iz][4] * cp[ix + 2][iz]
                       + gxl[iz][5] * pp[ix + 0][iz] + gxl[iz][6] * pp[ix + 1][iz] + gxl[iz][7] * pp[ix + 2][iz]) +
               0.5f * (gzb[ix][0] * fp[ix][iz - 1] + gzb[ix][1] * fp[ix][iz - 2]
                       + gzb[ix][2] * cp[ix][iz - 0] + gzb[ix][3] * cp[ix][iz - 1] + gzb[ix][4] * cp[ix][iz - 2]
                       + gzb[ix][5] * pp[ix][iz - 0] + gzb[ix][6] * pp[ix][iz - 1] + gzb[ix][7] * pp[ix][iz - 2]);
        fp[ix][iz] = fp[ix][iz] * ((float) (nb - 1 - ii)) / ((float) nb) + temp * ((float) (1 + ii)) / ((float) nb);

        ix = nxb - marg - nb + ii - 1;
        iz = nzb - marg - nb + ii - 1;
        temp = 0.5f * (gxr[iz][0] * fp[ix - 1][iz] + gxr[iz][1] * fp[ix - 2][iz]
                       + gxr[iz][2] * cp[ix - 0][iz] + gxr[iz][3] * cp[ix - 1][iz] + gxr[iz][4] * cp[ix - 2][iz]
                       + gxr[iz][5] * pp[ix - 0][iz] + gxr[iz][6] * pp[ix - 1][iz] + gxr[iz][7] * pp[ix - 2][iz]) +
               0.5f * (gzb[ix][0] * fp[ix][iz - 1] + gzb[ix][1] * fp[ix][iz - 2]
                       + gzb[ix][2] * cp[ix][iz - 0] + gzb[ix][3] * cp[ix][iz - 1] + gzb[ix][4] * cp[ix][iz - 2]
                       + gzb[ix][5] * pp[ix][iz - 0] + gzb[ix][6] * pp[ix][iz - 1] + gzb[ix][7] * pp[ix][iz - 2]);
        fp[ix][iz] = fp[ix][iz] * ((float) (nb - 1 - ii)) / ((float) nb) + temp * ((float) (1 + ii)) / ((float) nb);
    }

    return 1;
}
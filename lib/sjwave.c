// Author: Hou, Sian - sianhou1987@outlook.com

#include "sjwave.h"
#include "sjinc.h"

int sjricker1d(float *ricker, int nt, int t0, float dt, float fp, float amp) {

    int it = 0;

    float tmpf = (float) (pi * pi * fp * fp), tmpt = dt * dt;

    double tmpr;

    for (it = 0; it < nt; ++it) {
        tmpr = (1.0 - 2.0 * tmpf * (it - t0) * (it - t0) * tmpt) * exp(-tmpf * (it - t0) * (it - t0) * tmpt);
        ricker[it] = ((float) tmpr) * amp;
    }

    return 1;
}

void sjextend2d(float **z, int nx, int nz, int ex0, int ex1, int ez0, int ez1, float **x) {
    int ix, iz;

    //! Centering
    for (ix = 0; ix < nx; ++ix)
        for (iz = 0; iz < nz; ++iz)
            z[ix + ex0][iz + ez0] = x[ix][iz];

    //! Left
    for (ix = 0; ix < ex0; ++ix)
        for (iz = 0; iz < nz; ++iz)
            z[ix][iz + ez0] = x[0][iz];

    //! Right
    for (ix = 0; ix < ex1; ++ix)
        for (iz = 0; iz < nz; ++iz)
            z[ex0 + nx + ix][iz + ez0] = x[nx - 1][iz];

    //! Upper
    for (ix = 0; ix < ex0 + nx + ex1; ++ix)
        for (iz = 0; iz < ez0; ++iz)
            z[ix][iz] = z[ix][ez0];

    //! Below
    for (ix = 0; ix < ex0 + nx + ex1; ++ix)
        for (iz = 0; iz < ez1; ++iz)
            z[ix][ez0 + nz + iz] = z[ix][ez0 + nz - 1];
}

void sjextract2d(float **z, int x0, int z0, int nx, int nz, float **x) {
    int ix, iz;
    for (ix = 0; ix < nx; ++ix)
        for (iz = 0; iz < nz; ++iz)
            z[ix][iz] = x[x0 + ix][z0 + iz];
}

void sjfilter2dx(float **a, int n2, int n1, char *mode) {
    int ix, iz;
    if (strcmp(mode, "laplace") == 0) {
        float **p = (float **) sjalloc2d(n2, n1, sizeof(float));
        memcpy(p[0], a[0], n2 * n1 * sizeof(float));
        for (ix = 2; ix < n2 - 2; ++ix)
            for (iz = 2; iz < n1 - 2; ++iz)
                a[ix][iz] = -4.0f * p[ix][iz] + p[ix - 1][iz] + p[ix + 1][iz] + p[ix][iz - 1] + p[ix][iz + 1];

        sjmfree2d(p);
    }
}

void sjsetsurface(float **a, int n2, int n1, float val) {
    int ix, iz;
    for (ix = 0; ix < n2; ++ix)
        for (iz = 0; iz < n1; ++iz)
            a[ix][iz] = val;
}

/**********************************************************************************************/
/* ! Finite Difference                                                                        */
/**********************************************************************************************/

#define C50 ( 1.239407e+0f)
#define C51 (-1.105315e-1f)
#define C52 ( 2.496329e-2f)
#define C53 (-5.804879e-3f)
#define C54 ( 9.358680e-4f)

#define sjmsgfd2dn1(p, ix, iz) (C50*(p[ix][iz+1]-p[ix][iz-0]) + \
                                C51*(p[ix][iz+2]-p[ix][iz-1]) + \
                                C52*(p[ix][iz+3]-p[ix][iz-2]) + \
                                C53*(p[ix][iz+4]-p[ix][iz-3]) + \
                                C54*(p[ix][iz+5]-p[ix][iz-4]) )

#define sjmsgfd2dn2(p, ix, iz) (C50*(p[ix+1][iz]-p[ix-0][iz]) + \
                                C51*(p[ix+2][iz]-p[ix-1][iz]) + \
                                C52*(p[ix+3][iz]-p[ix-2][iz]) + \
                                C53*(p[ix+4][iz]-p[ix-3][iz]) + \
                                C54*(p[ix+5][iz]-p[ix-4][iz]) )

#define B60 (-2.982778e+0f)
#define B61 ( 1.714286e+0f)
#define B62 (-2.678571e-1f)
#define B63 ( 5.291005e-2f)
#define B64 (-8.928571e-3f)
#define B65 ( 1.038961e-3f)
#define B66 (-6.012506e-5f)

#define B611  0.562500000000f
#define B612 -0.112500000000f
#define B613  0.012500000000f
#define B622  0.022500000000f
#define B623 -0.002500000000f
#define B633  0.000277777778f

#define sjmfd2dn1(a, ix, iz)( B60* a[ix][iz]+ \
                            B61*(a[ix][iz+1]+a[ix][iz-1]) + \
                            B62*(a[ix][iz+2]+a[ix][iz-2]) + \
                            B63*(a[ix][iz+3]+a[ix][iz-3]) + \
                            B64*(a[ix][iz+4]+a[ix][iz-4]) + \
                            B65*(a[ix][iz+5]+a[ix][iz-5]) + \
                            B66*(a[ix][iz+6]+a[ix][iz-6]) )

#define sjmfd2dn2(a, ix, iz)( B60* a[ix][iz]+ \
                            B61*(a[ix+1][iz]+a[ix-1][iz]) + \
                            B62*(a[ix+2][iz]+a[ix-2][iz]) + \
                            B63*(a[ix+3][iz]+a[ix-3][iz]) + \
                            B64*(a[ix+4][iz]+a[ix-4][iz]) + \
                            B65*(a[ix+5][iz]+a[ix-5][iz]) + \
                            B66*(a[ix+6][iz]+a[ix-6][iz]) )

#define sjmfd2dnc(a, ix, iz)( B611*(a[ix+1][iz+1]-a[ix-1][iz+1]-a[ix+1][iz-1]+a[ix-1][iz-1]) + \
                            B612*(a[ix+1][iz+2]-a[ix-1][iz+2]-a[ix+1][iz-2]+a[ix-1][iz-2]+a[ix+2][iz+1]-a[ix-2][iz+1]-a[ix+2][iz-1]+a[ix-2][iz-1]) + \
                            B613*(a[ix+1][iz+3]-a[ix-1][iz+3]-a[ix+1][iz-3]+a[ix-1][iz-3]+a[ix+3][iz+1]-a[ix-3][iz+1]-a[ix+3][iz-1]+a[ix-3][iz-1]) + \
                            B622*(a[ix+2][iz+2]-a[ix-2][iz+2]-a[ix+2][iz-2]+a[ix-2][iz-2]) + \
                            B623*(a[ix+2][iz+3]-a[ix-2][iz+3]-a[ix+2][iz-3]+a[ix-2][iz-3]+a[ix+3][iz+2]-a[ix-3][iz+2]-a[ix+3][iz-2]+a[ix-3][iz-2]) + \
                            B633*(a[ix+3][iz+3]-a[ix-3][iz+3]-a[ix+3][iz-3]+a[ix-3][iz-3]))

/**********************************************************************************************/
/* ! Acousitc                                                                                 */
/**********************************************************************************************/

//! Two dimension constant density acoustic forward exploration
void sjafor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int it, ir, ix, iz;
    const int marg = 6;

    //------------------------ Option ------------------------//
    int nt = opt->nt;
    int k1 = opt->k1;
    int jsnap = opt->jsnap;
    int srcrange = opt->srcrange;
    int srctrunc = opt->srctrunc;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    float fp = opt->fp;
    float amp = opt->amp;
    float srcdecay = opt->srcdecay;
    int nb = opt->nb;
    float ds = opt->ds;
    int ycutdirect = opt->ycutdirect;
    float *wavelet = sjmflloc1d(nt);
    sjricker1d(wavelet, nt, k1, dt, fp, amp);

    //------------------------ Survey ------------------------//
    int nx = sur->nx;
    int nz = sur->nz;
    int sx = sur->sx + nb + marg;
    int sz = sur->sz + nb + marg;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = dt2 / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    sjveczerof(p2[0], nxb * nzb);
    sjveczerof(p1[0], nxb * nzb);
    sjveczerof(p0[0], nxb * nzb);

    //------------------------ Wavefield exploration ------------------------//
    for (it = 0; it < nt; it++) {

        //! Source
        for (ix = -srcrange; ix <= srcrange; ix++)
            for (iz = -srcrange; iz <= srcrange; iz++)
                p1[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));

        //! Calculate veloctiy
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz))
                             + 2.0f * p1[ix][iz] - p0[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
        for (ir = 0; ir < nr; ir++)
            wav->profz[ir][it] = p1[nb + marg + rx[ir]][nb + marg + rz[ir]];

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++)
                    wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] = p1[ix][iz];

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
    }

    //------------------------ Cut direct wav ------------------------//
    if (ycutdirect == 1) {
        //------------------------ Model ------------------------//
        for (ix = 0; ix < nxb; ix++)
            for (iz = 0; iz < nzb; iz++)
                cp[ix][iz] = geo->vp2d[sx - marg - nb][sz - marg - nb];

        //------------------------ Initialization ------------------------//
        //! Boundary condition
        sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);

        //! Model
        for (ix = 0; ix < nxb; ix++)
            for (iz = 0; iz < nzb; iz++)
                cp[ix][iz] = cp[ix][iz] * cp[ix][iz] * ids2;

        //! Wavefield
        memset(p2[0], 0, nxb * nzb * sizeof(float));
        memset(p1[0], 0, nxb * nzb * sizeof(float));
        memset(p0[0], 0, nxb * nzb * sizeof(float));

        //! Wavefield exploration
        for (it = 0; it < nt; it++) {

            //! Source
            if (it < srctrunc)
                for (ix = -srcrange; ix <= srcrange; ix++)
                    for (iz = -srcrange; iz <= srcrange; iz++)
                        p1[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));

            //! Calculate veloctiy
            for (ix = marg; ix < nxb - marg; ix++)
                for (iz = marg; iz < nzb - marg; iz++)
                    p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz))
                                 + 2.0f * p1[ix][iz] - p0[ix][iz];

            //! Boundary condition
            sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

            //! Record
            for (ir = 0; ir < nr; ir++)
                wav->profz[ir][it] -= p1[nb + marg + rx[ir]][nb + marg + rz[ir]];

            //! Update
            memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
            memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
        }
    }

    sjmfree1d(wavelet);

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p2);
    sjmfree2d(p1);
    sjmfree2d(p0);
}

//! Two dimension constant density acoustic scatter forward exploration
void sjasfor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int it, ir, ix, iz;
    const int marg = 6;

    //------------------------ Option ------------------------//
    int nt = opt->nt;
    int k1 = opt->k1;
    int jsnap = opt->jsnap;
    int srcrange = opt->srcrange;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    float fp = opt->fp;
    float amp = opt->amp;
    float srcdecay = opt->srcdecay;
    int nb = opt->nb;
    float ds = opt->ds;
    float *wavelet = sjmflloc1d(nt);
    sjricker1d(wavelet, nt, k1, dt, fp, amp);

    //------------------------ Survey ------------------------//
    int nx = sur->nx;
    int nz = sur->nz;
    int sx = sur->sx + nb + marg;
    int sz = sur->sz + nb + marg;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = dt2 / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    float **ipp = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);
    sjextend2d(ipp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->izz2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);
    float **s2 = sjmflloc2d(nxb, nzb);
    float **s1 = sjmflloc2d(nxb, nzb);
    float **s0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    sjveczerof(p2[0], nxb * nzb);
    sjveczerof(p1[0], nxb * nzb);
    sjveczerof(p0[0], nxb * nzb);
    sjveczerof(s2[0], nxb * nzb);
    sjveczerof(s1[0], nxb * nzb);
    sjveczerof(s0[0], nxb * nzb);

    //------------------------ Wavefield exploration ------------------------//
    for (it = 0; it < nt; ++it) {

        //! Stack primary source
        for (ix = -srcrange; ix <= srcrange; ++ix)
            for (iz = -srcrange; iz <= srcrange; ++iz)
                p1[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));

        //! Calculate primary wavefield
        for (ix = marg; ix < nxb - marg; ++ix)
            for (iz = marg; iz < nzb - marg; ++iz)
                p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz))
                             + 2.0f * p1[ix][iz] - p0[ix][iz];

        //! Stack scatter source
        for (ix = marg + 10; ix < nxb - marg - 10; ++ix)
            for (iz = marg + 10; iz < nzb - marg - 10; ++iz)
                s1[ix][iz] += (p2[ix][iz] - 2.0f * p1[ix][iz] + p0[ix][iz]) / dt2 * ipp[ix][iz];

        //! Calculate scatter wavefield
        for (ix = marg; ix < nxb - marg; ++ix)
            for (iz = marg; iz < nzb - marg; ++iz)
                s2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(s1, ix, iz) + sjmfd2dn2(s1, ix, iz))
                             + 2.0f * s1[ix][iz] - s0[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplythabc2d(s2, s1, s0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
        for (ir = 0; ir < nr; ++ir)
            wav->profz[ir][it] = s1[nb + marg + rx[ir]][nb + marg + rz[ir]];

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ++ix)
                for (iz = nb + marg; iz < nzb - nb - marg; ++iz) {
                    wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] = p1[ix][iz];
                    wav->fsz2d[it / jsnap][ix - nb - marg][iz - nb - marg] = s1[ix][iz];
                }

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
        memcpy(s0[0], s1[0], nxb * nzb * sizeof(float));
        memcpy(s1[0], s2[0], nxb * nzb * sizeof(float));
    }

    sjmfree1d(wavelet);

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p2);
    sjmfree2d(p1);
    sjmfree2d(p0);

    sjmfree2d(s2);
    sjmfree2d(s1);
    sjmfree2d(s0);
}

//! Two dimension constant density acoustic RTM backward exploration
void sjartmbac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //! Runtime
    int it, ir, ix, iz;
    const int marg = 6;

    //! Option
    int nt = opt->nt;
    int jsnap = opt->jsnap;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    int nb = opt->nb;
    float ds = opt->ds;

    //! Survey
    int nx = sur->nx;
    int nz = sur->nz;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = dt2 / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    sjveczerof(p2[0], nxb * nzb);
    sjveczerof(p1[0], nxb * nzb);
    sjveczerof(p0[0], nxb * nzb);

    //------------------------ Wavefield exploration ------------------------//
    for (it = nt - 1; it >= 0; --it) {
        //! Source#
        if (opt->ystacksrc == 1) {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] += wav->profz[ir][it];
        } else {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] = wav->profz[ir][it];
        }

        //! Calculate velocity
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz))
                             + 2.0f * p1[ix][iz] - p0[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Wavefield
        if ((it % jsnap) == 0) {
            for (ix = nb + marg; ix < nxb - nb - marg; ix++) {
                for (iz = nb + marg; iz < nzb - nb - marg; iz++) {
                    geo->izz2d[ix - nb - marg][iz - nb - marg] +=
                            wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] * p1[ix][iz];
                    geo->nzz2d[ix - nb - marg][iz - nb - marg] +=
                            wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] *
                            wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg];
                }
            }
        }

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
    }

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p0);
    sjmfree2d(p1);
    sjmfree2d(p2);
}

//! Two dimension constant density acoustic Time-Shift RTM backward exploration
void sjatsrtmbac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //! Runtime
    int it, ir, ix, iz, shift;
    const int marg = 6;

    //! Option
    int nt = opt->nt;
    int jsnap = opt->jsnap;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    int nb = opt->nb;
    float ds = opt->ds;

    //! Survey
    int nx = sur->nx;
    int nz = sur->nz;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = dt2 / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    sjveczerof(p2[0], nxb * nzb);
    sjveczerof(p1[0], nxb * nzb);
    sjveczerof(p0[0], nxb * nzb);

    //------------------------ Wavefield exploration ------------------------//
    for (it = nt - 1; it >= 0; --it) {
        //! Source#
        if (opt->ystacksrc == 1) {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] += wav->profz[ir][it];
        } else {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] = wav->profz[ir][it];
        }

        //! Calculate velocity
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz))
                             + 2.0f * p1[ix][iz] - p0[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Time-Shift imaging
        if (((it % jsnap) == 0) && ((it / jsnap) >= opt->maxshift) && ((it / jsnap) < (opt->nsnap - opt->maxshift)))
            for (shift = -opt->maxshift; shift <= opt->maxshift; ++shift)
                for (ix = nb + marg; ix < nxb - nb - marg; ++ix)
                    for (iz = nb + marg; iz < nzb - nb - marg; ++iz)
                        geo->izz3d[opt->maxshift + shift][ix - nb - marg][iz - nb - marg] +=
                                wav->fwz2d[it / jsnap + shift][ix - nb - marg][iz - nb - marg] * p1[ix][iz];

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
    }

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p0);
    sjmfree2d(p1);
    sjmfree2d(p2);
}

//! Two dimension constant density acoustic FWI backward exploration
void sjafwibac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //! Runtime
    int it, ir, ix, iz;
    const int marg = 6;

    //! Option
    int nt = opt->nt;
    int jsnap = opt->jsnap;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    int nb = opt->nb;
    float ds = opt->ds;

    //! Survey
    int nx = sur->nx;
    int nz = sur->nz;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = 1.0f / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    sjveczerof(p2[0], nxb * nzb);
    sjveczerof(p1[0], nxb * nzb);
    sjveczerof(p0[0], nxb * nzb);

    //------------------------ Wavefield exploration ------------------------//
    for (it = nt - 1; it >= 0; --it) {

        //! Source
        if (opt->ystacksrc == 1) {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] += wav->profz[ir][it];
        } else {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] = wav->profz[ir][it];
        }

        //! Calculate velocity
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p0[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz)) * dt2
                             + 2.0f * p1[ix][iz] - p2[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p0, p1, p2, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++) {
                    geo->gzz2d[ix - nb - marg][iz - nb - marg] +=
                            wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] * p1[ix][iz];
                    geo->nzz2d[ix - nb - marg][iz - nb - marg] +=
                            wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] *
                            wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg];
                }

        //! Update
        memcpy(p2[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p0[0], nxb * nzb * sizeof(float));
    }

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p0);
    sjmfree2d(p1);
    sjmfree2d(p2);
}

//! Two dimension constant density acoustic LSRTM forward simulation
void sjalsrtmfor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int it, ir, ix, iz;
    const int marg = 6;

    //------------------------ Option ------------------------//
    int nt = opt->nt;
    int k1 = opt->k1;
    int jsnap = opt->jsnap;
    int srcrange = opt->srcrange;
    int srctrunc = opt->srctrunc;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    float fp = opt->fp;
    float amp = opt->amp;
    float srcdecay = opt->srcdecay;
    int nb = opt->nb;
    float ds = opt->ds;
    float *wavelet = sjmflloc1d(nt);
    sjricker1d(wavelet, nt, k1, dt, fp, amp);

    //------------------------ Survey ------------------------//
    int nx = sur->nx;
    int nz = sur->nz;
    int sx = sur->sx + nb + marg;
    int sz = sur->sz + nb + marg;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = 1.0f / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    float **ipp = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);
    sjextend2d(ipp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->izz2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);
    float **s2 = sjmflloc2d(nxb, nzb);
    float **s1 = sjmflloc2d(nxb, nzb);
    float **s0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    memset(p2[0], 0, nxb * nzb * sizeof(float));
    memset(p1[0], 0, nxb * nzb * sizeof(float));
    memset(p0[0], 0, nxb * nzb * sizeof(float));
    memset(s2[0], 0, nxb * nzb * sizeof(float));
    memset(s1[0], 0, nxb * nzb * sizeof(float));
    memset(s0[0], 0, nxb * nzb * sizeof(float));

    //------------------------ Wavefield exploration ------------------------//
    for (it = 0; it < nt; it++) {
        //! Source
        if (it < srctrunc)
            for (ix = -srcrange; ix <= srcrange; ix++)
                for (iz = -srcrange; iz <= srcrange; iz++)
                    p1[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));

        //! Laplace operator
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz)) * dt2
                             + 2.0f * p1[ix][iz] - p0[ix][iz];

        //! Scatter source
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                s1[ix][iz] += (p2[ix][iz] - 2.0f * p1[ix][iz] + p0[ix][iz]) / dt2 * ipp[ix][iz];

        //! Scatter wavefield
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                s2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(s1, ix, iz) + sjmfd2dn2(s1, ix, iz)) * dt2
                             + 2.0f * s1[ix][iz] - s0[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplythabc2d(s2, s1, s0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
        for (ir = 0; ir < nr; ir++)
            wav->profz[ir][it] = s1[nb + marg + rx[ir]][nb + marg + rz[ir]];

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++)
                    wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] =
                            (p2[ix][iz] - 2.0f * p1[ix][iz] + p0[ix][iz]) / dt2;

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));

        memcpy(s0[0], s1[0], nxb * nzb * sizeof(float));
        memcpy(s1[0], s2[0], nxb * nzb * sizeof(float));
    }

    sjmfree1d(wavelet);

    sjmfree2d(cp);
    sjmfree2d(ipp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p2);
    sjmfree2d(p1);
    sjmfree2d(p0);
    sjmfree2d(s2);
    sjmfree2d(s1);
    sjmfree2d(s0);
}

//! Two dimension constant density acoustic LSRTM backward exploration
void sjalsrtmbac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //! Runtime
    int it, ir, ix, iz;
    const int marg = 6;

    //! Option
    int nt = opt->nt;
    int jsnap = opt->jsnap;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    int nb = opt->nb;
    float ds = opt->ds;

    //! Survey
    int nx = sur->nx;
    int nz = sur->nz;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = dt2 / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    sjveczerof(p2[0], nxb * nzb);
    sjveczerof(p1[0], nxb * nzb);
    sjveczerof(p0[0], nxb * nzb);

    //------------------------ Wavefield exploration ------------------------//
    for (it = nt - 1; it >= 0; --it) {
        //! Source#
        if (opt->ystacksrc == 1) {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] += wav->profz[ir][it];
        } else {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] = wav->profz[ir][it];
        }

        //! Calculate velocity
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++)
                p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz))
                             + 2.0f * p1[ix][iz] - p0[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Wavefield
        if ((it % jsnap) == 0) {
            for (ix = nb + marg; ix < nxb - nb - marg; ix++) {
                for (iz = nb + marg; iz < nzb - nb - marg; iz++) {
                    geo->gzz2d[ix - nb - marg][iz - nb - marg] +=
                            wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] * p1[ix][iz];
                    geo->nzz2d[ix - nb - marg][iz - nb - marg] +=
                            wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] *
                            wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg];
                }
            }
        }

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
    }

    sjmfree2d(cp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p0);
    sjmfree2d(p1);
    sjmfree2d(p2);
}

//! Two dimension constant density acoustic RTI backward exploration
void sjartibac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int it, ir, ix, iz, shift;
    const int marg = 6;

    //------------------------ Option ------------------------//
    int nt = opt->nt;
    int jsnap = opt->jsnap;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    int nb = opt->nb;
    float ds = opt->ds;

    //------------------------ Survey ------------------------//
    int nx = sur->nx;
    int nz = sur->nz;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = dt2 / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    float **ipp = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);
    sjextend2d(ipp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->izz2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);
    float **s2 = sjmflloc2d(nxb, nzb);
    float **s1 = sjmflloc2d(nxb, nzb);
    float **s0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    memset(p2[0], 0, nxb * nzb * sizeof(float));
    memset(p1[0], 0, nxb * nzb * sizeof(float));
    memset(p0[0], 0, nxb * nzb * sizeof(float));
    memset(s2[0], 0, nxb * nzb * sizeof(float));
    memset(s1[0], 0, nxb * nzb * sizeof(float));
    memset(s0[0], 0, nxb * nzb * sizeof(float));

    //------------------------ Wavefield exploration ------------------------//
    for (it = nt - 1; it >= 0; --it) {

        //! Source
        if (opt->ystacksrc == 1) {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] += wav->profz[ir][it];
        } else {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] = wav->profz[ir][it];
        }

        //! Primary wavefield
        for (ix = marg; ix < nxb - marg; ++ix)
            for (iz = marg; iz < nzb - marg; ++iz)
                p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz))
                             + 2.0f * p1[ix][iz] - p0[ix][iz];

        //! Scatter source
        for (ix = nb + marg + 10; ix < nxb - nb - marg - 10; ++ix)
            for (iz = nb + marg + 10; iz < nzb - nb - marg - 10; ++iz)
                s1[ix][iz] += (p2[ix][iz] - 2.0f * p1[ix][iz] + p0[ix][iz]) / dt2 * ipp[ix][iz];

        //! Scatter wavefield
        for (ix = marg; ix < nxb - marg; ++ix)
            for (iz = marg; iz < nzb - marg; ++iz)
                s2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(s1, ix, iz) + sjmfd2dn2(s1, ix, iz))
                             + 2.0f * s1[ix][iz] - s0[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplythabc2d(s2, s1, s0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Imaging
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ++ix)
                for (iz = nb + marg; iz < nzb - nb - marg; ++iz) {
                    geo->gzz2d[ix - nb - marg][iz - nb - marg] +=
                            wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] * s1[ix][iz] +
                            wav->fsz2d[it / jsnap][ix - nb - marg][iz - nb - marg] * p1[ix][iz];
                }

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
        memcpy(s0[0], s1[0], nxb * nzb * sizeof(float));
        memcpy(s1[0], s2[0], nxb * nzb * sizeof(float));
    }

    sjmfree2d(cp);
    sjmfree2d(ipp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p2);
    sjmfree2d(p1);
    sjmfree2d(p0);
    sjmfree2d(s2);
    sjmfree2d(s1);
    sjmfree2d(s0);
}

//! Two dimension constant density acoustic WTI backward exploration
void sjawtibac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int it, ir, ix, iz, shift;
    const int marg = 6;

    //------------------------ Option ------------------------//
    int nt = opt->nt;
    int jsnap = opt->jsnap;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    int nb = opt->nb;
    float ds = opt->ds;

    //------------------------ Survey ------------------------//
    int nx = sur->nx;
    int nz = sur->nz;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = dt2 / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    float **ipp = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);
    sjextend2d(ipp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->izz2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **p2 = sjmflloc2d(nxb, nzb);
    float **p1 = sjmflloc2d(nxb, nzb);
    float **p0 = sjmflloc2d(nxb, nzb);
    float **s2 = sjmflloc2d(nxb, nzb);
    float **s1 = sjmflloc2d(nxb, nzb);
    float **s0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cp, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    //! Wavefield
    memset(p2[0], 0, nxb * nzb * sizeof(float));
    memset(p1[0], 0, nxb * nzb * sizeof(float));
    memset(p0[0], 0, nxb * nzb * sizeof(float));
    memset(s2[0], 0, nxb * nzb * sizeof(float));
    memset(s1[0], 0, nxb * nzb * sizeof(float));
    memset(s0[0], 0, nxb * nzb * sizeof(float));

    //------------------------ Wavefield exploration ------------------------//
    for (it = nt - 1; it >= 0; --it) {

        //! Source
        if (opt->ystacksrc == 1) {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] += wav->profz[ir][it];
        } else {
            for (ir = 0; ir < nr; ir++)
                p1[nb + marg + rx[ir]][nb + marg + rz[ir]] = wav->profz[ir][it];
        }

        //! Calculate primary wavefield
        for (ix = marg; ix < nxb - marg; ++ix)
            for (iz = marg; iz < nzb - marg; ++iz)
                p2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(p1, ix, iz) + sjmfd2dn2(p1, ix, iz))
                             + 2.0f * p1[ix][iz] - p0[ix][iz];

        //! Stack scatter source
        for (ix = nb + marg + 10; ix < nxb - nb - marg - 10; ++ix)
            for (iz = nb + marg + 10; iz < nzb - nb - marg - 10; ++iz)
                s1[ix][iz] += (p2[ix][iz] - 2.0f * p1[ix][iz] + p0[ix][iz]) / dt2 * ipp[ix][iz];

        //! Calculate scatter wavefield
        for (ix = marg; ix < nxb - marg; ++ix)
            for (iz = marg; iz < nzb - marg; ++iz)
                s2[ix][iz] = cp[ix][iz] * (sjmfd2dn1(s1, ix, iz) + sjmfd2dn2(s1, ix, iz))
                             + 2.0f * s1[ix][iz] - s0[ix][iz];

        //! Boundary condition
        sjapplythabc2d(p2, p1, p0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplythabc2d(s2, s1, s0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Time-Shift imaging
        if (((it % jsnap) == 0) && ((it / jsnap) >= opt->maxshift) && ((it / jsnap) < (opt->nsnap - opt->maxshift)))
            for (shift = -opt->maxshift; shift <= opt->maxshift; ++shift)
                for (ix = nb + marg; ix < nxb - nb - marg; ++ix)
                    for (iz = nb + marg; iz < nzb - nb - marg; ++iz)
                        geo->gzz3d[opt->maxshift + shift][ix - nb - marg][iz - nb - marg] +=
                                wav->fwz2d[it / jsnap + shift][ix - nb - marg][iz - nb - marg] * s1[ix][iz] +
                                wav->fsz2d[it / jsnap + shift][ix - nb - marg][iz - nb - marg] * p1[ix][iz];

        //! Update
        memcpy(p0[0], p1[0], nxb * nzb * sizeof(float));
        memcpy(p1[0], p2[0], nxb * nzb * sizeof(float));
        memcpy(s0[0], s1[0], nxb * nzb * sizeof(float));
        memcpy(s1[0], s2[0], nxb * nzb * sizeof(float));
    }

    sjmfree2d(cp);
    sjmfree2d(ipp);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(p2);
    sjmfree2d(p1);
    sjmfree2d(p0);
    sjmfree2d(s2);
    sjmfree2d(s1);
    sjmfree2d(s0);
}

/**********************************************************************************************/
/* ! Elastic                                                                                  */
/**********************************************************************************************/

//! Two dimension constant density elastic forward exploration
void sjefor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int it, ir, ix, iz;
    const int marg = 6;

    //------------------------ Option ------------------------//
    int nt = opt->nt;
    int k1 = opt->k1;
    int jsnap = opt->jsnap;
    int srcrange = opt->srcrange;
    int srctrunc = opt->srctrunc;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    float fp = opt->fp;
    float amp = opt->amp;
    float srcdecay = opt->srcdecay;
    int nb = opt->nb;
    float ds = opt->ds;
    int ycutdirect = opt->ycutdirect;
    float *wavelet = sjmflloc1d(nt);
    sjricker1d(wavelet, nt, k1, dt, fp, amp);

    //------------------------ Survey ------------------------//
    int nx = sur->nx;
    int nz = sur->nz;
    int sx = sur->sx + nb + marg;
    int sz = sur->sz + nb + marg;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = dt2 / ds / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    float **cs = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);
    sjextend2d(cs, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vs2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **u2 = sjmflloc2d(nxb, nzb);
    float **u1 = sjmflloc2d(nxb, nzb);
    float **u0 = sjmflloc2d(nxb, nzb);
    float **w2 = sjmflloc2d(nxb, nzb);
    float **w1 = sjmflloc2d(nxb, nzb);
    float **w0 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cs, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
    sjvecmulf(cs[0], nxb * nzb, ids2, cs[0], cs[0]);
    //! Wavefield
    sjveczerof(u2[0], nxb * nzb);
    sjveczerof(u1[0], nxb * nzb);
    sjveczerof(u0[0], nxb * nzb);
    sjveczerof(w2[0], nxb * nzb);
    sjveczerof(w1[0], nxb * nzb);
    sjveczerof(w0[0], nxb * nzb);

    //------------------------ Wavefield exploration ------------------------//
    for (it = 0; it < nt; it++) {

        //! Source
        for (ix = -srcrange; ix <= srcrange; ix++)
            for (iz = -srcrange; iz <= srcrange; iz++) {
                u1[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
                w1[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
            }

        //! Calculate veloctiy
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++) {
                u2[ix][iz] = cp[ix][iz] * sjmfd2dn2(u1, ix, iz) + cs[ix][iz] * sjmfd2dn1(u1, ix, iz) +
                             (cp[ix][iz] - cs[ix][iz]) * sjmfd2dnc(w1, ix, iz) + 2.0f * u1[ix][iz] - u0[ix][iz];

                w2[ix][iz] = cp[ix][iz] * sjmfd2dn1(w1, ix, iz) + cs[ix][iz] * sjmfd2dn2(w1, ix, iz) +
                             (cp[ix][iz] - cs[ix][iz]) * sjmfd2dnc(u1, ix, iz) + 2.0f * w1[ix][iz] - w0[ix][iz];
            }

        //! Boundary condition
        sjapplythabc2d(u2, u1, u0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplythabc2d(w2, w1, w0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
        for (ir = 0; ir < nr; ir++) {
            wav->profx[ir][it] = u1[nb + marg + rx[ir]][nb + marg + rz[ir]];
            wav->profz[ir][it] = w1[nb + marg + rx[ir]][nb + marg + rz[ir]];
        }

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++) {
                    wav->fwx2d[it / jsnap][ix - nb - marg][iz - nb - marg] = u1[ix][iz];
                    wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] = w1[ix][iz];
                }

        //! Update
        memcpy(u0[0], u1[0], nxb * nzb * sizeof(float));
        memcpy(u1[0], u2[0], nxb * nzb * sizeof(float));
        memcpy(w0[0], w1[0], nxb * nzb * sizeof(float));
        memcpy(w1[0], w2[0], nxb * nzb * sizeof(float));
    }

    //------------------------ Cut direct wav ------------------------//
    if (ycutdirect == 1) {
        //------------------------ Model ------------------------//
        for (ix = 0; ix < nxb; ix++)
            for (iz = 0; iz < nzb; iz++) {
                cp[ix][iz] = geo->vp2d[sx - marg - nb][sz - marg - nb];
                cs[ix][iz] = geo->vs2d[sx - marg - nb][sz - marg - nb];
            }

        //------------------------ Initialization ------------------------//
        //! Boundary condition
        sjinitthabc2d(cp, cs, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
        //! Model
        sjvecmulf(cp[0], nxb * nzb, ids2, cp[0], cp[0]);
        sjvecmulf(cs[0], nxb * nzb, ids2, cs[0], cs[0]);
        //! Wavefield
        sjveczerof(u2[0], nxb * nzb);
        sjveczerof(u1[0], nxb * nzb);
        sjveczerof(u0[0], nxb * nzb);
        sjveczerof(w2[0], nxb * nzb);
        sjveczerof(w1[0], nxb * nzb);
        sjveczerof(w0[0], nxb * nzb);

        //! Wavefield exploration
        for (it = 0; it < nt; it++) {

            //! Source
            for (ix = -srcrange; ix <= srcrange; ix++)
                for (iz = -srcrange; iz <= srcrange; iz++) {
                    u1[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
                    w1[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
                }

            //! Calculate veloctiy
            for (ix = marg; ix < nxb - marg; ix++)
                for (iz = marg; iz < nzb - marg; iz++) {
                    u2[ix][iz] = cp[ix][iz] * sjmfd2dn2(u1, ix, iz) + cs[ix][iz] * sjmfd2dn1(u1, ix, iz) +
                                 (cp[ix][iz] - cs[ix][iz]) * sjmfd2dnc(w1, ix, iz) + 2.0f * u1[ix][iz] - u0[ix][iz];

                    w2[ix][iz] = cp[ix][iz] * sjmfd2dn1(w1, ix, iz) + cs[ix][iz] * sjmfd2dn2(w1, ix, iz) +
                                 (cp[ix][iz] - cs[ix][iz]) * sjmfd2dnc(u1, ix, iz) + 2.0f * w1[ix][iz] - w0[ix][iz];
                }

            //! Boundary condition
            sjapplythabc2d(u2, u1, u0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
            sjapplythabc2d(w2, w1, w0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

            //! Record
            for (ir = 0; ir < nr; ir++) {
                wav->profx[ir][it] -= u1[nb + marg + rx[ir]][nb + marg + rz[ir]];
                wav->profz[ir][it] -= w1[nb + marg + rx[ir]][nb + marg + rz[ir]];
            }

            //! Update
            memcpy(u0[0], u1[0], nxb * nzb * sizeof(float));
            memcpy(u1[0], u2[0], nxb * nzb * sizeof(float));
            memcpy(w0[0], w1[0], nxb * nzb * sizeof(float));
            memcpy(w1[0], w2[0], nxb * nzb * sizeof(float));
        }
    }

    sjmfree1d(wavelet);

    sjmfree2d(cp);
    sjmfree2d(cs);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(u2);
    sjmfree2d(u1);
    sjmfree2d(u0);
    sjmfree2d(w2);
    sjmfree2d(w1);
    sjmfree2d(w0);
}

//! Two dimension constant density elastic forward exploration with SGFD
void sjesgfor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int it, ir, ix, iz;
    const int marg = 6;

    //------------------------ Option ------------------------//
    int nt = opt->nt;
    int k1 = opt->k1;
    int jsnap = opt->jsnap;
    int srcrange = opt->srcrange;
    int srctrunc = opt->srctrunc;
    float dt = opt->dt;
    float dt2 = opt->dt * opt->dt;
    float fp = opt->fp;
    float amp = opt->amp;
    float srcdecay = opt->srcdecay;
    int nb = opt->nb;
    float ds = opt->ds;
    int ycutdirect = opt->ycutdirect;
    float *wavelet = sjmflloc1d(nt);
    sjricker1d(wavelet, nt, k1, dt, fp, amp);

    //------------------------ Survey ------------------------//
    int nx = sur->nx;
    int nz = sur->nz;
    int sx = sur->sx + nb + marg;
    int sz = sur->sz + nb + marg;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids2 = dt2 / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    float **cs = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);
    sjextend2d(cs, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vs2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 8);
    float **gxr = sjmflloc2d(nzb, 8);
    float **gzu = sjmflloc2d(nxb, 8);
    float **gzb = sjmflloc2d(nxb, 8);

    //------------------------ Wavefield ------------------------//
    float **u2 = sjmflloc2d(nxb, nzb);
    float **u1 = sjmflloc2d(nxb, nzb);
    float **u0 = sjmflloc2d(nxb, nzb);
    float **w2 = sjmflloc2d(nxb, nzb);
    float **w1 = sjmflloc2d(nxb, nzb);
    float **w0 = sjmflloc2d(nxb, nzb);
    float **txx = sjmflloc2d(nxb, nzb);
    float **tzz = sjmflloc2d(nxb, nzb);
    float **txz = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitthabc2d(cp, cs, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, 1.0f / ds, cp[0], cp[0]);
    sjvecmulf(cs[0], nxb * nzb, 1.0f / ds, cs[0], cs[0]);
    //! Wavefield
    sjveczerof(u2[0], nxb * nzb);
    sjveczerof(u1[0], nxb * nzb);
    sjveczerof(u0[0], nxb * nzb);
    sjveczerof(w2[0], nxb * nzb);
    sjveczerof(w1[0], nxb * nzb);
    sjveczerof(w0[0], nxb * nzb);
    sjveczerof(txx[0], nxb * nzb);
    sjveczerof(tzz[0], nxb * nzb);
    sjveczerof(txz[0], nxb * nzb);

    //------------------------ Wavefield exploration ------------------------//
    for (it = 0; it < nt; it++) {

        //! Source
        for (ix = -srcrange; ix <= srcrange; ix++)
            for (iz = -srcrange; iz <= srcrange; iz++) {
                txx[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
                tzz[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
            }

        //! Calculate u and w
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++) {
                u2[ix][iz] = (sjmsgfd2dn1(txz, ix, iz - 1) + sjmsgfd2dn2(txx, ix, iz)) * ids2 + 2.0f * u1[ix][iz] -
                             u0[ix][iz];
                w2[ix][iz] = (sjmsgfd2dn1(tzz, ix, iz) + sjmsgfd2dn2(txz, ix - 1, iz)) * ids2 + 2.0f * w1[ix][iz] -
                             w0[ix][iz];
            }

        //! Calculate txx, tzz and txz
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++) {
                txx[ix][iz] = cp[ix][iz] * sjmsgfd2dn2(u2, ix - 1, iz) +
                              (cp[ix][iz] - 2.0f * cs[ix][iz]) * sjmsgfd2dn1(w2, ix, iz - 1);
                tzz[ix][iz] = cp[ix][iz] * sjmsgfd2dn1(w2, ix, iz - 1) +
                              (cp[ix][iz] - 2.0f * cs[ix][iz]) * sjmsgfd2dn2(u2, ix - 1, iz);
                txz[ix][iz] = cs[ix][iz] * (sjmsgfd2dn1(u2, ix, iz) + sjmsgfd2dn2(w2, ix, iz));
            }

        //! Boundary condition
        sjapplythabc2d(u2, u1, u0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        sjapplythabc2d(w2, w1, w0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
        for (ir = 0; ir < nr; ir++) {
            wav->profx[ir][it] = u1[nb + marg + rx[ir]][nb + marg + rz[ir]];
            wav->profz[ir][it] = w1[nb + marg + rx[ir]][nb + marg + rz[ir]];
        }

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++) {
                    wav->fwx2d[it / jsnap][ix - nb - marg][iz - nb - marg] = u1[ix][iz];
                    wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] = w1[ix][iz];
                }

        //! Update
        memcpy(u0[0], u1[0], nxb * nzb * sizeof(float));
        memcpy(u1[0], u2[0], nxb * nzb * sizeof(float));
        memcpy(w0[0], w1[0], nxb * nzb * sizeof(float));
        memcpy(w1[0], w2[0], nxb * nzb * sizeof(float));
    }

    //------------------------ Cut direct wav ------------------------//
    if (ycutdirect == 1) {
        //------------------------ Model ------------------------//
        for (ix = 0; ix < nxb; ix++)
            for (iz = 0; iz < nzb; iz++) {
                cp[ix][iz] = geo->vp2d[sx - marg - nb][sz - marg - nb];
                cs[ix][iz] = geo->vs2d[sx - marg - nb][sz - marg - nb];
            }

        //------------------------ Initialization ------------------------//
        //! Boundary condition
        sjinitthabc2d(cp, cs, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
        //! Model
        sjvecmulf(cp[0], nxb * nzb, 1.0f / ds, cp[0], cp[0]);
        sjvecmulf(cs[0], nxb * nzb, 1.0f / ds, cs[0], cs[0]);
        //! Wavefield
        sjveczerof(u2[0], nxb * nzb);
        sjveczerof(u1[0], nxb * nzb);
        sjveczerof(u0[0], nxb * nzb);
        sjveczerof(w2[0], nxb * nzb);
        sjveczerof(w1[0], nxb * nzb);
        sjveczerof(w0[0], nxb * nzb);
        sjveczerof(txx[0], nxb * nzb);
        sjveczerof(tzz[0], nxb * nzb);
        sjveczerof(txz[0], nxb * nzb);

        //! Wavefield exploration
        for (it = 0; it < nt; it++) {

            //! Source
            for (ix = -srcrange; ix <= srcrange; ix++)
                for (iz = -srcrange; iz <= srcrange; iz++) {
                    txx[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
                    tzz[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
                }

            //! Calculate u and w
            for (ix = marg; ix < nxb - marg; ix++)
                for (iz = marg; iz < nzb - marg; iz++) {
                    u2[ix][iz] = (sjmsgfd2dn1(txz, ix, iz - 1) + sjmsgfd2dn2(txx, ix, iz)) * ids2 + 2.0f * u1[ix][iz] -
                                 u0[ix][iz];
                    w2[ix][iz] = (sjmsgfd2dn1(tzz, ix, iz) + sjmsgfd2dn2(txz, ix - 1, iz)) * ids2 + 2.0f * w1[ix][iz] -
                                 w0[ix][iz];
                }

            //! Calculate txx, tzz and txz
            for (ix = marg; ix < nxb - marg; ix++)
                for (iz = marg; iz < nzb - marg; iz++) {
                    txx[ix][iz] = cp[ix][iz] * sjmsgfd2dn2(u2, ix - 1, iz) +
                                  (cp[ix][iz] - 2.0f * cs[ix][iz]) * sjmsgfd2dn1(w2, ix, iz - 1);
                    tzz[ix][iz] = cp[ix][iz] * sjmsgfd2dn1(w2, ix, iz - 1) +
                                  (cp[ix][iz] - 2.0f * cs[ix][iz]) * sjmsgfd2dn2(u2, ix - 1, iz);
                    txz[ix][iz] = cs[ix][iz] * (sjmsgfd2dn1(u2, ix, iz) + sjmsgfd2dn2(w2, ix, iz));
                }

            //! Boundary condition
            sjapplythabc2d(u2, u1, u0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
            sjapplythabc2d(w2, w1, w0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

            //! Record
            for (ir = 0; ir < nr; ir++) {
                wav->profx[ir][it] -= u1[nb + marg + rx[ir]][nb + marg + rz[ir]];
                wav->profz[ir][it] -= w1[nb + marg + rx[ir]][nb + marg + rz[ir]];
            }

            //! Update
            memcpy(u0[0], u1[0], nxb * nzb * sizeof(float));
            memcpy(u1[0], u2[0], nxb * nzb * sizeof(float));
            memcpy(w0[0], w1[0], nxb * nzb * sizeof(float));
            memcpy(w1[0], w2[0], nxb * nzb * sizeof(float));
        }
    }

    sjmfree1d(wavelet);

    sjmfree2d(cp);
    sjmfree2d(cs);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(u2);
    sjmfree2d(u1);
    sjmfree2d(u0);
    sjmfree2d(w2);
    sjmfree2d(w1);
    sjmfree2d(w0);

    sjmfree2d(txx);
    sjmfree2d(tzz);
    sjmfree2d(txz);
}

//! Two dimension constant density elastic forward exploration with V-S SGFD
void sjevssgfor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int it, ir, ix, iz;
    const int marg = 6;

    //------------------------ Option ------------------------//
    int nt = opt->nt;
    int k1 = opt->k1;
    int jsnap = opt->jsnap;
    int srcrange = opt->srcrange;
    float dt = opt->dt;
    float fp = opt->fp;
    float amp = opt->amp;
    float srcdecay = opt->srcdecay;
    int nb = opt->nb;
    float ds = opt->ds;
    int ycutdirect = opt->ycutdirect;
    float *wavelet = sjmflloc1d(nt);
    sjricker1d(wavelet, nt, k1, dt, fp, amp);

    //------------------------ Survey ------------------------//
    int nx = sur->nx;
    int nz = sur->nz;
    int sx = sur->sx + nb + marg;
    int sz = sur->sz + nb + marg;
    int nr = sur->nr;
    int *rx = sur->rx;
    int *rz = sur->rz;

    //------------------------ Model ------------------------//
    int nxb = nx + 2 * marg + 2 * nb;
    int nzb = nz + 2 * marg + 2 * nb;
    float ids = dt / ds;
    float **cp = sjmflloc2d(nxb, nzb);
    float **cs = sjmflloc2d(nxb, nzb);
    sjextend2d(cp, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vp2d);
    sjextend2d(cs, nx, nz, nb + marg, nb + marg, nb + marg, nb + marg, geo->vs2d);

    //------------------------ Boundary condition ------------------------//
    float **gxl = sjmflloc2d(nzb, 3);
    float **gxr = sjmflloc2d(nzb, 3);
    float **gzu = sjmflloc2d(nxb, 3);
    float **gzb = sjmflloc2d(nxb, 3);

    //------------------------ Wavefield ------------------------//
    float **vx0 = sjmflloc2d(nxb, nzb);
    float **vx1 = sjmflloc2d(nxb, nzb);
    float **vz0 = sjmflloc2d(nxb, nzb);
    float **vz1 = sjmflloc2d(nxb, nzb);
    float **txx0 = sjmflloc2d(nxb, nzb);
    float **txx1 = sjmflloc2d(nxb, nzb);
    float **tzz0 = sjmflloc2d(nxb, nzb);
    float **tzz1 = sjmflloc2d(nxb, nzb);
    float **txz0 = sjmflloc2d(nxb, nzb);
    float **txz1 = sjmflloc2d(nxb, nzb);

    //------------------------ Initialization ------------------------//
    //! Boundary condition
    sjinitohabc2d(cp, cs, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
    //! Model
    sjvecmulf(cp[0], nxb * nzb, 1.0f / ds, cp[0], cp[0]);
    sjvecmulf(cs[0], nxb * nzb, 1.0f / ds, cs[0], cs[0]);
    //! Wavefield
    sjveczerof(vx0[0], nxb * nzb);
    sjveczerof(vx1[0], nxb * nzb);
    sjveczerof(vz0[0], nxb * nzb);
    sjveczerof(vz1[0], nxb * nzb);
    sjveczerof(txx0[0], nxb * nzb);
    sjveczerof(txx1[0], nxb * nzb);
    sjveczerof(tzz0[0], nxb * nzb);
    sjveczerof(tzz1[0], nxb * nzb);
    sjveczerof(txz0[0], nxb * nzb);
    sjveczerof(txz1[0], nxb * nzb);

    //------------------------ Wavefield exploration ------------------------//
    for (it = 0; it < nt; it++) {

        //! Source
        for (ix = -srcrange; ix <= srcrange; ix++)
            for (iz = -srcrange; iz <= srcrange; iz++) {
                txx0[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
                tzz0[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
            }

        //! Calculate vx and vz
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++) {
                vx1[ix][iz] = (sjmsgfd2dn1(txz0, ix, iz - 1) + sjmsgfd2dn2(txx0, ix, iz)) * ids + vx0[ix][iz];
                vz1[ix][iz] = (sjmsgfd2dn1(tzz0, ix, iz) + sjmsgfd2dn2(txz0, ix - 1, iz)) * ids + vz0[ix][iz];
            }

        //! Calculate txx, tzz and txz
        for (ix = marg; ix < nxb - marg; ix++)
            for (iz = marg; iz < nzb - marg; iz++) {
                txx1[ix][iz] = cp[ix][iz] * sjmsgfd2dn2(vx1, ix - 1, iz) +
                               (cp[ix][iz] - 2.0f * cs[ix][iz]) * sjmsgfd2dn1(vz1, ix, iz - 1) + txx0[ix][iz];
                tzz1[ix][iz] = cp[ix][iz] * sjmsgfd2dn1(vz1, ix, iz - 1) +
                               (cp[ix][iz] - 2.0f * cs[ix][iz]) * sjmsgfd2dn2(vx1, ix - 1, iz) + tzz0[ix][iz];
                txz1[ix][iz] = cs[ix][iz] * (sjmsgfd2dn1(vx1, ix, iz) + sjmsgfd2dn2(vz1, ix, iz)) + txz0[ix][iz];
            }

        //! Boundary condition
        //sjapplyohabc2d(vx1, vx0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        //sjapplyohabc2d(vz1, vz0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        //sjapplyohabc2d(txx1, txx0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        //sjapplyohabc2d(tzz1, tzz0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
        //sjapplyohabc2d(txz1, txz0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

        //! Record
        for (ir = 0; ir < nr; ir++) {
            wav->profx[ir][it] = vx1[nb + marg + rx[ir]][nb + marg + rz[ir]];
            wav->profz[ir][it] = vz1[nb + marg + rx[ir]][nb + marg + rz[ir]];
        }

        //! Wavefield
        if ((it % jsnap) == 0)
            for (ix = nb + marg; ix < nxb - nb - marg; ix++)
                for (iz = nb + marg; iz < nzb - nb - marg; iz++) {
                    wav->fwx2d[it / jsnap][ix - nb - marg][iz - nb - marg] = vx1[ix][iz];
                    wav->fwz2d[it / jsnap][ix - nb - marg][iz - nb - marg] = vz1[ix][iz];
                }

        //! Update
        memcpy(vx0[0], vx1[0], nxb * nzb * sizeof(float));
        memcpy(vz0[0], vz1[0], nxb * nzb * sizeof(float));
        memcpy(txx0[0], txx1[0], nxb * nzb * sizeof(float));
        memcpy(tzz0[0], tzz1[0], nxb * nzb * sizeof(float));
        memcpy(txz0[0], txz1[0], nxb * nzb * sizeof(float));
    }

    //------------------------ Cut direct wav ------------------------//
    if (ycutdirect == 1) {
        //------------------------ Model ------------------------//
        for (ix = 0; ix < nxb; ix++)
            for (iz = 0; iz < nzb; iz++) {
                cp[ix][iz] = geo->vp2d[sx - marg - nb][sz - marg - nb];
                cs[ix][iz] = geo->vs2d[sx - marg - nb][sz - marg - nb];
            }

        //------------------------ Initialization ------------------------//
        //! Boundary condition
        sjinitohabc2d(cp, cs, ds, dt, nxb, nzb, gxl, gxr, gzu, gzb);
        //! Model
        sjvecmulf(cp[0], nxb * nzb, 1.0f / ds, cp[0], cp[0]);
        sjvecmulf(cs[0], nxb * nzb, 1.0f / ds, cs[0], cs[0]);
        //! Wavefield
        sjveczerof(vx0[0], nxb * nzb);
        sjveczerof(vx1[0], nxb * nzb);
        sjveczerof(vz0[0], nxb * nzb);
        sjveczerof(vz1[0], nxb * nzb);
        sjveczerof(txx0[0], nxb * nzb);
        sjveczerof(txx1[0], nxb * nzb);
        sjveczerof(tzz0[0], nxb * nzb);
        sjveczerof(tzz1[0], nxb * nzb);
        sjveczerof(txz0[0], nxb * nzb);
        sjveczerof(txz1[0], nxb * nzb);

        //------------------------ Wavefield exploration ------------------------//
        for (it = 0; it < nt; it++) {

            //! Source
            for (ix = -srcrange; ix <= srcrange; ix++)
                for (iz = -srcrange; iz <= srcrange; iz++) {
                    txx0[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
                    tzz0[ix + sx][iz + sz] += wavelet[it] * expf(-srcdecay * (ix * ix + iz * iz));
                }

            //! Calculate vx and vz
            for (ix = marg; ix < nxb - marg; ix++)
                for (iz = marg; iz < nzb - marg; iz++) {
                    vx1[ix][iz] = (sjmsgfd2dn1(txz0, ix, iz - 1) + sjmsgfd2dn2(txx0, ix, iz)) * ids + vx0[ix][iz];
                    vz1[ix][iz] = (sjmsgfd2dn1(tzz0, ix, iz) + sjmsgfd2dn2(txz0, ix - 1, iz)) * ids + vz0[ix][iz];
                }

            //! Calculate txx, tzz and txz
            for (ix = marg; ix < nxb - marg; ix++)
                for (iz = marg; iz < nzb - marg; iz++) {
                    txx1[ix][iz] = cp[ix][iz] * sjmsgfd2dn2(vx1, ix - 1, iz) +
                                   (cp[ix][iz] - 2.0f * cs[ix][iz]) * sjmsgfd2dn1(vz1, ix, iz - 1) + txx0[ix][iz];
                    tzz1[ix][iz] = cp[ix][iz] * sjmsgfd2dn1(vz1, ix, iz - 1) +
                                   (cp[ix][iz] - 2.0f * cs[ix][iz]) * sjmsgfd2dn2(vx1, ix - 1, iz) + tzz0[ix][iz];
                    txz1[ix][iz] = cs[ix][iz] * (sjmsgfd2dn1(vx1, ix, iz) + sjmsgfd2dn2(vz1, ix, iz)) + txz0[ix][iz];
                }

            //! Boundary condition
            //sjapplyohabc2d(vx1, vx0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
            //sjapplyohabc2d(vz1, vz0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
            //sjapplyohabc2d(txx1, txx0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
            //sjapplyohabc2d(tzz1, tzz0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);
            //sjapplyohabc2d(txz1, txz0, gxl, gxr, gzu, gzb, nxb, nzb, nb, marg);

            //! Record
            for (ir = 0; ir < nr; ir++) {
                wav->profx[ir][it] -= vx1[nb + marg + rx[ir]][nb + marg + rz[ir]];
                wav->profz[ir][it] -= vz1[nb + marg + rx[ir]][nb + marg + rz[ir]];
            }

            //! Update
            memcpy(vx0[0], vx1[0], nxb * nzb * sizeof(float));
            memcpy(vz0[0], vz1[0], nxb * nzb * sizeof(float));
            memcpy(txx0[0], txx1[0], nxb * nzb * sizeof(float));
            memcpy(tzz0[0], tzz1[0], nxb * nzb * sizeof(float));
            memcpy(txz0[0], txz1[0], nxb * nzb * sizeof(float));
        }
    }

    sjmfree1d(wavelet);

    sjmfree2d(cp);
    sjmfree2d(cs);

    sjmfree2d(gxl);
    sjmfree2d(gxr);
    sjmfree2d(gzu);
    sjmfree2d(gzb);

    sjmfree2d(vx0);
    sjmfree2d(vx1);
    sjmfree2d(vz0);
    sjmfree2d(vz1);
    sjmfree2d(txx0);
    sjmfree2d(txx1);
    sjmfree2d(tzz0);
    sjmfree2d(tzz1);
    sjmfree2d(txz0);
    sjmfree2d(txz1);
}
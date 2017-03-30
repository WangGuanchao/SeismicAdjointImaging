//
// Authors: Hou, Sian - sianhou1987@outlook.com
//          Wang, Guangchao - wgcupc@163.com

#ifndef SJI_SJIMAGE_H
#define SJI_SJIMAGE_H

#include "../lib/sjinc.h"
#include <mpi.h>

//! Two dimension constant density acoustic RTM
void sjartm2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int is = 0;
    double tstart, tend;
    int rankid, nrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //------------------------ Model ------------------------//
    geo->gnzz2d = sjmflloc2d(sur->gnx, sur->gnz);
    MPI_Bcast(geo->gvp2d[0], sur->gnx * sur->gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //------------------------ Calculating gradient ------------------------//
    //! Informaiton
    if (rankid == 0) {
        tstart = (double) clock();
        printf("Calculating RTM -                     ");
    }

    for (is = rankid; is < sur->ns; is += nrank) {
        //! Set survey
        sjssurvey_readis(sur, is);
        //! Memory
        geo->vp2d = sjmflloc2d(sur->nx, sur->nz);
        geo->izz2d = sjmflloc2d(sur->nx, sur->nz);
        geo->nzz2d = sjmflloc2d(sur->nx, sur->nz);
        wav->profz = sjmflloc2d(sur->nr, opt->nt);
        wav->fwz2d = sjmflloc3d(opt->nsnap, sur->nx, sur->nz);
        //! Set model
        sjextract2d(geo->vp2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->gvp2d);
        //! Forward exploration
        opt->ycutdirect = 0;
        sjafor2d(sur, geo, wav, opt);
        //! Read record
        sjreadsu(wav->profz[0], sur->nr, opt->nt, sizeof(float), sur->tr, 0, wav->profzfile);
        //! Backward exploration
        opt->ystacksrc = 0;
        sjartmbac2d(sur, geo, wav, opt);
        //! Stacking
        sjvecaddf(&geo->gizz2d[sur->x0][sur->z0], sur->nx * sur->nz, 1.0f, &geo->gizz2d[sur->x0][sur->z0], 1.0f,
                  &geo->izz2d[0][0]);
        sjvecaddf(&geo->gnzz2d[sur->x0][sur->z0], sur->nx * sur->nz, 1.0f, &geo->gnzz2d[sur->x0][sur->z0], 1.0f,
                  &geo->nzz2d[0][0]);
        //! Free
        sjmfree2d(geo->vp2d);
        sjmfree2d(geo->izz2d);
        sjmfree2d(geo->nzz2d);
        sjmfree2d(wav->profz);
        sjmfree3d(wav->fwz2d);
        //! Informaiton
        if (rankid == 0) {
            tend = (double) clock();
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("                     ");
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf(" %4d/ %4d - %6.2fs", (is + 1) * nrank, sur->ns, (tend - tstart) / CLOCKS_PER_SEC);
        }
    }

    if (rankid == 0) {
        //! Communication
        MPI_Reduce(MPI_IN_PLACE, geo->gizz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, geo->gnzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        //! Source
        sjvecdivf(geo->gizz2d[0], sur->gnx * sur->gnz, 1.0, geo->gizz2d[0], geo->gnzz2d[0], 0.00001f);
        //! Cut source
        sjsetsurface(geo->ggzz2d, sur->gnx, 30, 0.0);
        //! Information
        tend = (double) clock();
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("                     ");
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf(" complete - %6.2fs.\n", (tend - tstart) / CLOCKS_PER_SEC);
    } else {
        //! Communication
        MPI_Reduce(geo->gizz2d[0], geo->gizz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(geo->gnzz2d[0], geo->gnzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    //! Free
    sjmfree2d(geo->gnzz2d);
}

//! Two dimension constant density acoustic Time-Shift RTM
void sjatsrtm2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int is = 0, shift;
    double tstart, tend;
    int rankid, nrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //------------------------ Model ------------------------//
    MPI_Bcast(geo->gvp2d[0], sur->gnx * sur->gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //------------------------ Calculating gradient ------------------------//
    //! Informaiton
    if (rankid == 0) {
        tstart = (double) clock();
        printf("Calculating Time-Shift RTM -                     ");
    }

    for (is = rankid; is < sur->ns; is += nrank) {

        //! Set survey
        sjssurvey_readis(sur, is);

        //! Memory
        geo->vp2d = sjmflloc2d(sur->nx, sur->nz);
        wav->profz = sjmflloc2d(sur->nr, opt->nt);
        wav->fwz2d = sjmflloc3d(opt->nsnap, sur->nx, sur->nz);
        wav->bwz2d = sjmflloc3d(opt->nsnap, sur->nx, sur->nz);
        geo->izz3d = sjmflloc3d(2 * opt->maxshift + 1, sur->nx, sur->nz);

        //! Set model
        sjextract2d(geo->vp2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->gvp2d);

        //! Forward exploration

        sjafor2d(sur, geo, wav, opt);

        //! Read record
        sjreadsu(wav->profz[0], sur->nr, opt->nt, sizeof(float), sur->tr, 0, wav->profzfile);

        //! Backward exploration
        opt->ystacksrc = 0;
        sjatsrtmbac2d(sur, geo, wav, opt);

        //! Time-Shift stack
        for (shift = -opt->maxshift; shift <= opt->maxshift; ++shift) {
            sjvecaddf(&geo->gizz3d[shift + opt->maxshift][sur->x0][sur->z0], sur->nx * sur->nz, 1.0f,
                      &geo->gizz3d[shift + opt->maxshift][sur->x0][sur->z0], 1.0f,
                      &geo->izz3d[shift + opt->maxshift][0][0]);
        }

        //! Free
        sjmfree2d(geo->vp2d);
        sjmfree2d(wav->profz);
        sjmfree3d(wav->fwz2d);
        sjmfree3d(wav->bwz2d);
        sjmfree3d(geo->izz3d);

        //! Informaiton
        if (rankid == 0) {
            tend = (double) clock();
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("                     ");
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf(" %4d/ %4d - %6.2fs", (is + 1) * nrank, sur->ns, (tend - tstart) / CLOCKS_PER_SEC);
        }
    }

    if (rankid == 0) {
        //! Communication
        MPI_Reduce(MPI_IN_PLACE, geo->gizz3d[0][0], (2 * opt->maxshift + 1) * sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM,
                   0, MPI_COMM_WORLD);

        //! Information
        tend = (double) clock();
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("                     ");
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf(" complete - %6.2fs.\n", (tend - tstart) / CLOCKS_PER_SEC);
    } else {
        //! Communication
        MPI_Reduce(geo->gizz3d[0][0], geo->gizz3d[0][0], (2 * opt->maxshift + 1) * sur->gnx * sur->gnz, MPI_FLOAT,
                   MPI_SUM, 0, MPI_COMM_WORLD);
    }
}

//! Two dimension constant density acoustic FWI gradient
void sjafwig2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int is = 0;
    double tstart, tend;
    int rankid, nrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //------------------------ Bcast model ------------------------//
    geo->gnzz2d = sjmflloc2d(sur->gnx, sur->gnz);
    MPI_Bcast(geo->gvp2d[0], sur->gnx * sur->gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //------------------------ Calculating gradient ------------------------//
    //! Informaiton
    if (rankid == 0) {
        tstart = (double) clock();
        printf("Calculating AFWI gradient -                     ");
    }

    //! Optimization
    for (is = rankid; is < sur->ns; is += nrank) {

        //! Memory
        geo->vp2d = sjmflloc2d(sur->nx, sur->nz);
        wav->profz = sjmflloc2d(sur->nr, opt->nt);
        geo->gzz2d = sjmflloc2d(sur->nx, sur->nz);
        geo->nzz2d = sjmflloc2d(sur->nx, sur->nz);
        float **recz = sjmflloc2d(sur->nr, opt->nt);
        wav->fwz2d = sjmflloc3d(opt->nsnap, sur->nx, sur->nz);

        //! Set survey
        sjssurvey_readis(sur, is);

        //! Set model
        sjextract2d(geo->vp2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->gvp2d);

        //! Simulation
        sjafor2d(sur, geo, wav, opt);

        //! Difference wavefield
        sjreadsu(recz[0], sur->nr, opt->nt, sizeof(float), sur->tr, 0, wav->profzfile);
        sjvecsubf(wav->profz[0], sur->nr * opt->nt, 1.0f, recz[0], 1.0f, wav->profz[0]);

        //! Adjoint image
        sjafwibac2d(sur, geo, wav, opt);

        //! Stacking
        sjvecaddf(geo->ggzz2d[sur->x0], sur->nx * sur->nz, 1.0f, geo->ggzz2d[sur->x0], 1.0f, geo->gzz2d[0]);
        sjvecaddf(geo->gnzz2d[sur->x0], sur->nx * sur->nz, 1.0f, geo->gnzz2d[sur->x0], 1.0f, geo->nzz2d[0]);

        //! Free
        sjmfree2d(geo->vp2d);
        sjmfree2d(wav->profz);
        sjmfree2d(geo->gzz2d);
        sjmfree2d(geo->nzz2d);
        sjmfree2d(recz);
        sjmfree3d(wav->fwz2d);

        //! Informaiton
        if (rankid == 0) {
            tend = (double) clock();
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("                     ");
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf(" %4d/ %4d - %6.2fs", (is + 1) * nrank, sur->ns, (tend - tstart) / CLOCKS_PER_SEC);
        }
    }

    if (rankid == 0) {
        //! Communication
        MPI_Reduce(MPI_IN_PLACE, geo->ggzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, geo->gnzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

        //! Source
        sjvecdivf(geo->ggzz2d[0], sur->gnx * sur->gnz, 1.0, geo->ggzz2d[0], geo->gnzz2d[0], 0.00001f);

        //! Cut source
        sjsetsurface(geo->ggzz2d, sur->gnx, 30, 0.0);

        //! Information
        tend = (double) clock();
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("                     ");
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf(" complete - %6.2fs.\n", (tend - tstart) / CLOCKS_PER_SEC);
    } else {
        //! Communication
        MPI_Reduce(geo->ggzz2d[0], geo->ggzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(geo->gnzz2d[0], geo->gnzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    sjmfree2d(geo->gnzz2d);
}

//! Two dimension constant density acoustic LSRTM gradient
void sjalsrtmg2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int is = 0;
    double tstart, tend;
    int rankid, nrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //------------------------ Bcast model ------------------------//
    geo->gnzz2d = sjmflloc2d(sur->gnx, sur->gnz);
    MPI_Bcast(geo->gvp2d[0], sur->gnx * sur->gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(geo->gizz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //------------------------ Calculating gradient ------------------------//
    //! Informaiton
    if (rankid == 0) {
        tstart = (double) clock();
        printf("Calculating ALSRTM gradient -                     ");
    }

    //! Optimization
    for (is = rankid; is < sur->ns; is += nrank) {
        //! Memory
        geo->vp2d = sjmflloc2d(sur->nx, sur->nz);
        geo->izz2d = sjmflloc2d(sur->nx, sur->nz);
        geo->gzz2d = sjmflloc2d(sur->nx, sur->nz);
        geo->nzz2d = sjmflloc2d(sur->nx, sur->nz);
        wav->profz = sjmflloc2d(sur->nr, opt->nt);
        wav->fwz2d = sjmflloc3d(opt->nsnap, sur->nx, sur->nz);
        float **recz = sjmflloc2d(sur->nr, opt->nt);

        //! Set survey
        sjssurvey_readis(sur, is);

        //! Set model
        sjextract2d(geo->vp2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->gvp2d);
        sjextract2d(geo->izz2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->gizz2d);

        //! Simulation
        sjalsrtmfor2d(sur, geo, wav, opt);

        //! Difference
        sjreadsu(recz[0], sur->nr, opt->nt, sizeof(float), sur->tr, 0, wav->profzfile);
        sjvecsubf(wav->profz[0], sur->nr * opt->nt, 1.0f, wav->profz[0], 1.0f, recz[0]);

        //! Imaging
        sjalsrtmbac2d(sur, geo, wav, opt);

        //! Stacking
        sjvecaddf(geo->ggzz2d[sur->x0], sur->nx * sur->nz, 1.0f, geo->ggzz2d[sur->x0], 1.0f, geo->gzz2d[0]);
        sjvecaddf(geo->gnzz2d[sur->x0], sur->nx * sur->nz, 1.0f, geo->gnzz2d[sur->x0], 1.0f, geo->nzz2d[0]);

        //! Free
        sjmfree2d(geo->vp2d);
        sjmfree2d(geo->izz2d);
        sjmfree2d(geo->gzz2d);
        sjmfree2d(geo->nzz2d);
        sjmfree2d(wav->profz);
        sjmfree2d(recz);
        sjmfree3d(wav->fwz2d);

        //! Informaiton
        if (rankid == 0) {
            tend = (double) clock();
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("                     ");
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf(" %4d/ %4d - %6.2fs", (is + 1) * nrank, sur->ns, (tend - tstart) / CLOCKS_PER_SEC);
        }
    }

    if (rankid == 0) {
        //! Communication
        MPI_Reduce(MPI_IN_PLACE, geo->ggzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, geo->gnzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

        //! Source
        sjvecdivf(geo->ggzz2d[0], sur->gnx * sur->gnz, 1.0, geo->ggzz2d[0], geo->gnzz2d[0], 0.00001f);

        //! Cut source
        sjsetsurface(geo->ggzz2d, sur->gnx, 30, 0.0);

        //! Information
        tend = (double) clock();
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("                     ");
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf(" complete - %6.2fs.\n", (tend - tstart) / CLOCKS_PER_SEC);
    } else {
        //! Communication
        MPI_Reduce(geo->ggzz2d[0], geo->ggzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(geo->gnzz2d[0], geo->gnzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    sjmfree2d(geo->gnzz2d);
}

//! Two dimension constant density acoustic RTI sfhit
void sjartishift(float **pobs, float **pcal, int nr, int nt, int maxshift, float dt) {

    int ir, it, its;

    float *cor = sjmflloc1d(2 * maxshift + 1), temp;
    int index = 0;

    for (ir = 0; ir < nr; ++ir) {

        //! Calculate timeshift correlation
        sjveczerof(cor, 2 * maxshift + 1);
        for (its = -maxshift; its <= maxshift; ++its)
            for (it = abs(its); it < nt - abs(its); ++it) {
                cor[its + maxshift] += pobs[ir][it + its] * pcal[ir][it];
            }

        //! Calculate timeshift maxmium
        temp = cor[0];
        index = -maxshift;
        for (its = -maxshift; its <= maxshift; ++its)
            if (temp > cor[its + maxshift]) {
                temp = cor[its + maxshift];
                index = its;
            }

        sjveczerof(pcal[ir], nt);
        for (it = abs(index) + 1; it < nt - abs(index) - 1; ++it)
            pcal[ir][it] = index * dt * (pobs[ir][it + index + 1] - pobs[ir][it + index - 1]) / (2.0f * dt);
    }

    sjmfree1d(cor);
}

//! Two dimension constant density acoustic RTI gradient
void sjartig2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt) {

    //------------------------ Runtime ------------------------//
    int is = 0;
    double tstart, tend;
    int rankid, nrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //------------------------ Model ------------------------//
    MPI_Bcast(geo->gizz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(geo->gvp2d[0], sur->gnx * sur->gnz, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //------------------------ Calculating gradient ------------------------//
    //! Informaiton
    if (rankid == 0) {
        tstart = (double) clock();
        printf("Calculating ARTI Gradient -                     ");
    }

    for (is = rankid; is < sur->ns; is += nrank) {

        //! Set survey
        sjssurvey_readis(sur, is);

        //! Memory
        geo->vp2d = sjmflloc2d(sur->nx, sur->nz);
        geo->izz2d = sjmflloc2d(sur->nx, sur->nz);
        geo->gzz2d = sjmflloc2d(sur->nx, sur->nz);
        float **profz = sjmflloc2d(sur->nr, opt->nt);
        wav->profz = sjmflloc2d(sur->nr, opt->nt);
        wav->fwz2d = sjmflloc3d(opt->nsnap, sur->nx, sur->nz);
        wav->fsz2d = sjmflloc3d(opt->nsnap, sur->nx, sur->nz);

        //! Set model
        sjextract2d(geo->vp2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->gvp2d);
        sjextract2d(geo->izz2d, sur->x0, sur->z0, sur->nx, sur->nz, geo->gizz2d);

        //! Forward exploration
        sjasfor2d(sur, geo, wav, opt);

        //! Calculate time shift
        sjreadsu(profz[0], sur->nr, opt->nt, sizeof(float), sur->tr, 0, wav->profzfile);
        sjartishift(profz, wav->profz, sur->nr, opt->nt, opt->maxshift, opt->dt);

        //! Backward exploration
        opt->ystacksrc = 1;
        sjartibac2d(sur, geo, wav, opt);

        //! Stacking
        sjvecaddf(geo->ggzz2d[sur->x0], sur->nx * sur->nz, 1.0f, geo->ggzz2d[sur->x0], 1.0f, geo->gzz2d[0]);

        //! Free
        sjmfree2d(geo->vp2d);
        sjmfree2d(geo->izz2d);
        sjmfree2d(geo->gzz2d);
        sjmfree2d(profz);
        sjmfree2d(wav->profz);
        sjmfree3d(wav->fwz2d);
        sjmfree3d(wav->fsz2d);

        //! Informaiton
        if (rankid == 0) {
            tend = (double) clock();
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf("                     ");
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            printf(" %4d/ %4d - %6.2fs", (is + 1) * nrank, sur->ns, (tend - tstart) / CLOCKS_PER_SEC);
        }
    }

    if (rankid == 0) {
        //! Communication
        MPI_Reduce(MPI_IN_PLACE, geo->ggzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

        //! Information
        tend = (double) clock();
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("                     ");
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf(" complete - %6.2fs.\n", (tend - tstart) / CLOCKS_PER_SEC);
    } else {
        //! Communication
        MPI_Reduce(geo->ggzz2d[0], geo->ggzz2d[0], sur->gnx * sur->gnz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
}

#endif //SJI_SJIMAGE_H

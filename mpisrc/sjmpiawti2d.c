//
// Authors: Hou, Sian - sianhou1987@outlook.com
//          Wang, Guangchao - wgcupc@163.com

#include <mpi.h>
#include "sjimage.h"
#include "../lib/sjinc.h"






int main(int argc, char *argv[]) {

    //! Runtime
    int iter = 0, flag = 1;
    double Tstart, Tend;

    //! MPI
    int rankid, nrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //! Survey
    sjssurvey sur;
    flag &= sjssurvey_init(&sur);
    flag &= sjssurvey_getparas(&sur, argc, argv);

    //! Model
    sjsgeology geo;
    flag &= sjsgeo_init(&geo);
    flag &= sjsgeo_getparas2d(&geo, argc, argv, "vp");
    flag &= sjsgeo_getparas2d(&geo, argc, argv, "izz");

    //! Wave
    sjswave wav;
    flag &= sjswave_init(&wav);
    flag &= sjswave_getparas(&wav, argc, argv, "profz");

    //! Option
    sjsoption opt;
    flag &= sjsoption_init(&opt);
    flag &= sjsoption_getparas(&opt, argc, argv);

    //------------------------ AWTI2D ------------------------//
    if (flag) {
        //! Time
        if (rankid == 0) {
            Tstart = (double) clock();
            printf("------------------------ 2D Acoustic WTI start ------------------------\n");
        }

        opt.niter = 0;
        opt.maxshift = 200;

        //! Allocate memory
        geo.gvp2d = sjmflloc2d(sur.gnx, sur.gnz);
        geo.ggzz2d = sjmflloc2d(sur.gnx, sur.gnz);
        geo.gizz3d = sjmflloc3d(2 * opt.maxshift + 1, sur.gnx, sur.gnz);
        geo.ggzz3d = sjmflloc3d(2 * opt.maxshift + 1, sur.gnx, sur.gnz);

        //! Read model
        sjveczerof(geo.gvp2d[0], sur.gnx * sur.gnz);
        sjreadsuall(geo.gvp2d[0], sur.gnx, sur.gnz, geo.vpfile);

        //! Process optimization
        float *coe = sjmflloc1d(2 * opt.maxshift + 1);
        float **g0 = sjmflloc2d(sur.gnx, sur.gnz);
        float **cg = sjmflloc2d(sur.gnx, sur.gnz);
/*
        //! Inversion
        do {
            //! Calculate Time-Shift RTM image
            opt.ycutdirect = 0;
            opt.ystacksrc = 0;
            sjveczerof(coe, 2 * opt.maxshift + 1);
            sjveczerof(geo.gizz3d[0][0], (2 * opt.maxshift + 1) * sur.gnx * sur.gnz);
            sjatsrtm2d(&sur, &geo, &wav, &opt);

            //! Calculate RTI gradient
            sjveczerof(geo.ggzz2d[0], sur.gnx * sur.gnz);
            sjveczerof(geo.ggzz3d[0][0], (2 * opt.maxshift + 1) * sur.gnx * sur.gnz);
            sjawti2d(&sur, &geo, &wav, &opt);

            //! Optimization
            if (rankid == 0) {

                sjwritesu(geo.gizz3d[0][0], (2 * opt.maxshift + 1) * sur.gnx, sur.gnz, sizeof(float), opt.ds, 0,
                          "tsipp.su");

                sjwritesu(geo.ggzz3d[0][0], (2 * opt.maxshift + 1) * sur.gnx, sur.gnz, sizeof(float), opt.ds, 0,
                          "tsjad.su");


                int shift, ix, iz;
                for (shift = -opt.maxshift; shift < opt.maxshift; ++shift)
                    for (ix = 0; ix < sur.gnx; ++ix)
                        for (iz = 0; iz < sur.gnz; ++iz)
                            coe[shift + opt.maxshift] += shift * geo.gizz3d[shift + opt.maxshift][ix][iz] *
                                                         geo.gizz3d[shift + opt.maxshift][ix][iz];


                sjwritesu(coe, 1, 2 * opt.maxshift + 1, sizeof(float), opt.ds, 0, "coe.su");

                for (shift = -opt.maxshift; shift < opt.maxshift; ++shift)
                    for (ix = 0; ix < sur.gnx; ++ix)
                        for (iz = 0; iz < sur.gnz; ++iz)
                            geo.ggzz2d[ix][iz] += coe[shift + opt.maxshift] * geo.ggzz3d[shift + opt.maxshift][ix][iz];


                sjwritesu(geo.ggzz2d[0], sur.gnx, sur.gnz, sizeof(float), opt.ds, 0, "ggzz2d.su");

                //! Information
                Tend = (double) clock();
                printf("Acoustic AWTI complete - %2d/%2d - time=%6.2fs.\n",
                       iter + 1, opt.niter + 1, (Tend - Tstart) / CLOCKS_PER_SEC);
            }

            iter += 1;

        } while (iter < opt.niter);
*/
        sjmfree2d(geo.gvp2d);
        sjmfree2d(geo.ggzz2d);
        sjmfree3d(geo.gizz3d);
        sjmfree3d(geo.ggzz3d);
        sjmfree2d(cg);
        sjmfree2d(g0);

        if (rankid == 0) {
            printf("Acoustic WTI completed.\n\n");
        }

/*
        //! Inversion
        do {
            //! Optimization
            if (rankid == 0) {

                int ix, iz;
                float maxamp = 0.0;
                for (ix = 0; ix < sur.gnx; ix++) {
                    for (iz = 0; iz < sur.gnz; iz++) {
                        if (fabs(geo.ggzz2d[ix][iz]) > maxamp) {
                            maxamp = fabs(geo.ggzz2d[ix][iz]);
                        }
                    }
                }
                for (ix = 0; ix < sur.gnx; ix++) {
                    for (iz = 0; iz < sur.gnz; iz++) {
                        geo.ggzz2d[ix][iz] = geo.ggzz2d[ix][iz] / maxamp * 10.0f;
                    }
                }

                //! CG
                sjcgsolver(geo.gvp2d[0], sur.gnx * sur.gnz, cg[0], geo.ggzz2d[0], g0[0], iter);

                char *file = (char *) malloc(1024 * sizeof(char));

                sprintf(file, "%s-gvp.su", geo.izzfile);
                sjwritesu(geo.gvp2d[0], sur.gnx, sur.gnz, sizeof(float), opt.ds, iter, file);
                sprintf(file, "%s-grad.su", geo.izzfile);
                sjwritesu(geo.ggzz2d[0], sur.gnx, sur.gnz, sizeof(float), opt.ds, iter, file);
                sprintf(file, "%s-cg.su", geo.izzfile);
                sjwritesu(cg[0], sur.gnx, sur.gnz, sizeof(float), opt.ds, iter, file);
                sprintf(file, "%s-rtm.su", geo.izzfile);
                sjwritesu(geo.gizz2d[0], sur.gnx, sur.gnz, sizeof(float), opt.ds, iter, file);
            }
        } while (iter < opt.niter);
*/
    } else {
        if (rankid == 0) {
            printf("\nExamples:   sjmpilsartm2d survey=survey.su vp=vp.su profz=profz.su ipp=lsipp.su\n");

            sjbasicinformation();

        }
    }

    //------------------------ MPI finish ------------------------//
    MPI_Finalize();

    return 0;
}


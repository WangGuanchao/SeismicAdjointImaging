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

    //------------------------ LSARTM2D ------------------------//
    if (flag) {
        //! Time
        if (rankid == 0) {
            Tstart = (double) clock();
            printf("------------------------ 2D Acoustic LSRTM start ------------------------\n");
        }

        //! Set model
        geo.gvp2d = sjmflloc2d(sur.gnx, sur.gnz);
        geo.gizz2d = sjmflloc2d(sur.gnx, sur.gnz);
        sjreadsuall(geo.gvp2d[0], sur.gnx, sur.gnz, geo.vpfile);

        //! Process migration
        geo.ggzz2d = sjmflloc2d(sur.gnx, sur.gnz);
        float **g0 = sjmflloc2d(sur.gnx, sur.gnz);
        float **cg = sjmflloc2d(sur.gnx, sur.gnz);

        //! Inversion
        do {
            //! Gradient
            if (iter == 0) {
                opt.ystacksrc = 0;
                sjalsrtmg2d(&sur, &geo, &wav, &opt);
            } else {
                opt.ystacksrc = 1;
                sjalsrtmg2d(&sur, &geo, &wav, &opt);
            }

            //! Optimization
            if (rankid == 0) {
                //! CG
                sjcgsolver(geo.gizz2d[0], sur.gnx * sur.gnz, cg[0], geo.ggzz2d[0], g0[0], iter);

                //! Output details
                if(opt.ydetails==1) {
                    char *file = (char *) malloc(1024 * sizeof(char));
                    sprintf(file, "%s-cg", geo.izzfile);
                    sjwritesu(cg[0], sur.gnx, sur.gnz, sizeof(float), opt.ds, iter,file);
                    sprintf(file, "%s-izz", geo.izzfile);
                    sjwritesu(geo.gizz2d[0], sur.gnx, sur.gnz, sizeof(float), opt.ds, iter,file);
                }

                //! Information
                Tend = (double) clock();
                printf("Acoustic LSRTM complete - %2d/%2d - time=%6.2fs.\n",
                       iter + 1, opt.niter, (Tend - Tstart) / CLOCKS_PER_SEC);
            }
            iter += 1;
        } while (iter < opt.niter);

        //! Output
        if (rankid == 0) {
            sjwritesuall(geo.gizz2d[0], sur.gnx, sur.gnz, opt.ds, geo.izzfile);
            printf("Acoustic LSRTM completed.\n\n");
        }

        sjmfree2d(geo.gvp2d);
        sjmfree2d(geo.gizz2d);
        sjmfree2d(geo.ggzz2d);
        sjmfree2d(cg);
        sjmfree2d(g0);
    } else {
        if (rankid == 0) {
            printf("\nExamples:   sjmpialsrtm2d survey=survey.su vp=vp.su profz=profz.su izz=lsipp.su\n");
            sjbasicinformation();
        }
    }

    //------------------------ MPI finish ------------------------//
    MPI_Finalize();

    return 0;
}


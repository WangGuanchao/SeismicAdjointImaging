//
// Authors: Hou, Sian - sianhou1987@outlook.com
//          Wang, Guangchao - wgcupc@163.com

#include <mpi.h>
#include "../lib/sjinc.h"

int main(int argc, char *argv[]) {

    //! Runtime
    int flag = 1, is = 0;
    double tstart, tend, Tstart, Tend;

    //! MPI
    int mpiid, rankid, nrank;
    MPI_Status stauts;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    //! Survey
    sjssurvey sur;
    flag &= sjssurvey_init(&sur);
    flag &= sjssurvey_getparas(&sur, argc, argv);

    //! Geology
    sjsgeology geo;
    flag &= sjsgeo_init(&geo);
    flag &= sjsgeo_getparas2d(&geo, argc, argv, "vp");

    //! Wavefield
    sjswave wav;
    flag &= sjswave_init(&wav);
    flag &= sjswave_getparas(&wav, argc, argv, "profz");

    //! Option
    sjsoption opt;
    flag &= sjsoption_init(&opt);
    flag &= sjsoption_getparas(&opt, argc, argv);

    if (flag) {
        //! Time
        if (rankid == 0) {
            Tstart = (double) clock();
            printf("------------------------ 2D Acoustic simulation start  ------------------------\n");
        }

        //! Model
        geo.gvp2d = (float **) sjalloc2d(sur.gnx, sur.gnz, sizeof(float));
        sjreadsuall(geo.gvp2d[0], sur.gnx, sur.gnz, geo.vpfile);

        //! Simulation
        for (is = rankid; is < sur.ns; is += nrank) {
            //! Time
            tstart = (double) clock();

            //! Survey
            sjssurvey_readis(&sur, is);

            //! Model
            geo.vp2d = (float **) sjalloc2d(sur.nx, sur.nz, sizeof(float));
            sjextract2d(geo.vp2d, sur.x0, sur.z0, sur.nx, sur.nz, geo.gvp2d);

            //! Wavefield
            wav.profz = (float **) sjalloc2d(sur.nr, opt.nt, sizeof(float));
            wav.fwz2d = (float ***) sjalloc3d(opt.nsnap, sur.nx, sur.nz, sizeof(float));

            //! Simulation
            sjafor2d(&sur, &geo, &wav, &opt);

            //! Output
            if (rankid == 0) {
                sjwritesu(wav.profz[0], sur.nr, opt.nt, sizeof(float), opt.dt, is, wav.profzfile);
                tend = (double) clock();
                printf("Single shot simulation complete - %d/%d - time=%fs.\n", is + 1, sur.ns,
                       (tend - tstart) / CLOCKS_PER_SEC);
                printf("Rankid=%d, sx=%d, sz=%d, rx=%d to %d, rz=%d to %d.\n", rankid,
                       sur.sx + sur.x0, sur.sz + sur.z0,
                       sur.rx[0] + sur.x0, sur.rx[sur.nr - 1] + sur.x0,
                       sur.rz[0] + sur.z0, sur.rz[sur.nr - 1] + sur.z0);

                for (mpiid = 1; mpiid < nrank; mpiid++)
                    if ((is + mpiid) < sur.ns) {
                        MPI_Recv(wav.profz[0], sur.nr * opt.nt, MPI_FLOAT, mpiid, 99, MPI_COMM_WORLD, &stauts);
                        //! Output in rank != 0
                        sjwritesu(wav.profz[0], sur.nr, opt.nt, sizeof(float), opt.dt, is + mpiid, wav.profzfile);
                    }
            } else {
                MPI_Send(wav.profz[0], sur.nr * opt.nt, MPI_FLOAT, 0, 99, MPI_COMM_WORLD);
                tend = (double) clock();
                printf("Single shot simulation complete - %d/%d - time=%fs.\n", is + 1, sur.ns,
                       (tend - tstart) / CLOCKS_PER_SEC);
                printf("Rankid=%d, sx=%d, sz=%d, rx=%d to %d, rz=%d to %d.\n", rankid,
                       sur.sx + sur.x0, sur.sz + sur.z0,
                       sur.rx[0] + sur.x0, sur.rx[sur.nr - 1] + sur.x0,
                       sur.rz[0] + sur.z0, sur.rz[sur.nr - 1] + sur.z0);
            }

            //! Free
            sjmfree2d(geo.vp2d);
            sjmfree2d(wav.profz);
            sjmfree3d(wav.fwz2d);
        }

        //! Free Global Memory
        sjmfree2d(geo.gvp2d);

        //------------------------ Information ------------------------//
        if (rankid == 0) {
            Tend = (double) clock();
            printf("2D acoustic simulation complete - time=%fs.\n\n", (Tend - Tstart) / CLOCKS_PER_SEC);
        }

    } else {
        printf("\nExamples:   sjmpiawfd2d survey=survey.su vp=vp.su profz=profz.su nt=3001\n");
        sjbasicinformation();
    }


    MPI_Finalize();

    return 0;
}

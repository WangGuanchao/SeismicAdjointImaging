//
// Created by hsa on 09/12/16.
//

#include "../lib/sjinc.h"

int main(int argc, char *argv[]) {
    if (argc == 1) {
        printf("\nCreative a uniform 2d survey for adjoint imaging.\n\n");
        printf("ns:      Number of shots in survey.\n");
        printf("nr:      Number of receivers in each survey.\n");

        printf("vel:     Input a velocity file.\n");
        printf("x0:      Begin of local model in inline(x or n2) direction.\n");
        printf("nx:      Length of local model  in inline(x or n2) direction.\n");
        printf("dx0:     Interval of local model in inline(x or n2) direction.\n");

        printf("sx0:     Begin of shot (in local model) in inline(x or n2) direction.\n");
        printf("sz0:     Begin of shot (in local model) in vertical(z or n1) direction.\n");
        printf("dsx:     Interval of each shot (in local model) in inline(x or n2) direction.\n");

        printf("rx0:     Begin of receiver in inline(x or n2) direction.\n");
        printf("rz0:     Begin of receiver in vertical(z or n1) direction.\n");

        printf("drx:     Interval of receiver in inline(x or n2) direction in each shot.\n");
        printf("drz:     Interval of receiver in vertical(z or n1) direction in each shot.\n");


        printf("survey:  Output name of survey file.\n\n");
        printf("Explain: for (is = 0; is < ns; ++is)\n");
        printf("             x0 = x0 + is*dx0\n");
        printf("             z0 = z0\n");
        printf("             sx = sx0 + is*dsx\n");
        printf("             sz = sz0\n");
        printf("             for (ir = 0; ir < nr; ++ir)\n");
        printf("                 rx = rx0 + ir*drx\n");
        printf("                 rz = rz0 + ir*drz\n\n");
        printf("Example: sjsurvey2d vel=vel.su sx0=11 sz0=5 ns=2 dsx=10 dsz=0\n");
        printf("                    rx0=51 rz0=5 nr=5 drx=2 drz=0 x0=0 lxl=201\n");
        printf("                    survey=svy.su\n");
        sjbasicinformation();
    }
    //! Read parameters

    //! Number
    int ns, nr;
    if (!sjmgeti("ns", ns)) {
        printf("ERROR: Should set ns in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("nr", nr)) {
        printf("ERROR: Should set nr in program sjsurvey2d!\n");
        exit(0);
    };

    //! Model
    int x0, nx, dx0;
    char *inputname;
    if (!sjmgets("vel", inputname)) {
        printf("ERROR: Should set model in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("x0", x0)) {
        printf("ERROR: Should set x0 in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("nx", nx)) {
        printf("ERROR: Should set nx in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("dx0", dx0)) {
        printf("ERROR: Should set dx0 in program sjsurvey2d!\n");
        exit(0);
    };

    //! Source
    int sx0, sz0, dsx;
    if (!sjmgeti("sx0", sx0)) {
        printf("ERROR: Should set sx0 in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("sz0", sz0)) {
        printf("ERROR: Should set sz0 in program sjsurvey2d!\n");
        exit(0);
    };

    if (!sjmgeti("dsx", dsx)) {
        printf("ERROR: Should set dsx in program sjsurvey2d!\n");
        exit(0);
    };

    //! Receiver
    int rx0, rz0, drx, drz;
    if (!sjmgeti("rx0", rx0)) {
        printf("ERROR: Should set rx0 in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("rz0", rz0)) {
        printf("ERROR: Should set rz0 in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("drx", drx)) {
        printf("ERROR: Should set drx in program sjsurvey2d!\n");
        exit(0);
    };
    if (!sjmgeti("drz", drz)) {
        printf("ERROR: Should set drz in program sjsurvey2d!\n");
        exit(0);
    };

    //! Survey
    char *outputname;
    if (!sjmgets("survey", outputname)) {
        printf("ERROR: Should set survey in program sjsurvey2d!\n");
        exit(0);
    };

    int is, ir;
    //! Get Size of model
    int vpnx, vpnz;
    vpnx = sjgetsun2(sizeof(float), inputname);
    vpnz = sjgetsun1(sizeof(float), inputname);

    //! Check parameters
    for (is = 0; is < ns; ++is) {
        //! Local model
        if ((x0 + is * dx0) < 0 || (x0 + is * dx0) > (vpnx - 1)) {
            printf("ERROR: The begin of local model exceed the velocity model!\n");
            exit(0);
        }
        if ((x0 + is * dx0 + nx - 1) < 0 || (x0 + is * dx0 + nx - 1) > (vpnx - 1)) {
            printf("ERROR: The end of local model exceed the velocity model!\n");
            exit(0);
        }
        //! Source
        if ((sx0 + is*dsx) < 0 || (sx0 + is*dsx) > (nx - 1)) {
            printf("ERROR: The sx exceed the local model!\n");
            exit(0);
        }
        if (sz0 < 0 || sz0 > (vpnz - 1)) {
            printf("ERROR: The sz exceed the local model!\n");
            exit(0);
        }
        //! Receiver
        if (rx0 < 0 || rx0 > (nx - 1)) {
            printf("ERROR: The rx exceed the local model!\n");
            exit(0);
        }
        if ((rx0 + (nr-1) * drx) < 0 || (rx0 + (nr-1) * drx) > (nx - 1)) {
            printf("ERROR: The rx exceed the local model!\n");
            exit(0);
        }
        if (rz0 < 0 || rz0 > (vpnz - 1)) {
            printf("ERROR: The rz exceed the local model!\n");
            exit(0);
        }
        if ((rz0 + (nr - 1) * drz) < 0 || (rz0 + (nr - 1) * drz) > (vpnz - 1)) {
            printf("ERROR: The rz exceed the local model!\n");
            exit(0);
        }
    }

    //! Output survey
    sjssurvey survey;
    sjssurvey_init(&survey);
    survey.surveyfile = outputname;
    survey.maxnr = nr;
    survey.ry = (int *) sjalloc1d(nr, sizeof(int));
    survey.rx = (int *) sjalloc1d(nr, sizeof(int));
    survey.rz = (int *) sjalloc1d(nr, sizeof(int));
    for (is = 0; is < ns; ++is) {
        //! y0,x0,z0
        survey.y0 = 0;
        survey.x0 = x0 + is * dx0;
        survey.z0 = 0;
        //! sy,sx,sz
        survey.sy = 0;
        survey.sx = sx0 + is*dsx;
        survey.sz = sz0;
        //! yl,xl,zl
        survey.ny = 1;
        survey.nx = nx;
        survey.nz = vpnz;
        //! gyl,gxl,gzl
        survey.gny = 1;
        survey.gnx = vpnx;
        survey.gnz = vpnz;
        //! nr
        survey.nr = nr;
        //! tr
        if (is != 0) {
            survey.tr += nr;
        } else {
            survey.tr = 0;
        }
        //! receiver
        for (ir = 0; ir < nr; ++ir) {
            survey.ry[ir] = 0;
            survey.rx[ir] = rx0 + ir * drx;
            survey.rz[ir] = rz0 + ir * drz;
        }
        //! Output
        sjssurvey_write((void *) &survey, is);
    }

    return 0;
}
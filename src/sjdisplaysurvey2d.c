
//
// Created by hsa on 16/12/16.
//

#include "../lib/sjinc.h"
#include "../lib/sjfile.h"

int main(int argc, char *argv[]) {

    sjssurvey survey;
    sjssurvey_init(&survey);

    if (argc == 1) {
        printf("\nDisplay 2D survey.\n\n");
        sjssurvey_getparas(&survey, argc, argv);
        printf("  is:         The sequence number of shot to display, default = 1.\n");
        printf("              1<= is <= ns.\n\n");
        printf("  flag:       1 - Display local information, default = 1.\n");
        printf("  Example:    sjdisplaysurvey2d survey=svy.su\n");
        sjbasicinformation();
    }

    int is, ir, flag;

    //! Get parameters
    if (!sjmgeti("is", is)) is = 1;
    if (!sjmgeti("flag", flag)) flag = 1;

    sjssurvey_getparas(&survey, argc, argv);

    //! Check
    if (is > survey.ns) {
        printf("ERROR: Shot sequence number exceed ns!\n");
        exit(0);
    }

    sjssurvey_readis(&survey, is - 1);
    if (flag) {
        printf("\n");
        printf("Display local survey information(in local coordinate)\n");
        printf("\n");
        printf("Total shot number: %d\n", survey.ns);
        printf("Single shot sequence: %d\n", is);
        printf("Single shot position: sx=%d sz=%d\n", survey.sx, survey.sz);
        printf("\n");
        printf("Local model size: nx=%d nz=%d\n", survey.nx, survey.nz);
        printf("Local model begin: x0=%d z0=%d\n", survey.x0, survey.z0);
        printf("Global model size: gnx=%d gnz=%d\n", survey.gnx, survey.gnz);
        printf("\n");
        printf("Receiver number: nr=%d\n", survey.nr);
        printf("Receiver position (in local): rx=%-4d ", survey.rx[0]);
        for (ir = 1; ir < 3; ++ir) {
            if (ir < survey.nr)
                if (ir != survey.nr - 1)
                    printf("%-4d ", survey.rx[ir]);
        }
        printf("... %-4d\n", survey.rx[survey.nr - 1]);
        printf("Receiver position (in local): rz=%-4d ", survey.rz[0]);
        for (ir = 1; ir < 3; ++ir) {
            if (ir < survey.nr)
                if (ir != survey.nr - 1)
                    printf("%-4d ", survey.rz[ir]);
        }
        printf("... %-4d\n", survey.rz[survey.nr - 1]);
    } else {
        printf("\n");
        printf("Display global survey information(in global coordinate)");
        printf("\n");
        printf("Total shot number: %d\n", survey.ns);
        printf("Single shot sequence: %d\n", is);
        printf("Single shot position: sx=%d sz=%d\n", survey.sx + survey.x0, survey.sz + survey.z0);
        printf("\n");
        printf("Local model size: nx=%d nz=%d\n", survey.nx, survey.nz);
        printf("Local model begin: x0=%d z0=%d\n", survey.x0, survey.z0);
        printf("Global model size: gnx=%d gnz=%d\n", survey.gnx, survey.gnz);
        printf("\n");
        printf("Receiver number: nr=%d\n", survey.nr);
        printf("Receiver position (in global): rx=%-4d ", survey.rx[0] + survey.x0);
        for (ir = 1; ir < 3; ++ir) {
            if (ir < survey.nr)
                if (ir != survey.nr - 1)
                    printf("%-4d ", survey.rx[ir] + survey.x0);
        }
        printf("... %-4d\n", survey.rx[survey.nr - 1] + survey.x0);
        printf("Receiver position (in global): rz=%-4d ", survey.rz[0] + survey.z0);
        for (ir = 1; ir < 3; ++ir) {
            if (ir < survey.nr)
                if (ir != survey.nr - 1)
                    printf("%-4d ", survey.rz[ir] + survey.z0);
        }
        printf("... %-4d\n", survey.rz[survey.nr - 1] + survey.z0);
    }

    return 0;
}
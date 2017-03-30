//
// Created by hsa on 07/12/16.
//

#include "../lib/sjinc.h"

int main(int argc, char *argv[]) {

    //!  data

    float **d1 = (float **) sjalloc2d(460, 201, sizeof(float));
    float **d2 = (float **) sjalloc2d(460 + 50, 201 + 50, sizeof(float));
    sjreadsuall(d1[0],460,201,"/home/hsa/ClionProjects/SeismicAdjointImaging/vp.segy");
    sjextend2d(d2,460,201,25,25,25,25,d1);
    sjwritesuall(d2[0],460+50,201+50,10.0,"/home/hsa/ClionProjects/SeismicAdjointImaging/vp2.segy");
    sjreadsuall(d1[0],460,201,"/home/hsa/ClionProjects/SeismicAdjointImaging/rho.segy");
    sjextend2d(d2,460,201,25,25,25,25,d1);
    sjwritesuall(d2[0],460+50,201+50,10.0,"/home/hsa/ClionProjects/SeismicAdjointImaging/rho2.segy");

    return 0;
}

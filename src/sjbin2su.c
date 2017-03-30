//
// Created by hsa on 15/12/16.
//

#include "../lib/sjinc.h"

int main(int argc, char *argv[]) {

    if (argc == 1) {
        printf("\nConvert binary file to SU file.\n\n");
        printf("  E-mail:     sianhou1987@outlook.com\n");

        printf("* binary:     Input filename of binary file.\n");
        printf("* n2:         Number of samples in n2.\n");
        printf("* n1:         Number of samples in n1.\n");
        printf("  d1:         n1 interval, default = 0.001.\n");
        printf("* su:         Output filename of SU file.\n");
        printf("Example:      sjbin2su binary=dat.bin n2=301 n1=201 su=dat.su\n");
        sjbasicinformation();
    } else {
        //! Define parameters
        char *inputfile, *outputfile;
        int n2, n1;
        float d1;
        //! Read parameters
        if (!sjmgets("binary", inputfile)) {
            printf("ERROR: Should input binary file in program sjbin2su!\n");
            exit(0);
        }
        if (!sjmgets("su", outputfile)) {
            printf("ERROR: Should output SU file in program sjbin2su!\n");
            exit(0);
        }
        if (!sjmgeti("n2", n2)) {
            printf("ERROR: Should input n2 in program sjbin2su!\n");
            exit(0);
        }
        if (!sjmgeti("n1", n1)) {
            printf("ERROR: Should input n1 in program sjbin2su!\n");
            exit(0);
        }
        if (!sjmgetf("d1", d1)) d1 = 0.001f;

        FILE *fp;
        //! Allocate memory
        float **ptr = (float **) sjalloc2d(n2, n1, sizeof(float));
        //! Read binary
        if ((fp = fopen(inputfile, "rb")) != NULL) {
            //! Read data
            fread(ptr[0], sizeof(float), n2 * n1, fp);
            //! Write data
            sjwritesuall(ptr[0], n2, n1, d1, outputfile);
        } else {
            printf("ERROR: Cannot open file in program sjbin2su.\n");
            exit(0);
        }
        fclose(fp);
        sjmfree2d(ptr);
    };
}
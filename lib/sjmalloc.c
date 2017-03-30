// Author: Hou, Sian - sianhou1987@outlook.com

#include "sjmalloc.h"

//! function for allocate the 1 demension array and initialize with zero
void *sjalloc1d(int n1, int size) {
    void *p;

    if ((p = (void *) calloc(n1, size)) == NULL) {
        printf("ERROR:Can not allocate the 1d memory!\n");
        return NULL;
    }

    return p;
}

//! function for allocate the 2 demension array and initialize with zero
void **sjalloc2d(int n2, int n1, int size) {
    void **p;
    int ii;

    if ((p = (void **) calloc(n2, sizeof(void *))) == NULL) {
        printf("ERROR:Can not allocate the 2d memory\n");
        return NULL;
    }

    if ((p[0] = (void *) calloc(n2 * n1, size)) == NULL) {
        printf("ERROR:Can not allocate the 2d memory\n");
        return NULL;
    }

    for (ii = 0; ii < n2; ii++)
        p[ii] = (char *) p[0] + ii * n1 * size;
    // Note:
    // In ANCI C, the void * tpye pointer can not be used for operating.
    // So p[ii]=p[0]+ii*nj*size is wrong, p[ii]=(char *)p[0]+ii*nj*size is right
    //
    // IN GUN C, the void * tpye pointer can be used for operating, just like the char * type pointer.
    // So p[ii]=p[0]+ii*nj*size is right
    return p;
}

//! function for allocate the 3 demension array and initialize with zero
void ***sjalloc3d(int n3, int n2, int n1, int size) {
    void ***p;
    int ii, jj;

    //allocate pointer pointer, n3 size
    if ((p = (void ***) calloc(n3, sizeof(void **))) == NULL) {
        printf("ERROR: Can not allocate the 3d memory\n");
        return NULL;
    }

    //allocate pointer, n3*n2 size
    if ((p[0] = (void **) calloc(n3 * n2, sizeof(void *))) == NULL) {
        printf("ERROR: Can not allocate the 3d memory\n");
        return NULL;
    }

    //allocate data, n3*n2*n1 size
    if ((p[0][0] = (void *) calloc(n3 * n2 * n1, size)) == NULL) {
        printf("ERROR: Can not allocate the 3d memory\n");
        return NULL;
    }

    //pointer management
    for (ii = 0; ii < n3; ii++) {
        p[ii] = p[0] + ii * n2;
        for (jj = 0; jj < n2; jj++) {
            p[ii][jj] = (char *) p[0][0] + (ii * n2 + jj) * n1 * size;
        }
    }
    // Note:
    // In ANCI C, the void * tpye pointer can not be used for operating.
    // So p[ii]=p[0]+ii*nj*size is wrong, p[ii]=(char *)p[0]+ii*nj*size is right
    //
    // IN GUN C, the void * tpye pointer can be used for operating, just like the char * type pointer.
    // So p[ii]=p[0]+ii*nj*size is right
    return p;
}

//! function for free the 1 demension array
void sjfree1d(void *p) {
    free(p);
}

void sjcheckfree1d(void *p) {
    if (p != NULL) {
        free(p);
    }
}

//! function for free the 2 demension array
void sjfree2d(void **p) {
    free(p[0]);
    free(p);
}

void sjcheckfree2d(void **p) {
    if (p != NULL) {
        free(p[0]);
        free(p);
    }
}

//! function for free the 3 demension array
void sjfree3d(void ***p) {
    free(p[0][0]);
    free(p[0]);
    free(p);
}

void sjcheckfree3d(void ***p) {
    if (p != NULL) {
        free(p[0][0]);
        free(p[0]);
        free(p);
    }
}

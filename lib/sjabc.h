// Author: Hou, Sian - sianhou1987@outlook.com

#ifndef SJI_SJABC_H
#define SJI_SJABC_H

#include "sjinc.h"

//! Initialize two order time hybrid abc for 1 demension simulation
int sjinitthabc1d(float *vp, float *vs, float ds, float dt, int nxb, float *gxl, float *gxr);

//! Apply two order time hybrid abc for 1 demension simulation
int sjapplythabc1d(float *fp, float *cp, float *pp, float *gxl, float *gxr, int nxb, int nb, int marg);

//! Initialize one order time hybrid abc for 2 demension simulation
int sjinitohabc2d(float **vp, float **vs, float ds, float dt, int nxb, int nzb,
                  float **gxl, float **gxr, float **gzu, float **gzb);

//! Apply one order time hybrid abc for 2 demension simulation
int sjapplyohabc2d(float **fp, float **cp, float **gxl, float **gxr, float **gzu, float **gzb,
                   int nxb, int nzb, int nb, int marg);

//! Initialize two order time hybrid abc for 2 demension simulation
int sjinitthabc2d(float **vp, float **vs, float ds, float dt, int nxb, int nzb,
                  float **gxl, float **gxr, float **gzu, float **gzb);

//! Apply two order time hybrid abc for 2 demension simulation
int sjapplythabc2d(float **fp, float **cp, float **pp, float **gxl, float **gxr, float **gzu, float **gzb,
                   int nxb, int nzb, int nb, int marg);

#endif //SJI_SJABC_H

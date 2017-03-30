// Author: Hou, Sian - sianhou1987@outlook.com

#ifndef SJI_SJWAVE_H
#define SJI_SJWAVE_H

#include "sjinc.h"

int sjricker1d(float *ricker, int nt, int t0, float dt, float fp, float amp);

void sjextend2d(float **z, int nx, int nz, int ex0, int ex1, int ez0, int ez1, float **x);

void sjextract2d(float **z, int x0, int z0, int nx, int nz, float **x);

void sjfilter2dx(float **a, int n2, int n1, char *mode);

void sjsetsurface(float **a, int n2, int n1, float val);

/**********************************************************************************************/
/* ! Acoustic                                                                                 */
/**********************************************************************************************/

//! Two dimension constant density acoustic forward exploration
void sjafor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension constant density acoustic scatter forward exploration
void sjasfor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension constant density acoustic RTM backward exploration
void sjartmbac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension constant density acoustic Time-Shift RTM backward exploration
void sjatsrtmbac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension constant density acoustic FWI backward exploration
void sjafwibac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension constant density acoustic LSRTM forward simulation
void sjalsrtmfor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension constant density acoustic LSRTM backward exploration
void sjalsrtmbac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension constant density acoustic RTI backward exploration
void sjartibac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension constant density acoustic WTI backward exploration
void sjawtibac2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);


/**********************************************************************************************/
/* ! Elastic                                                                                  */
/**********************************************************************************************/

//! Two dimension constant density elastic forward exploration
void sjefor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension constant density elastic forward exploration with SGFD
void sjesgfor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);

//! Two dimension constant density elastic forward exploration with V-S SGFD
void sjevssgfor2d(sjssurvey *sur, sjsgeology *geo, sjswave *wav, sjsoption *opt);
#endif //SJI_SJWAVE_H

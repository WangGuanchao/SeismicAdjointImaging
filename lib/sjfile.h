// Author: Hou, Sian - sianhou1987@outlook.com

#ifndef SJI_SJFILE_H
#define SJI_SJFILE_H

#include "sjinc.h"

//! Show information
void sjbasicinformation();

//! Process standard input
int sjgetint(int argc, char **argv, char *option, int *value);

int sjgetfloat(int argc, char **argv, char *option, float *value);

int sjgetchar(int argc, char **argv, char *option, char *value);

int sjgetstring(int argc, char **argv, char *option, char **value);

#define sjmgeti(name, value) sjgetint(argc, argv, name, &value)

#define sjmgetf(name, value) sjgetfloat(argc, argv, name, &value)

#define sjmgetc(name, value) sjgetchar(argc, argv, name, &value)

#define sjmgets(name, value) sjgetstring(argc, argv, name, &value)

//! Process file

long int sggetfilesize(void *filename);

//! Read and write sample SU file

int sjgetsun1(size_t size, char *inputname);

int sjgetsun2(size_t size, char *inputname);

int sjwritesu(void *ptr, size_t n2, size_t n1, size_t size, float d1, int ifappend, char *outputname);

int sjreadsu(void *ptr, size_t n2, size_t n1, size_t size, size_t offsetn2, size_t offsetn1, char *inputname);

int sjwritesuall(float *ptr, int n2, int n1, float dt, char *outputname);

int sjreadsuall(float *ptr, int n2, int n1, char *inputname);

//! Survey
int sjssurvey_init(sjssurvey *ptr);

int sjssurvey_display(sjssurvey *ptr);

int sjssurvey_readis(sjssurvey *ptr, int is);

int sjssurvey_write(sjssurvey *ptr, int ifappend);

int sjssurvey_getparas(sjssurvey *ptr, int argc, char **argv);

//! model
int sjsgeo_init(sjsgeology *ptr);

int sjsgeo_display(sjsgeology *ptr);

int sjsgeo_getparas2d(sjsgeology *ptr, int argc, char **argv, char *info);

//! wave
int sjswave_init(sjswave *ptr);

int sjswave_display(sjswave *ptr);

int sjswave_getparas(sjswave *ptr, int argc, char **argv, char *info);

//! Option
int sjsoption_init(sjsoption *ptr);

int sjsoption_display(sjsoption *ptr);

int sjsoption_getparas(sjsoption *ptr, int argc, char **argv);

#endif //SJI_SJFILE_H

#! /bin/bash

mpdboot

binpath=../../../bin
outpath=.
nthd=64

#-----------------------------------------------------------------------------
# Convert the data
#-----------------------------------------------------------------------------

# $binpath/sjbin2su binary=sb2as0.rsf@ n2=2067 n1=401 su=vel0.su
# $binpath/sjbin2su binary=sb2as7r10.rsf@ n2=2067 n1=401 su=vel1.su

#-----------------------------------------------------------------------------
# Creative a survey
#-----------------------------------------------------------------------------

# $binpath/sjsurvey2d ns=108 nr=301 vel=vel0.su x0=341 nx=321 dx0=10 sx0=160 sz0=5 dsx=0 rx0=10 rz0=5 drx=1 drz=0 survey=survey.su

#-----------------------------------------------------------------------------
# Simulation
#-----------------------------------------------------------------------------

# mpirun -np $nthd $binpath/sjmpiawfd2d survey=survey.su vp=vel0.su recz=recz.su nt=4201 dt=0.001

#-----------------------------------------------------------------------------
# RTM
#-----------------------------------------------------------------------------

# mpirun -np $nthd $binpath/sjmpiartm2d survey=survey.su vp=vel0.su recz=recz.su nt=4201 dt=0.001 ipp=mig.su

#-----------------------------------------------------------------------------
# LSRTM
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpilsartm2d survey=$outpath/survey.su vp=$outpath/vel0.su recz=$outpath/recz.su dt=0.001 nt=4201 niter=50 ydetails=1 izz=$outpath/lsmig.su



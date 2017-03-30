#! /bin/bash

binpath=../../../bin
outpath=.
nthd=30

#-----------------------------------------------------------------------------
# Convert binary
#-----------------------------------------------------------------------------

$binpath/sjbin2su binary=vp.bin n2=301 n1=301 su=$outpath/vp.su

#-----------------------------------------------------------------------------
# Creative survey
#-----------------------------------------------------------------------------

$binpath/sjsurvey2d ns=30 nr=291 vel=$outpath/vp.su x0=0 nx=301 dx0=0 sx0=5 sz0=5 dsx=10 rx0=5 rz0=5 drx=1 drz=0 survey=$outpath/survey.su

#-----------------------------------------------------------------------------
# Simulation
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiawfd2d survey=$outpath/survey.su vp=$outpath/vp.su profz=$outpath/recz.su nt=1251 dt=0.002

#-----------------------------------------------------------------------------
# RTM
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiartm2d survey=$outpath/survey.su vp=$outpath/vp.su profz=$outpath/recz.su nt=1251 dt=0.002 izz=$outpath/mig.su

#-----------------------------------------------------------------------------
# LSRTM
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpilsartm2d survey=$outpath/survey.su vp=$outpath/vp.su profz=$outpath/recz.su nt=1251 dt=0.002 niter=50 ydetails=1 izz=$outpath/lsmig.su


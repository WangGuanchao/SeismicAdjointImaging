#! /bin/bash

binpath=../../../bin
outpath=.
nthd=37

#-----------------------------------------------------------------------------
# Convert binary
#-----------------------------------------------------------------------------

$binpath/sjbin2su binary=vp.bin n2=460 n1=201 su=$outpath/vp.su

#-----------------------------------------------------------------------------
# Creative survey
#-----------------------------------------------------------------------------

$binpath/sjsurvey2d ns=37 nr=441 vel=$outpath/vp.su x0=0 nx=460 dx0=0 sx0=50 sz0=5 dsx=10 rx0=10 rz0=5 drx=1 drz=0  survey=$outpath/survey.su

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

mpirun -np $nthd $binpath/sjmpialsrtm2d survey=$outpath/survey.su vp=$outpath/vp.su profz=$outpath/recz.su nt=1251 dt=0.002 niter=50 ydetails=1 izz=$outpath/lsmig.su


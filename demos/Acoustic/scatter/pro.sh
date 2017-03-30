#! /bin/bash

binpath=../../../bin
outpath=.
nthd=41

#-----------------------------------------------------------------------------
# Convert binary
#-----------------------------------------------------------------------------

$binpath/sjbin2su binary=vp.bin n2=501 n1=184 su=$outpath/vp.su

#-----------------------------------------------------------------------------
# Creative survey
#-----------------------------------------------------------------------------

$binpath/sjsurvey2d ns=41 nr=482 vel=$outpath/vp.su x0=0 nx=501 dx0=0 sx0=50 sz0=5 dsx=10 rx0=10 rz0=5 drx=1 drz=0  survey=$outpath/survey.su

#-----------------------------------------------------------------------------
# Simulation
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiawfd2d survey=$outpath/survey.su vp=$outpath/vp.su profz=$outpath/recz.su nt=2501 dt=0.001

#-----------------------------------------------------------------------------
# RTM
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiartm2d survey=$outpath/survey.su vp=$outpath/vp.su profz=$outpath/recz.su nt=2501 dt=0.001 izz=$outpath/mig.su

#-----------------------------------------------------------------------------
# LSRTM
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpilsartm2d survey=$outpath/survey.su vp=$outpath/vp.su profz=$outpath/recz.su nt=2501 dt=0.001 niter=50 ydetails=1 izz=$outpath/lsmig.su


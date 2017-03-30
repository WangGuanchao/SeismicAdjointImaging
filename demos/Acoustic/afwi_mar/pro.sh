#! /bin/bash

binpath=../../../bin
outpath=.
nthd=41

#-----------------------------------------------------------------------------
# Creative survey
#-----------------------------------------------------------------------------

$binpath/sjsurvey2d ns=41 nr=482 vel=$outpath/vp.su x0=0 nx=501 dx0=0 sx0=50 sz0=5 dsx=10 rx0=10 rz0=5 drx=1 drz=0  survey=$outpath/survey.su

#-----------------------------------------------------------------------------
# Simulation
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiawfd2d survey=$outpath/survey.su vp=$outpath/vp.su profz=$outpath/recz.su nt=2501 dt=0.001 fp=5

#-----------------------------------------------------------------------------
# AFWI
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiafwi2d survey=$outpath/survey.su vp=$outpath/inivp.su profz=$outpath/recz.su nt=2501 dt=0.001 fp=5 jsnap=1 niter=5 ydetails=1 izz=$outpath/final_result.su


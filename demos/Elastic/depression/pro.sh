#! /bin/bash

binpath=../../../bin
outpath=.
nthd=37

#-----------------------------------------------------------------------------
# Convert binary
#-----------------------------------------------------------------------------

$binpath/sjbin2su binary=vp.bin n2=460 n1=201 su=$outpath/vp.su
$binpath/sjbin2su binary=vs.bin n2=460 n1=201 su=$outpath/vs.su

#-----------------------------------------------------------------------------
# Creative survey
#-----------------------------------------------------------------------------

$binpath/sjsurvey2d ns=37 nr=441 vel=$outpath/vp.su x0=0 nx=460 dx0=0 sx0=50 sz0=5 dsx=10 rx0=10 rz0=15 drx=1 drz=0  survey=$outpath/survey.su

#-----------------------------------------------------------------------------
# Simulation
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiewfd2d survey=$outpath/survey.su vp=$outpath/vp.su vs=$outpath/vs.su profx=$outpath/recx.su profz=$outpath/recz.su nt=2000 dt=0.001

mpirun -np $nthd $binpath/sjmpiewsgfd2d survey=$outpath/survey.su vp=$outpath/vp.su vs=$outpath/vs.su profx=$outpath/recx_sg.su profz=$outpath/recz_sg.su nt=2000 dt=0.001

mpirun -np $nthd $binpath/sjmpiewvssgfd2d survey=$outpath/survey.su vp=$outpath/vp.su vs=$outpath/vs.su profx=$outpath/recx_vssg.su profz=$outpath/recz_vssg.su nt=2000 dt=0.001 ycutdirect=0
#! /bin/bash

binpath=../../../bin
outpath=.

nthd=37

#-----------------------------------------------------------------------------
# Convert binary
#-----------------------------------------------------------------------------

$binpath/sjbin2su binary=vp0.bin n2=460 n1=201 su=$outpath/vp0.su

$binpath/sjbin2su binary=vps.bin n2=460 n1=201 su=$outpath/vps.su

$binpath/sjbin2su binary=vpf.bin n2=460 n1=201 su=$outpath/vpf.su

#-----------------------------------------------------------------------------
# Creative survey
#-----------------------------------------------------------------------------

$binpath/sjsurvey2d ns=37 nr=441 vel=$outpath/vp0.su x0=0 nx=460 dx0=0 sx0=50 sz0=5 dsx=10 rx0=10 rz0=5 drx=1 drz=0 survey=$outpath/survey.su

#-----------------------------------------------------------------------------
# Simulation
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiawfd2d survey=$outpath/survey.su vp=$outpath/vp0.su profz=$outpath/recz.su nt=1500 k1=50 dt=0.002

#-----------------------------------------------------------------------------
# Inversion
#-----------------------------------------------------------------------------

mpirun -np $nthd $binpath/sjmpiarti2d survey=$outpath/survey.su vp=vps.su profz=$outpath/recz.su nt=1500 k1=50 dt=0.002 niter=100 ydetails=1 izz=$outpath/ipps.su

# mpirun -np $nthd $binpath/sjmpiarti2d survey=$outpath/survey.su vp=vp0.su profz=$outpath/recz.su nt=2001 k1=200 dt=0.002 ipp=$outpath/ipp0.su

# mpirun -np $nthd $binpath/sjmpiarti2d survey=$outpath/survey.su vp=vpf.su profz=$outpath/recz.su nt=2001 k1=200 dt=0.002 ipp=$outpath/ippf.su
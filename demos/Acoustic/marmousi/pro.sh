#! /bin/bash

#-----------------------------------------------------------------------------
# Creative a survey data
#-----------------------------------------------------------------------------

../../bin/sgsurvey2d sgin=marvel.su sx0=361 sz0=5 ns=2 dsx=10 dsz=0 rx0=211 rz0=5 nr=301 drx=1 drz=0 lx0=201 lxl=321 sgot=marsvy.su

#-----------------------------------------------------------------------------
# Simulation with GFDXYZ
#-----------------------------------------------------------------------------

# MPI
mpirun -np 2 ../../bin/sjmpiawsgfd2d svy=marsvy.su vp=marvel.su rec=marrec.su nt=3001 dt=0.001

#-----------------------------------------------------------------------------
# RTM with GFDXYZ
#-----------------------------------------------------------------------------

# MPI
mpirun -np 2 ../../bin/sjmpiartm2d svy=marsvy.su vp=marvel.su rec=marrec.su izz=marmig_mpi.su

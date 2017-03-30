#! /bin/bash

sfspike n1=131 n2=460 mag=2300.0 > vpu.rsf
sfspike n1=70 n2=460 mag=2300.0 > vpb.rsf
sfcat < vpu.rsf vpb.rsf axis=1 > vpf.rsf
cp /var/tmp/vpf.rsf@ vpf.bin

sfrm *.rsf

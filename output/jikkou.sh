#!/bin/sh
./compile2
time ./ficticious
sh gnu.sh
#ifort f_to_d.f90 -I/usr/local/include -lfftw3
#./a.out

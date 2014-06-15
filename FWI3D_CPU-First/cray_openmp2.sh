#!/bin/bash
#============ LSF Options ============
#QSUB -q gr10042a
#QSUB -W 72:00
#QSUB -A p=1:t=32:c=32:m=2G

#============ Shell Script ============
set -x

# automatically
export OMP_NUM_THREADS=$LSB_THREADS

#aprun -n $LSB_PROCS -d $LSB_CPUS -N $LSB_PPN ./fwi2d_posi.out
time aprun -n $LSB_PROCS -d $LSB_CPUS -N $LSB_PPN sh run2.sh
#aprun -n $LSB_PROCS -d $LSB_CPUS -N $LSB_PPN ./fwi2d_wave.out

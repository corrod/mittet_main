#!/bin/bash

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'hz(x0,y0,z0+10).eps'
p 'hz(x0,y0,z0+10).d' w l title '10 from shot'
EOF

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'hz(x0,y0,z0+20).eps' 
p 'hz(x0,y0,z0+20).d' w l title '20 from shot'

EOF

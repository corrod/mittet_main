#!/bin/bash

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'hz1010.eps'
p 'hz1010.d' w l title '10 from shot'
EOF

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'hz1020.eps' 
p 'hz1020.d' w l title '20 from shot'

EOF

#!/bin/bash

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'hz1010.eps'
p 'hz1010.d' every 5 u 1:2 w l title '10 from shot'
EOF

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'hz1020.eps' 
p 'hz1020.d' every 5 u 1:2 w l title '20 from shot'
EOF

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'hz1030.eps' 
p 'hz1030.d' every 5 u 1:2 w l title '30 from shot'
EOF
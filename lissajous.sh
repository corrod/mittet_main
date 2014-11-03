gnuplot<<EOF
set term posts enhanced color
set out 'lissajous2_23.eps'
set xlabel 'nocrack'
set ylabel 'crack'
set title 'Lissajous 23'
plot "< paste code_ver0023grid_nocrack2/pattern2_23diff.d code_ver0023grid_crack2/pattern2_23diff.d" every ::1::4000 using 2:4 w l lw 2
EOF

gnuplot<<EOF
set term posts enhanced color
set out 'lissajous2_12.eps'
set xlabel 'nocrack'
set ylabel 'crack'
set title 'Lissajous 12'
plot "< paste code_ver0023grid_nocrack2/pattern2_12diff.d code_ver0023grid_crack2/pattern2_12diff.d" every ::1::4000 using 2:4 w l lw 2
EOF
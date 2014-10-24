パターン 2-3
# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'pattern2_differential.eps'
# p "< paste patern2_2.d patern2_3.d" using 1:($2-$4) w l
# EOF

gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'pattern2_differential_all.eps'
p "< paste patern2_2.d patern2_3.d" using 1:($2-$4) w l, "" using 1:2 w l, "" using 1:4 w l
EOF
#!/bin/sh

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'pattern2_12diff.eps'
p 'pattern2_12diff.d' w l lw 2
EOF

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'pattern2_23diff.eps'
p 'pattern2_23diff.d' w l lw 2
EOF

# パターン1 ①−②
gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'pattern1_12differential.eps'
p "< paste pattern1_1.d pattern1_2.d" using 1:(\$2-\$5) w l
EOF

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'pattern1_12differential_all.eps'
p "< paste pattern1_1.d pattern1_2.d" using 1:(\$2-\$5) w l, "" using 1:2 w l, "" using 1:5 w l
EOF

# パターン1 ②−③
gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'pattern1_23differential.eps'
p "< paste pattern1_2.d pattern1_3.d" using 1:(\$2-\$5) w l
EOF

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'pattern1_23differential_all.eps'
p "< paste pattern1_2.d pattern1_3.d" using 1:(\$2-\$5) w l, "" using 1:2 w l, "" using 1:5 w l
EOF


# # パターン２ ①ー②
# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'pattern2_12differential.eps'
# p "< paste pattern2_1.d pattern2_2.d" using 1:(\$2-\$5) w l
# EOF

# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'pattern2_12differential_all.eps'
# p "< paste pattern2_1.d pattern2_2.d" using 1:(\$2-\$5) w l, "" using 1:2 w l, "" using 1:5 w l
# EOF

# # パターン２ ②−③
# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'pattern2_23differential.eps'
# p "< paste pattern2_2.d pattern2_3.d" using 1:(\$2-\$5) w l
# EOF

# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'pattern2_23differential_all.eps'
# p "< paste pattern2_2.d pattern2_3.d" using 1:(\$2-\$5) w l, "" using 1:2 w l, "" using 1:5 w l
# EOF
#!/bin/bash
# HZの伝播確認
gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'hz_from_shot.eps'
p 'hz1010.d' every 10 u 1:2 w l title '10 from shot', \
'hz1020.d' every 10 u 1:2 w l title '20 from shot', \
'hz1030.d' every 10 u 1:2 w l title '30 from shot'
EOF

# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'hz1010.eps'
# p 'hz1010.d' every 10 u 1:2 w l title '10 from shot'
# EOF

# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'hz1020.eps'
# p 'hz1020.d' every 10 u 1:2 w l title '20 from shot'
# EOF

# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'hz1030.eps'
# p 'hz1030.d' every 10 u 1:2 w l title '30 from shot'
# EOF


# モデルの出力
gnuplot<<EOF
set term postscript enhanced color
set out 'model_sig_xz.eps'
set dgrid3d 101,101,4
set pm3d map
set size ratio 1
set palette rgbformulae 33,11,10
set xlabel "x"
set ylabel "z"
set cblabel "[S/m]"
set xrange[1:91]
set yrange[91:1]
splot 'model_sig_xz.dat' u 1:3:4 with pm3d
EOF

gnuplot<<EOF
set term postscript enhanced color
set out 'model_sig_xy.eps'
set dgrid3d 101,101, 4
set pm3d map
set samples 100, 100
set isosamples 10, 10
set mapping cartesian
set size ratio 1
set style data pm3d
set style function pm3d
set xlabel "x"
set xrange [ 1.00000 : 101.0000 ] noreverse nowriteback
set ylabel "y"
set yrange [ 1.00000 : 41.0000 ] noreverse nowriteback
#set zrange [ 3.20000 : 7.50000e+06 ] noreverse nowriteback
set cblabel "[S/m]"
set cbrange [ 3.20000 : 7.50000e+06 ] noreverse nowriteback
set palette rgbformulae 33, 13, 10
splot 'model_sig_xy.dat' u 1:2:4
EOF

gnuplot<<EOF
set term postscript enhanced color
set out 'model_myu_xz.eps'
set dgrid3d 101,101,4
set pm3d map
#set size square
set size ratio 1
set palette rgbformulae 33,11,10
set xlabel "x"
set ylabel "z"
set cblabel "[H/m]"
set xrange[1:101]
set yrange[41:1]
splot 'model_myu_xz.dat' u 1:3:4
EOF

gnuplot<<EOF
set term postscript enhanced color
set out 'model_myu_xy.eps'
set dgrid3d 101,101,4
set pm3d map
set size ratio 1
set samples 100, 100
set isosamples 10, 10
set mapping cartesian
unset hidden3d
set xlabel "x"
set xrange [ 1.00000 : 101.0000 ] noreverse nowriteback
set ylabel "y"
set yrange [ 1.00000 : 41.0000 ] noreverse nowriteback
set zrange [ 1.25000e-06 : 0.00510000 ] noreverse nowriteback
set cblabel "[H/m]"
set cbrange [ 1.25000e-06 : 0.00510000 ] noreverse nowriteback
set palette rgbformulae 33, 13, 10
splot "model_myu_xy.dat" u 1:2:4
EOF


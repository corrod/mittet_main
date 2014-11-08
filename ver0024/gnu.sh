#!/bin/bash
# HZの伝播確認

# real
# ficticious
gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'ex_from_shot.eps'
p 'ex1010.d' every 10 u 1:2 w l lw 2 title '10 from shot ex', \
'ex1020.d' every 10 u 1:2 w l lw 2 title '20 from shot ex', \
'ex1030.d' every 10 u 1:2 w l lw 2 title '30 from shot ex'
EOF

gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'hz_from_shot.eps'
p 'hz1010.d' every 10 u 1:2 w l lw 2 title '10 from shot hz',\
'hz1020.d' every 10 u 1:2 w l lw 2 title '20 from shot hz', \
'hz1030.d' every 10 u 1:2 w l lw 2 title '30 from shot hz'
EOF


# diffusive freqency
gnuplot<<EOF
set xlabel 'f[Hz]'
set ylabel 'amp'
set term postscript enhanced color
set out 'ex_wfrom_shot.eps'
p 'EX_w1010.d' every 10 u 1:2 w l lw 2 title '10 from shot ex_w', \
'EX_w1020.d' every 10 u 1:2 w l lw 2 title '20 from shot ex_w', \
'EX_w1030.d' every 10 u 1:2 w l lw 2 title '30 from shot ex_w'
EOF

gnuplot<<EOF
set xlabel 'f[Hz]'
set ylabel 'amp'
set term postscript enhanced color
set out 'hz_wfrom_shot.eps'
p 'HZ_w1010.d' every 10 u 1:2 w l lw 2 title '10 from shot hz_w', \
'HZ_w1020.d' every 10 u 1:2 w l lw 2 title '20 from shot hz_w', \
'HZ_w1030.d' every 10 u 1:2 w l lw 2 title '30 from shot hz_w'
EOF

# diffusive time
gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'ex_tfrom_shot.eps'
p 'EX_t1010.d' every 10 u 1:2 w l lw 2 title '10 from shot ex_t', \
'EX_t1020.d' every 10 u 1:2 w l lw 2 title '20 from shot ex_t', \
'EX_t1030.d' every 10 u 1:2 w l lw 2 title '30 from shot ex_t'
EOF

# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'hz_tfrom_shot.eps'
# p 'HZ_t1010.d' every 10 u 1:2 w l title '10 from shot hz_t',\
# 'HZ_t1020.d' every 10 u 1:2 w l title '20 from shot hz_t', \
# 'HZ_t1030.d' every 10 u 1:2 w l title '30 from shot hz_t'
# EOF




# imaginary
# ficticious  imaginary は ない
# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'eximag_from_shot.eps'
# p 'ex1010.d' every 10 u 1:3 w l title '10 from shot ex imag', \
# 'ex1020.d' every 10 u 1:3 w l title '20 from shot ex imag', \
# 'ex1030.d' every 10 u 1:3 w l title '30 from shot ex imag'
# EOF

# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'hzimag_from_shot.eps'
# p 'hz1010.d' every 10 u 1:3 w l title '10 from shot hz imag',\
# 'hz1020.d' every 10 u 1:3 w l title '20 from shot hz imag', \
# 'hz1030.d' every 10 u 1:3 w l title '30 from shot hz imag'
# EOF


# diffusive freqency
gnuplot<<EOF
set xlabel 'f[Hz]'
set ylabel 'amp'
set term postscript enhanced color
set out 'eximag_wfrom_shot.eps'
p 'EX_w1010.d' every 10 u 1:3 w l lw 2 title '10 from shot ex_w imag', \
'EX_w1020.d' every 10 u 1:3 w l lw 2 title '20 from shot ex_w imag', \
'EX_w1030.d' every 10 u 1:3 w l lw 2 title '30 from shot ex_w imag'
EOF

gnuplot<<EOF
set xlabel 'f[Hz]'
set ylabel 'amp'
set term postscript enhanced color
set out 'hzimag_wfrom_shot.eps'
p 'HZ_w1010.d' every 10 u 1:3 w l lw 2 title '10 from shot hz_w imag', \
'HZ_w1020.d' every 10 u 1:3 w l lw 2 title '20 from shot hz_w imag', \
'HZ_w1030.d' every 10 u 1:3 w l lw 2 title '30 from shot hz_w imag'
EOF

# diffusive time
gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'eximag_tfrom_shot.eps'
p 'EX_t1010.d' every 10 u 1:3 w l lw 2 title '10 from shot ex_t imag', \
'EX_t1020.d' every 10 u 1:3 w l lw 2 title '20 from shot ex_t imag', \
'EX_t1030.d' every 10 u 1:3 w l lw 2 title '30 from shot ex_t imag'
EOF

# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'hzimag_tfrom_shot.eps'
# p 'HZ_t1010.d' every 10 u 1:3 w l title '10 from shot hz_t imag',\
# 'HZ_t1020.d' every 10 u 1:3 w l title '20 from shot hz_t imag', \
# 'HZ_t1030.d' every 10 u 1:3 w l title '30 from shot hz_t imag'
# EOF



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
set xrange[1:101]
set yrange[1:101]
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
set xrange [ 1.00000 : 50.0000 ] noreverse nowriteback
set ylabel "y"
set yrange [ 1.00000 : 50.0000 ] noreverse nowriteback
#set zrange [ 3.20000 : 7.50000e+06 ] noreverse nowriteback
set cblabel "[S/m]"
# set cbrange [ 3.20000 : 7.50000e+06 ] noreverse nowriteback
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
set xrange[1:50]
set yrange[1:166]
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
set xrange [ 1.00000 : 50.0000 ] noreverse nowriteback
set ylabel "y"
set yrange [ 1.00000 : 50.0000 ] noreverse nowriteback
# set zrange [ 1.25000e-06 : 0.00510000 ] noreverse nowriteback
set cblabel "[H/m]"
# set cbrange [ 1.25000e-06 : 0.00510000 ] noreverse nowriteback
set palette rgbformulae 33, 13, 10
splot "model_myu_xy.dat" u 1:2:4
EOF


# 差分プローブ用出力
# パターン2 1-2
gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'pattern2_12diff.eps'
p '41pattern2_12diff.d' w l lw 2
EOF
# パターン2 2-3
gnuplot<<EOF
set xlabel 't[s]'
set ylabel 'amp'
set term postscript enhanced color
set out 'pattern2_23diff.eps'
p '41pattern2_23diff.d' w l lw 2
EOF

# # パターン1 ①−②
# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'pattern1_12differential.eps'
# p "< paste pattern1_1.d pattern1_2.d" using 1:(\$2-\$5) w l
# EOF

# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'pattern1_12differential_all.eps'
# p "< paste pattern1_1.d pattern1_2.d" using 1:(\$2-\$5) w l, "" using 1:2 w l, "" using 1:5 w l
# EOF

# # パターン1 ②−③
# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'pattern1_23differential.eps'
# p "< paste pattern1_2.d pattern1_3.d" using 1:(\$2-\$5) w l
# EOF

# gnuplot<<EOF
# set xlabel 't[s]'
# set ylabel 'amp'
# set term postscript enhanced color
# set out 'pattern1_23differential_all.eps'
# p "< paste pattern1_2.d pattern1_3.d" using 1:(\$2-\$5) w l, "" using 1:2 w l, "" using 1:5 w l
# EOF

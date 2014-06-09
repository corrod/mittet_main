#!/usr/bin/sh
i=10001
max=10005

while [ $i -le $max ] ; do
gnuplot<<EOF

set xlabel 'x[m]'
set ylabel 'y[m]'
set xrange [1:100]
set yrange [1:100]
set zrange [-1:1]
splot '$i.d' u 2:3:4 w l
set term postscript enhanced color
set out '$i.eps'
replot

EOF

i=`expr $i + 1`

echo "$i"

done
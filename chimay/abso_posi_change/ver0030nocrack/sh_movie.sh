#!/bin/sh
#HZの波動伝播movie
i=10001
j=`expr $i - 10000`
max=70

while [ $j -le $max ]; do
gnuplot <<EOF
 set terminal postscript color
 set out "$i.ps"
 set font "Helvetica,20"
 set dgrid3d 100,100
 set pm3d map
 set size square
 # set palette rgbformulae 33,13,10
 # set cbrange[-1e+15:1e+15]
 # set zr[-6.6e+13:6.6e+13]
 splot "./out1/hz$i.d" u 2:3:4
EOF
i=`expr $i + 5`
j=`expr $j + 5`
done

convert 1*.ps hzfield_xz.gif
rm 1*ps


#Ex
i=10001
j=`expr $i - 10000`
max=70

while [ $j -le $max ]; do
gnuplot <<EOF
 set terminal postscript color
 set out "$i.ps"
 set font "Helvetica,20"
 set dgrid3d 100,100
 set pm3d map
 set size square
 set palette rgbformulae 33,13,10
 # set cbrange[-9e+06:9e+06]
 # set zr[-9e+05:9e+05]
 splot "./out2/ex$i.d" u 2:3:4
EOF
i=`expr $i + 5`
j=`expr $j + 5`
done

convert 1*.ps exfield_xy.gif
rm 1*ps


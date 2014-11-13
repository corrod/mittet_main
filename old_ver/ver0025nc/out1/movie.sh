#!/bin/sh
#HZの波動伝播movie
i=10001
j=`expr $i - 10000`
max=80

while [ $j -le $max ]; do
gnuplot <<EOF
 set terminal postscript color
 set out "$i.ps"
 set font "Helvetica,20"
 set dgrid3d 100,100
 set pm3d map
 set size square
 set palette rgbformulae 33,13,10
 # set cbrange[-1e+11:1e+11]
 # set zr[-1e+11:1e+11]
 # set cbrange[-6e+11:6e+11]
 # set zr[-6e+11:6e+11]
 splot "hz$i.d" u 2:3:4
EOF
i=`expr $i + 1`
j=`expr $j + 1`
done

convert 1*.ps hzfield.gif
rm 1*ps


# #Ex
# i=10001
# j=`expr $i - 10000`
# max=10

# while [ $j -le $max ]; do
# gnuplot <<EOF
#  set terminal postscript color
#  set out "$i.ps"
#  set font "Helvetica,20"
#  set dgrid3d 100,100
#  set pm3d map
#  set size square
#  set palette rgbformulae 33,13,10
#  set cbrange[-9e+06:9e+06]
#  set zr[-9e+05:9e+05]
#  splot "ex$i.d" u 2:3:4
# EOF
# i=`expr $i + 1`
# j=`expr $j + 1`
# done

# convert 1*.ps exfield.gif
# rm 1*ps


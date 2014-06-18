#!/bin/sh
#HZの波動伝播movie
i=10001
j=`expr $i - 10000`
max=66

while [ $j -le $max ]; do
gnuplot <<EOF
 set terminal postscript color
 set out "$i.ps"
 set font "Helvetica,20"
 set dgrid3d 100,100
 set pm3d map
 set size square
 set palette rgbformulae 33,13,10
 set cbrange[-0.05:0.05]#[-0.0002:0.0002]
 set zr[-0.05:0.05]#[-0.0002:0.0002]
 splot "hz$i.d" u 2:3:4
EOF

i=`expr $i + 1`
j=`expr $j + 1`

done

convert 1*.ps hzfield.gif
rm 1*ps

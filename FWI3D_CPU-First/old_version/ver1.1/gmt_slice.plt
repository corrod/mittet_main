#!/usr/local/Cellar/gnuplot/4.6.0/bin/gnuplot

#set size 1.25
#set format z "%6.3e"
set terminal png
set xlabel "X axis"
set ylabel "E-field"
set xrange[0:500]
set yrange[-1.e-5:1.e-5]
set format y "%e"
#set ytics 0,150,300

if (exist("n")==0 || n<0) n=0 #変数の初期化

file0(n) = sprintf("./data2/time_slice_%05d.dat",n) #入力ファイル名
outfile(n) = sprintf("./data3/output_%05d.png",n)  #出力ファイル名
title(n) = sprintf("t = %d",n)  #タイトル名

#unset label 
#set label title(n)  font 'Times,20' at 1800.0,  400.0, 0.0001

set output outfile(n)
plot file0(n) u 1:4 w l

n=n+10
if (n<600); reread

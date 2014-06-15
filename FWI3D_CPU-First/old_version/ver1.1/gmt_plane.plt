#!/usr/local/Cellar/gnuplot/4.6.0/bin/gnuplot

set size 1.25
set format z "%6.3e"
set terminal png
set xlabel "X axis"
set ylabel "Z axis"
set zrange[-1.5e-10:1.5e-10]
set ytics 0,150,300

if (exist("n")==0 || n<0) n=0 #変数の初期化

file0(n) = sprintf("./data2/time_plane_%05d.dat",n) #入力ファイル名
outfile(n) = sprintf("./data3/output_%05d.png",n)  #出力ファイル名
title(n) = sprintf("t = %d",n)  #タイトル名

unset label 
set view 80,80
set view equal xy
set label title(n)  font 'Times,20' at 1800.0,  400.0, 0.0001

set output outfile(n)
splot file0(n) u 1:3:4

n=n+10
if (n<2000); reread

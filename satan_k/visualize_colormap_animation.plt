if(i==0)
    set size square;
    set pm3d;
    set pm3d map;
    set pm3d interpolate 5,5;
    set cbrange[-0.0001:0.0001];
#    set palette model HSV defined (0 0.7 1 1,1 0 1 1)
    set palette gray

plot_title = sprintf("3D FDTD step = %d",i)
set title plot_title
splot "fieldzx.dat" index i title "E" with pm3d

i=i+20
pause 0.01
if (i<900) reread


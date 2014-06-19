# gnuplot script for gif animation 
# It is necessary that visualize_colormap.plt exists in the same directory.

input_data="field.dat" # The name of input data file
data_title="vecabs(Ex,Ey)"
plot_title="3D FDTD"
steps = 400  # The number of steps of animation

set size square;
set pm3d;
set pm3d map;
set pm3d interpolate 5,5;
set cbrange[-0.0002:0.0002];
#set palette model HSV defined ( 0 0.7 1 1, 1 0 1 1 )
set palette gray

set terminal gif animate optimize size 960, 720
set output "visualize_colormap.gif"

i = 0
load "visualize_colormap_animation.plt"
i = 0

set output

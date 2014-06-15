#!/bin/bash

./gmt_plane.rb
./gmt_plane.gmt

convert $HOME/work/FWI3D_CPU/data3/output* $HOME/work/FWI3D_CPU/data3/movie.gif

rm $HOME/work/FWI3D_CPU/data3/output*

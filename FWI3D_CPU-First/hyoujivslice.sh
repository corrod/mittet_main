#!/bin/bash -f
# GMT script for making bathymetry map (Japan)

#gmtdefaults -D > .gmtdefaults
#gmtset LABEL_FONT_SIZE 12 HEADER_FONT_SIZE 18
#read num
set Page Orientation landscape
rm tmp*

a=-100
b=0
c=-380
d=380
e=-400
f=0

#####################################################################################################

I=1
hma=-4
hmi=-6

awk '{print  $3,$2, $1, $4 ,0.2,0.5}' < estk3.dat > tmpe.dat
awk '{print  $3,$2, $1 }' < estk3.dat > tmpe2.dat

makecpt -Csealand -Z  -Q -T$hmi/$hma/1.e-2 > tmp.cpt
psxyz tmpe.dat -R$a/$b/$c/$d/$e/$f -E45/25 -JX2/8  -Jz0.04 -Sr -Ctmp.cpt -Ba50f50g50:"Z[m]":/a200f200g200:"Y[m]":/a200f200g200:"X[m]":wsNEZ+ -Xa1.6i -Ya1.i   -K > tmp.ps
#psxyz tmpe2.dat -R$a/$b/$c/$d/$e/$f -E45/25 -JX2/8  -Jz0.04 -Sc0.1 -Gblack -Xa1.6i -Ya1.i   -K -O >> tmp.ps


a=-100
b=0
c=-380
d=380
e=0
f=400

awk '{print  $3,$2, $1, $4 ,0.2,0.5}' < estk3.dat > tmpe.dat
makecpt -Csealand -Z  -Q -T$hmi/$hma/1.e-2 > tmp.cpt
psxyz tmpe.dat -R$a/$b/$c/$d/$e/$f -E45/25 -JX2/8 -Jz0.04 -Sr -Ctmp.cpt -Ba50f50g50:"Z[m]":/a200f200g200:"Y[m]":/a200f200g200:"X[m]":wsNEZ+ -Xa5.6i -Ya1.i -P -K -O >> tmp.ps

#mogrify -rotate 90 tmp.ps 
#grdmask tmpe.dat -R$a/$b/$c/$d -I40/10 -NNaN/1/1 -Gmask.grd
#grdmath tmpe.grd mask.grd OR = tmp.grd
#grdimage tmp.grd  -JX7/3 -R$a/$b/$c/$d -E45/45 -Ctmp.cpt -Xa1.i -Ya6.7i -P -K > tmp.ps
psscale -Ctmp.cpt -Q -D21.9/10.8/4.0/0.2 -B5/::  -O >> tmp.ps

gv tmp.ps
mv tpm.ps model.ps

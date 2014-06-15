#!/bin/bash -f
# GMT script for making bathymetry map (Japan)

#gmtdefaults -D > .gmtdefaults
#gmtset LABEL_FONT_SIZE 12 HEADER_FONT_SIZE 18
#read num
rm tmp*

file="data"   # directory including plane data
name="sig_000.dat"
#posi="iran2"
posi="./data/posi4gmt.dat"
waku="./data/waku.dat"
a=250
b=1250
c=250
d=1250
e=500
f=800

#####################################################################################################

hmi=0.4
hma=1.0

awk '{print  $1,$2, $3, $4}' < ${file}/${name} > tmpe.dat
#awk '{print  $1,$2, $3}' < model.dat > tmpe2.dat
makecpt -Cdrywet -Z  -T$hmi/$hma/1.e-2 > tmp.cpt
psxyz tmpe.dat -R$a/$b/$c/$d/$e/$f -E45/15 -JX7/4 -Jz0.050 -Ss0.5 -Ctmp.cpt -Ba200f200g200:"X[m]":/a200f200g200:"Y[m]":/a200f200g200:"Height[m]":wsNEZ+ -Xa1.i -Ya1.1i -P -K > tmp.ps
#psxyz tmpe2.dat -R$a/$b/$c/$d/$e/$f -E45/15 -JX7/5 -Jz0.25 -Sc0.05 -Gblack  -Xa1.i -Ya1.1i -P -K -O >> tmp.ps

#psxyz -R$a/$b/$c/$d/$e/$f -E45/20 -JX7/5 -Jz0.25 -Ss1 -Ctmp.cpt -Ba200f200g200:"X[m]":/a200f200g200:"Y[m]":/a25f25g25:"Height[m]":wsNEZ -Xa1.i -Ya1.7i -P -K  > tmp.ps

awk '{print $1, $2, $3}' <$posi | psxyz -E45/15 -JX7/4 -R -Jz0.050 -Sc0.2 -G0 -Xa1.i -Ya1.1i -P -K -O >> tmp.ps
awk '{print $1, $2, $3}' <$waku | psxyz -E45/15 -JX7/4 -R -W5 -Xa1.i -Ya1.1i -P -K -O >> tmp.ps


#grdmask tmpe.dat -R$a/$b/$c/$d -I40/10 -NNaN/1/1 -Gmask.grd
#grdmath tmpe.grd mask.grd OR = tmp.grd
#grdimage tmp.grd  -JX7/3 -R$a/$b/$c/$d -E45/45 -Ctmp.cpt -Xa1.i -Ya6.7i -P -K > tmp.ps
psscale -Ctmp.cpt -D1.9/2.4/4.0/0.2h -B5/:"[m/s]": -P  -O >> tmp.ps

gv tmp.ps
mv tpm.ps model.ps

#DFT,IDFT,IFFTスクリプト(H)::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# cp jh_fic.d inp2.dat #source J'
cp signal.d inp2.dat

i=1010
j=`expr $i - 1000`
max=30

while [ $j -le $max ]; do
cp hz$i.d inp1.dat #H'

ifort -parallel -o f_to_d parameter_wafe.f90 window_function.f90 f_to_d_h2.f90 -lfftw3 #-check all -CB -traceback -g -r8
./f_to_d

cp out1.dat HZ_w$i.d
cp out2.dat JZ_w$i.d
cp out3.dat GXh_w$i.d
cp out4.dat absHZ_w$i.d
cp out5.dat absJZ_w$i.d
cp out6.dat absGXh_w$i.d

cp conjg_hzw.dat conjg_hzw$i.d
cp conjg_jzw.dat conjg_jzw$i.d
cp conjg_gxhw.dat conjg_gxhw$i.d

cp invGH.dat HZ_t$i.d
cp invGJ.dat JZ_t$i.d
cp invGG.dat GXh_t$i.d
cp absHZ_t.dat absHZ_t$i.d
cp absJZ_t.dat absJZ_t$i.d
cp absGXh_t.dat absGXh_t$i.d
i=`expr $i + 10`
j=`expr $j + 10`
done

rm out*.dat invG*.dat conjg*.dat abs*.dat








#DFT,IDFT,IFFTスクリプト(E)::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# cp jh_fic.d inp2.dat #source J'
cp signal.d inp2.dat

i=1010
j=`expr $i - 1000`
max=30

while [ $j -le $max ]; do
cp ex$i.d inp1.dat #E'

ifort -parallel -o f_to_d parameter_wafe.f90 window_function.f90 f_to_d_e2.f90 -lfftw3 #-check all -CB -traceback -g -r8  #-I/usr/local/include -lfftw3 -check all -CB -traceback -g -r8
./f_to_d
# freq
cp out1.dat EX_w$i.d
cp out2.dat JZ_w$i.d
cp out3.dat GXe_w$i.d
cp out4.dat absEX_w$i.d
cp out5.dat absJZ_w$i.d
cp out6.dat absGXe_w$i.d
# # conjg_ver
cp conjg_exw.dat conjg_exw$i.d
cp conjg_jzw.dat conjg_jzw$i.d
cp conjg_gxew.dat conjg_gxew$i.d
# # time
cp invGE.dat EX_t$i.d
cp invGJ.dat JZ_t$i.d
cp invGG.dat GXe_t$i.d

cp absEX_t.dat absEX_t$i.d
cp absJZ_t.dat absJZ_t$i.d
cp absGXe_t.dat absGXe_t$i.d
i=`expr $i + 10`
j=`expr $j + 10`
done

rm out*.dat invG*.dat conjg*.dat abs*.dat

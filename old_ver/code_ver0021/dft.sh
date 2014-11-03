#DFTスクリプト(E)::::::::::::::::::::::::::::::
cp jh_fic.d inp2.dat #source

i=1010
j=`expr $i - 1000`
max=30

while [ $j -le $max ]; do
cp ex$i.d inp1.dat
ifort -o f_to_d parameter.f90 taper.f90 f_to_d.f90 #-I/usr/local/include -lfftw3 #-check all -CB -traceback -g -r8
./f_to_d
cp out1.dat EX_w$i.dat
cp out2.dat JX_w$i.dat
cp out3.dat GX_w$i.dat
i=`expr $i + 10`
j=`expr $j + 10`
done

rm out*.dat
make posi model wave rec
export OMP_NUM_THREADS=4
./fwi3d_posi.out
./fwi3d_model.out
./fwi3d_wave.out
./fwi3d_rec.out

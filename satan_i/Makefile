# Name of compiled program
PROG1=./fwi3d_main.out
PROG2=./fwi3d_wave.out
PROG3=./fwi3d_posi.out
PROG4=./fwi3d_rec.out
PROG5=./fwi3d_noise.out
PROG6=./fwi3d_convert.out
PROG8=./fwi3d_quadrup.out
PROG7=./fwi3d_inv.out
PROG9=./fwi3d_freq.out
PROG10=./fwi3d_model.out
# Name of compiler
LOCAL  =gcc -fopenmp -O3 -lfftw3 -lm
REMOTE =icc -openmp -parallel -O3 -lfftw3
SUPER=cc -h omp -h autothread -O3 -h msgs -h pic -dynamic
CC =$(LOCAL)
CC2=gfortran

# current code Version
VERSION = 1.4

help:
	@echo "Available makes:"
	@echo "posi => model => wave => rec => noise => inv => 3Cinv, clean, cleandata"

# Object files for the executable
cdriver		=   fwi3d_cpu_cpumain_v$(VERSION).o
function		= fwi3d_cpu_function_v$(VERSION).o
lattice		=   fwi3d_cpu_lattice_v$(VERSION).o
rec2rec		=   fwi3d_cpu_rec2rec_v$(VERSION).o
posiinf		=   fwi3d_cpu_posiinf_v$(VERSION).o
wavetra		=   fwi3d_cpu_waveform_v$(VERSION).o
convert		=   fwi3d_cpu_convert_v$(VERSION).o
quadrup		=   fwi3d_cpu_quadrup_v$(VERSION).o
inv_fwi		=   fwi3d_cpu_inv_fwi_v$(VERSION).o
inv_3C		=   fwi3d_cpu_3C_inv_v$(VERSION).o
cpmlabc		=   fwi3d_cpu_cpmlabc_v$(VERSION).o
frequency		= fwi3d_cpu_makewave_v$(VERSION).o
sigmodel		= fwi3d_cpu_sigmodel_v$(VERSION).o

3Cinv:		$(function) $(lattice) $(cpmlabc) $(inv_3C)
	$(CC) $(function) $(lattice) $(cpmlabc) $(inv_3C) -o $(PROG7)

inv:		$(function) $(lattice) $(cpmlabc) $(inv_fwi)
	$(CC) $(function) $(lattice) $(cpmlabc) $(inv_fwi) -o $(PROG7)

conv:		$(convert)
	$(CC) $(convert) -o $(PROG6)

model:		$(sigmodel)
	$(CC) $(sigmodel) -o $(PROG10)

freq:		$(frequency)
	$(CC) $(frequency) -o $(PROG9)

quad:		$(quadrup)
	$(CC2) $(quadrup) -o $(PROG8)

rec:		$(function) $(lattice) $(cpmlabc) $(rec2rec)
	$(CC) $(function) $(lattice) $(cpmlabc) $(rec2rec) -o $(PROG4)

wave:		$(wavetra)
	$(CC) $(wavetra) -o $(PROG2)

posi:		$(posiinf)
	$(CC) $(posiinf) -o $(PROG3)

clean:
	rm -f *.out *.o *.mod

cleandata:
	rm -f *.out *.o *.mod ./data/* ./data2/*

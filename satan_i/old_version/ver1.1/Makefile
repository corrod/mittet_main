# Name of compiled program
PROG1=./fwi3d_main.out
PROG2=./fwi3d_wave.out
PROG3=./fwi3d_posi.out
PROG4=./fwi3d_rec.out
PROG5=./fwi3d_noise.out
PROG6=./fwi3d_convert.out
PROG8=./fwi3d_quadrup.out
PROG7=./fwi3d_inv.out
# Name of compiler
LOCAL  =gcc -O3
REMOTE =icc -openmp -O3
SUPER=cc -h omp -h autothread -O3 -h msgs -h pic -dynamic
CC =$(LOCAL)
CC2=gfortran

# current code Version
VERSION = 1.1

help:
	@echo "Available makes:"
	@echo "posi => wave => rec => noise => main, clean, cleandata"

# Object files for the executable
cdriver		=   fwi3d_cpu_cpumain_v$(VERSION).o
function		= fwi3d_cpu_function_v$(VERSION).o
lattice		=   fwi3d_cpu_lattice_v$(VERSION).o
eh_mur1		=   fwi3d_cpu_eh_mur1_v$(VERSION).o
cerjan1		=   fwi3d_cpu_cerjan1_v$(VERSION).o
rec2rec		=   fwi3d_cpu_rec2rec_v$(VERSION).o
posiinf		=   fwi3d_cpu_posiinf_v$(VERSION).o
wavetra		=   fwi3d_cpu_wavetra_v$(VERSION).o
convert		=   fwi3d_cpu_convert_v$(VERSION).o
quadrup		=   fwi3d_cpu_quadrup_v$(VERSION).o
inv_fwi		=   fwi3d_cpu_inv_fwi_v$(VERSION).o
cpmlabc		=   fwi3d_cpu_cpmlabc_v$(VERSION).o

inv:		$(function) $(cerjan1) $(lattice) $(inv_fwi)
	$(CC) $(function) $(cerjan1) $(lattice) $(inv_fwi) -o $(PROG7)

conv:		$(convert)
	$(CC) $(convert) -o $(PROG6)

quad:		$(quadrup)
	$(CC2) $(quadrup) -o $(PROG8)

rec:		$(function) $(cerjan1) $(eh_mur1) $(lattice) $(cpmlabc) $(rec2rec) 
	$(CC) $(function) $(cerjan1) $(eh_mur1) $(lattice) $(cpmlabc) $(rec2rec) -o $(PROG4)

wave:		$(wavetra)
	$(CC) $(wavetra) -o $(PROG2)

posi:		$(posiinf)
	$(CC) $(posiinf) -o $(PROG3)

clean:	
	rm -f $(PROG1) $(PROG2) $(PROG3) $(PROG4) $(PROG5) $(PROG6) $(PROG7) *.o *.mod

cleandata:	
	rm -f $(PROG1) $(PROG2) $(PROG3) $(PROG4) $(PROG5) $(PROG6) $(PROG7) *.o *.mod ./data/* ./data2/*

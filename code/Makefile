# makefile for fwd simulation
TARGET = mittet
OBJECTS = parameter.o set_condition.o waveform.o model.o output.o e_field4.o h_field4.o cpml.o cpml2.o mur2.o cerjan.o pml.o main2.o
MOD_FILES  = parameter.mod
#FC = gfortran
FC = ifort

FFLAGS =
LDFLAGS =

# for gfortran
ifeq (${FC},gfortran)
	FFLAGS  += -fimplicit-none -fbounds-check
	LDFLAGS += -fopenmp -llapack -lblas
endif

# for ifort
 ifeq (${FC},ifort)
# 	MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
# 	FFLAGS  += -I${MKLROOT}/include/ia32 -I${MKLROOT}/include #-check all -CB -traceback -g -r8
	LDFLAGS += -lfftw3
# 	LDFLAGS += -L${MKLROOT}/lib/ia32 ${MKLROOT}/lib/ia32/libmkl_blas95.a
# 	LDFLAGS += ${MKLROOT}/lib/ia32/libmkl_lapack95.a
# 	LDFLAGS += -lmkl_intel -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm
 endif


.SUFFIXES : .o .f90
.f90.o:
	${FC} -c $<

${TARGET} : ${OBJECTS}
	${FC} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS}


.PHONY: clean
clean:
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES}
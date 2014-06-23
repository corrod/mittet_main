TARGET = mittet
OBJECTS = parameter.o main2.o set_condition.o waveform.o model.o \
	e_field4.o e_cpml4.o h_field4.o h_cpml4.o output.o
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
	MKLROOT = /opt/intel/composer_xe_2011_sp1.7.256/mkl
	FFLAGS  += -I${MKLROOT}/include/ia32 -I${MKLROOT}/include
	LDFLAGS += -L${MKLROOT}/lib/ia32 ${MKLROOT}/lib/ia32/libmkl_blas95.a
	LDFLAGS += ${MKLROOT}/lib/ia32/libmkl_lapack95.a
	LDFLAGS += -lmkl_intel -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm
endif


.SUFFIXES : .o .f90
.f90.o:
	${FC} -c $<

${TARGET} : ${OBJECTS}
	${FC} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS}


.PHONY: clean
clean:
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES}
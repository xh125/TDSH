#compiler
FCC	    =	ifort -fpp # -check all -pg -traceback
FFLAGS      = -g -O2

#scr directory
VPATH       = ../src:../utility/nomashifwannier

#MKL libraries
MKLLIB      = -L${MKLROOT}/lib/intel64
MKLINCLUDE  = -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
FCCFLAG	    = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -liomp5\
              -lmkl_blas95_lp64 -lmkl_lapack95_lp64  

PROG	    = SHexciton.x
#SRC
SRC      = constants.f90 io.f90 random.f90 utility.f90 parameters.f90 hamiltonian.f90 \
           dynamic.f90 main.f90
#objects
OBJ       = $(SRC:.f90=.o)

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(FCC) $(FFLAGS) ${MKLLIB} ${MKLINCLUDE} ${FCCFLAG} -c $<

##
${PROG} : ${OBJ}	
	${FCC} -o ${PROG} ${OBJ} ${MKLINCLUDE} ${MKLLIB} ${FCCFLAG}
	cp -f ${PROG} ../bin
  
#Compiler utility
Nshift.x : setposcar.o
	${FCC} -o setposcar.o
setposcar.o : setposcar.f90
	${FCC} -c setposcar.f90
	cp Nshift.x ../bin

all:${PROG} Nshift.x
	  
#make clean
.PHONY : clean all
clean:
	-rm ${PROG} ${OBJ} *.mod ../bin/*
	

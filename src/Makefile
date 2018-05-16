#compiler
FCC = ifort # -check all -pg -traceback

#scr directory
#VPATH       = src:utility/nomashifwannier

#MKL libraries
MKLLIB	= -L${MKLROOT}/lib/intel64
MKLINCLUDE = -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
FCCFLAG	= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -liomp5 \
          -lmkl_blas95_lp64 -lmkl_lapack95_lp64  

#SRC
SRC = constants.f90 io.f90 parameters.f90 hamiltonian.f90 dynamic.f90 \
              random.f90 utility.f90 main.f90
#objects
OBJ = $(SRC:.f90=.o)
PROG = SHexciton.x

##
${PROG}: ${OBJ}
	${FCC} -o ${OBJ}
#	cp -f ${PROG}.x ../bin
main.o : main.f90 random.o constants.o io.o parameters.o hamiltonian.o dynamic.o random.o utility.o
	${FCC} -c main.f90
random.o : random.f90
	${FCC} -c random.o
constants.o :constants.f90
	${FCC} -c constants.f90
sh.o : constants.o
	${FCC} -c sh.f90 constants.o
utility.o : sh.o constants.o
	${FCC} -c utility.f90 sh.o constants.o

	  
#Compiler utility
Nshift : setposcar.o
	${FCC} -o setposcar.o
setposcar.o : setposcar.f90
	${FCC} -c setposcar.f90
  
all:${PROG} Nshift
  
#make clean
.PHONY : clean all
clean:
	-rm ${PROG} ${OBJ} *.mod ../bin/*
	

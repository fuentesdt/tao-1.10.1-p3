# $Id$ 

ALL: rosenbrock1

CFLAGS		 = 
FFLAGS		 = 
CPPFLAGS         =
FPPFLAGS         =
LOCDIR		 = examples/

include ${TAO_DIR}/bmake/tao_common

rosenbrock1: rosenbrock1.o  tao_chkopts
	-${CLINKER} -o rosenbrock1 rosenbrock1.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} rosenbrock1.o

rosenbrock1f: rosenbrock1f.o  tao_chkopts
	-${FLINKER} -o rosenbrock1f rosenbrock1f.o ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_SNES_LIB}
	${RM} rosenbrock1f.o

minsurf1: minsurf1.o  tao_chkopts
	-${CLINKER} -o minsurf1 minsurf1.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} minsurf1.o

plate2: plate2.o  tao_chkopts
	-${CLINKER} -o plate2 plate2.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} plate2.o



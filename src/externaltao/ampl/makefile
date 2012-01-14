all: driver
#AMPL_SOLVERS_PATH = /home/aotools/blmvmneos/blmvm_neos_ampl/src/linux/amplsolvers
AMPL_SOLVERS_PATH=/home/neosotc/software/solvers
AMPL_LIB = ${AMPL_SOLVERS_PATH}/amplsolver.a
#MPI_INCLUDE = ${PETSC_DIR}/src/sys/src/mpiuni 

CFLAGS		 =
FFLAGS		 = 
CPPFLAGS         = -I${AMPL_SOLVERS_PATH} -I${TAO_DIR}/include -I. -I${PETSC_DIR}/include -I${PETSC_DIR}/bmake/linux ${MPI_INCLUDE}
FPPFLAGS         =
LOCDIR		 = 

SOURCEC = blmvm_ampl.c
SOURCEH = 
OBJSC = blmvm_ampl.o
include ${TAO_DIR}/bmake/tao_common

# See $TAO_DIR/bmake/tao_common_basic for definitions of
# example categories.
# include ${TAO_DIR}/bmake/tao_common

driver: blmvm_ampl.o tao_chkopts
	-${CLINKER} ${CPPFLAGS} -o blmvm blmvm_ampl.o ${TAO_LIB} ${PETSC_SNES_LIB} ${AMPL_LIB} -ldl
	rm -f *.o

staticdriver: blmvm_ampl.o tao_chkopts
	-${CLINKER} ${CPPFLAGS} -static -o blmvm_static blmvm_ampl.o ${TAO_LIB} ${PETSC_SNES_LIB} ${AMPL_LIB} -ldl
	rm -f *.o

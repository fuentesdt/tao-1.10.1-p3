all: lennard-jones_ga

MANSEC    = gaapp
FFLAGS    =
CFLAGS    = ${GA_INCLUDE} 
CPPFLAGS  =  
SOURCEC   = 
SOURCEF   = 
EXAMPLESC = lennard-jones_ga.c lennard-jones_ga2.c rosenbrock_ga.c
LOCDIR    = src/externaltao/globalarraytao/examples/
TUTORIALS = lennard-jones_ga2.c lennard-jones_ga.c rosenbrock_ga.c 

include ${TAO_DIR}/bmake/tao_common

#Uncomment to avoid need for PETSc headers
#PETSC_CCPPFLAGS=-DTAO_NO_PETSC ${MPI_INCLUDE}


lennard-jones_ga: lennard-jones_ga.o tao_chkopts
	-${FLINKER} -o lennard-jones_ga lennard-jones_ga.o ${GA_LIB} 
	${RM} lennard-jones_ga.o 

lennard-jones_ga2: lennard-jones_ga2.o tao_chkopts
	-${FLINKER} -o lennard-jones_ga2 lennard-jones_ga2.o ${GA_LIB} 
	${RM} lennard-jones_ga2.o

rosenbrock_ga: rosenbrock_ga.o tao_chkopts
	-${CLINKER} -o rosenbrock_ga rosenbrock_ga.o ${GA_LIB} 
	${RM} rosenbrock_ga.o



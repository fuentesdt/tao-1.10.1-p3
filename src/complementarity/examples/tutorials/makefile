
ALL: minsurf1 blackscholes

CFLAGS	        =
FFLAGS	        =
CPPFLAGS        =
FPPFLAGS        =	 
LOCDIR          = src/complementarity/examples/tutorials/
TUTORIALS	= minsurf1.c blackscholes.c
MANSEC		= Complementarity

TESTEXAMPLES_C  = minsurf1.tao runminsurf1 minsurf1.rm blackscholes.tao runblackscholes blackscholes.rm

EXAMPLESC       = minsurf1.c blackscholes.c
EXAMPLESF       = 
DIRS	        = 

TMPROSC          = ../../../taouni/simple.c
TMPROSH          = ../../../taouni/simple.h

include ${TAO_DIR}/bmake/tao_common

minsurf1: minsurf1.o  tao_chkopts
	-${CLINKER} -o minsurf1 minsurf1.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} minsurf1.o

blackscholes: blackscholes.o  tao_chkopts
	-${CLINKER} -o blackscholes blackscholes.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} blackscholes.o

runminsurf1:
	-@${MPIEXEC} -np 1 ./minsurf1 -tao_smonitor \
          -tao_method tao_ssils > minsurf1_1.tmp 2>&1;   \
          if (diff output/minsurf1_1.out minsurf1_1.tmp) then true; \
          else echo "Possible problem with minsurf1_1, diffs above"; fi; \
          ${RM} -f minsurf1_1.tmp

runblackscholes:
	-@${MPIEXEC} -np 1 ./blackscholes -mt 3 -tao_smonitor \
          -tao_method tao_ssils > blackscholes_1.tmp 2>&1;   \
          if (diff output/blackscholes_1.out blackscholes_1.tmp) then true; \
          else echo "Possible problem with blackscholes_1, diffs above"; fi; \
          ${RM} -f blackscholes_1.tmp


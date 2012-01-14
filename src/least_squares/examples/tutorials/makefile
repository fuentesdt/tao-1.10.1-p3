

ALL: enzreac1

MANSEC		 = LeastSquares
CFLAGS		 = 
FFLAGS		 = 
CPPFLAGS         =
FPPFLAGS         =
LOCDIR		 = src/least_squares/examples/tutorials/
#TUTORIALS        = coating.c enzreac1.c chebyq1.c bratu2.c chwirut.c

# See $TAO_DIR/bmake/tao_common_basic for definitions of
# example categories.
TESTEXAMPLES_1	 = enzreac1.tao runenzreac1 \
	           enzreac1.rm\
		   chebyq1.tao \
		   chebyq1.rm \
		   bratu2.tao runbratu2 \
	           bratu2.rm \
		   coating.tao runcoating coating.rm \
		   clean

TESTEXAMPLES_4	 = enzreac1.tao runenzreac1 \
	           enzreac1.rm\
		   chebyq1.tao \
		   chebyq1.rm \
		   coating.tao runcoating coating.rm \
		   clean
BUILDEXAMPLES_1	 = enzreac1.tao enzreac1.rm\
		   chebyq1.tao chebyq1.rm \
		   bratu2.tao bratu2.rm \
		   coating.tao coating.rm  clean

BUILDEXAMPLES_4	 = enzreac1.tao enzreac1.rm\
		   chebyq1.tao chebyq1.rm \
		   coating.tao coating.rm clean

EXAMPLESC	 = enzreac1.c chebyq1.c coating.c bratu2.c chwirut.c
EXAMPLESF	 = 
EXAMPLESCH	 = 
EXAMPLESFH       = 
DOCS		 = 
DIRS		 = 

include ${TAO_DIR}/bmake/tao_common

enzreac1: ./enzreac1.o  tao_chkopts
	-${CLINKER} -o ./enzreac ./enzreac1.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} ./enzreac1.o

chebyq: chebyq.o  tao_chkopts
	-${CLINKER} -o chebyq chebyq.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} chebyq.o

chebyq1: chebyq1.o  tao_chkopts
	-${CLINKER} -o chebyq1 chebyq1.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} chebyq1.o

coating: coating.o  tao_chkopts
	-${CLINKER} -o coating coating.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} coating.o

bratu2: ./bratu2.o tao_chkopts
	-${CLINKER} -o ./bratu2 ./bratu2.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} ./bratu2.o


chwirut: chwirut.o tao_chkopts
	-${CLINKER} -o chwirut chwirut.o ${TAO_LIB} ${PETSC_SNES_LIB}


runenzreac1:
	-@${MPIEXEC} -np 1 enzreac -tao_smonitor \
           > enzreac1_1.tmp 2>&1;   \
	   if (diff output/enzreac1_1.out enzreac1_1.tmp) then true; \
	   else echo "Possible problem with enzreac1_1, diffs above"; fi; \
	   ${RM} -f enzreac1_1.tmp

runchebyq1:
	-@${MPIEXEC} -np 1 chebyq1 \
           > chebyq1_1.tmp 2>&1;   \
	   if (diff output/chebyq1_1.out chebyq1_1.tmp) then true; \
	   else echo "Possible problem with chebyq1_1, diffs above"; fi; \
	   ${RM} -f chebyq1_1.tmp

runcoating:
	-@${MPIEXEC} -np 1 coating -tao_smonitor  \
           > coating_1.tmp 2>&1;   \
	   if (diff output/coating_1.out coating_1.tmp) then true; \
	   else echo "Possible problem with coating_1, diffs above"; fi; \
	   ${RM} -f coating_1.tmp

runbratu2:
	-@${MPIEXEC} -np 1 bratu2 -tao_smonitor \
           > bratu2_1.tmp 2>&1;   \
	   if (diff output/bratu2_1.out bratu2_1.tmp) then true; \
	   else echo "Possible problem with bratu2_1, diffs above"; fi; \
	   ${RM} -f bratu2_1.tmp







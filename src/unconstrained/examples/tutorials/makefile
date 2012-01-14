# $Id$ 

ALL: eptorsion1

MANSEC		 = Unconstrained
CFLAGS		 = 
FFLAGS		 = 
CPPFLAGS         =
FPPFLAGS         =
LOCDIR		 = src/unconstrained/examples/tutorials/
TUTORIALS	 = rosenbrock1.c eptorsion1.c minsurf1.c eptorsion2.c\
	           minsurf2.c eptorsion2f.F

# See $PETSC_DIR/bmake/test for definitions of example categories.

TESTEXAMPLES_C	 = eptorsion1.tao runeptorsion1 runeptorsion1_2 \
	           runeptorsion1_3 runeptorsion1_4 eptorsion1.rm\
		   eptorsion2.tao \
		   runeptorsion2_5 runeptorsion2_6 eptorsion2.rm \
		   minsurf2.tao runminsurf2 runminsurf2_2 runminsurf2_3 \
	           runminsurf2_4 minsurf2.rm \
		   rosenbrock1.tao runrosenbrock1 rosenbrock1.rm \
		   clean
TESTEXAMPLES_FORTRAN	 = eptorsion2f.tao runeptorsion2f_3 runeptorsion2f_3\
	           runeptorsion2f_3 eptorsion2f.rm \
                   rosenbrock1f.tao runrosenbrock1f rosenbrock1f.rm clean
TESTEXAMPLES_C_X11_MPIUNI	 = eptorsion1.tao runeptorsion1 runeptorsion1_2 \
		   runeptorsion1_3 runeptorsion1_4 eptorsion1.rm \
		   rosenbrock1.tao runrosenbrock1 rosenbrock1.rm \
		   minsurf1.tao runminsurf1 minsurf1.rm \
		   clean
TESTEXAMPLES_FORTRAN_MPIUNI	 = rosenbrock1f.tao runrosenbrock1f rosenbrock1f.rm clean
EXAMPLESC	 = eptorsion1.c eptorsion2.c minsurf2.c rosenbrock1.c
EXAMPLESF	 = eptorsion2f.F rosenbrock1f.F
EXAMPLESCH	 = 
EXAMPLESFH       = eptorsion2f.h rosenbrock1f.h
DOCS		 = 
DIRS		 = 

include ${TAO_DIR}/bmake/tao_common

eptorsion1: eptorsion1.o  tao_chkopts
	-${CLINKER} -o eptorsion1 eptorsion1.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} eptorsion1.o

eptorsion2: eptorsion2.o  tao_chkopts
	-${CLINKER} -o eptorsion2 eptorsion2.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} eptorsion2.o

eptorsion2f: eptorsion2f.o  tao_chkopts
	-${FLINKER} -o eptorsion2f eptorsion2f.o ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_SNES_LIB}
	${RM} eptorsion2f.o

minsurf1: minsurf1.o  tao_chkopts
	-${CLINKER} -o minsurf1 minsurf1.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} minsurf1.o

minsurf2: minsurf2.o  tao_chkopts
	-${CLINKER} -o minsurf2 minsurf2.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} minsurf2.o

rosenbrock1a: rosenbrock1a.o  tao_chkopts
	-${CLINKER} -o rosenbrock1a rosenbrock1a.o ${TAO_BLIB} ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} rosenbrock1a.o

rosenbrock1: rosenbrock1.o  tao_chkopts
	-${CLINKER} -o rosenbrock1 rosenbrock1.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} rosenbrock1.o

rosenbrock1f: rosenbrock1f.o  tao_chkopts
	-${FLINKER} -o rosenbrock1f rosenbrock1f.o ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_SNES_LIB}
	${RM} rosenbrock1f.o

runeptorsion1:
	-@${MPIEXEC} -np 1 ./eptorsion1 -tao_smonitor \
           -tao_method tao_lmvm > eptorsion1_1.tmp 2>&1;   \
	   if (diff output/eptorsion1_1.out eptorsion1_1.tmp) then true; \
	   else echo "Possible problem with eptorsion1_1, diffs above"; fi; \
	   ${RM} -f eptorsion1_1.tmp

runeptorsion1_2:
	-@${MPIEXEC} -np 1 ./eptorsion1 -tao_smonitor \
           -tao_method tao_nls > eptorsion1_2.tmp 2>&1;   \
	   if (diff output/eptorsion1_2.out eptorsion1_2.tmp) then true; \
	   else echo "Possible problem with eptorsion1_2, diffs above"; fi; \
	   ${RM} -f eptorsion1_2.tmp

runeptorsion1_3:
	-@${MPIEXEC} -np 1 ./eptorsion1 -tao_smonitor \
           -tao_method tao_cg -tao_cg_type prp > eptorsion1_3.tmp 2>&1;   \
	   if (diff output/eptorsion1_3.out eptorsion1_3.tmp) then true; \
	   else echo "Possible problem with eptorsion1_3, diffs above"; fi; \
	   ${RM} -f eptorsion1_3.tmp

runeptorsion1_4:
	-@${MPIEXEC} -np 1 ./eptorsion1 -tao_smonitor \
           -tao_method tao_ntr > eptorsion1_4.tmp 2>&1;   \
	   if (diff output/eptorsion1_4.out eptorsion1_4.tmp) then true; \
	   else echo "Possible problem with eptorsion1_4, diffs above"; fi; \
	   ${RM} -f eptorsion1_4.tmp

runeptorsion1f:
	-@${MPIEXEC} -np 1 ./eptorsion1f -tao_smonitor > eptorsion1f_1.tmp 2>&1;\
	   if (diff output/eptorsion1f_1.out eptorsion1f_1.tmp) then true; \
	   else echo "Possible problem with eptorsion1f_1, diffs above"; fi; \
	   ${RM} -f eptorsion1f_1.tmp

runeptorsion2:
	-@${MPIEXEC} -np 1 ./eptorsion2 -tao_smonitor \
	   -tao_method tao_nls > eptorsion2_1.tmp 2>&1;	\
	   if (diff output/eptorsion2_1.out eptorsion2_1.tmp) then true; \
	   else echo "Possible problem with eptorsion2_1, diffs above"; fi; \
	   ${RM} -f eptorsion2_1.tmp

runeptorsion2_2:
	-@${MPIEXEC} -np 2 ./eptorsion2 -tao_smonitor \
	   -tao_method tao_nls > eptorsion2_2.tmp 2>&1;	\
	   if (diff output/eptorsion2_2.out eptorsion2_2.tmp) then true; \
	   else echo "Possible problem with eptorsion2_2, diffs above"; fi; \
	   ${RM} -f eptorsion2_2.tmp

runeptorsion2_3:
	-@${MPIEXEC} -np 1 ./eptorsion2 -tao_smonitor \
	   -tao_method tao_ntr > eptorsion2_3.tmp 2>&1;	\
	   if (diff output/eptorsion2_3.out eptorsion2_3.tmp) then true; \
	   else echo "Possible problem with eptorsion2_3, diffs above"; fi; \
	   ${RM} -f eptorsion2_3.tmp

runeptorsion2_4:
	-@${MPIEXEC} -np 3 ./eptorsion2 -mx 16 -my 16 -tao_smonitor \
	   -tao_method tao_ntr > eptorsion2_4.tmp 2>&1;	\
	   if (diff output/eptorsion2_4.out eptorsion2_4.tmp) then true; \
	   else echo "Possible problem with eptorsion2_4, diffs above"; fi; \
	   ${RM} -f eptorsion2_4.tmp

runeptorsion2_5:
	-@${MPIEXEC} -np 1 ./eptorsion2 -tao_smonitor \
	   -tao_method tao_lmvm > eptorsion2_5.tmp 2>&1;	\
	   if (diff output/eptorsion2_5.out eptorsion2_5.tmp) then true; \
	   else echo "Possible problem with eptorsion2_5, diffs above"; fi; \
	   ${RM} -f eptorsion2_5.tmp

runeptorsion2_6:
	-@${MPIEXEC} -np 3 ./eptorsion2 -mx 16 -my 16 -tao_smonitor \
	   -tao_method tao_lmvm > eptorsion2_6.tmp 2>&1;	\
	   if (diff output/eptorsion2_6.out eptorsion2_6.tmp) then true; \
	   else echo "Possible problem with eptorsion2_6, diffs above"; fi; \
	   ${RM} -f eptorsion2_6.tmp

runeptorsion2f:
	-@${MPIEXEC} -np 1 ./eptorsion2f -tao_smonitor \
	   -tao_method tao_nls > eptorsion2f_1.tmp 2>&1;	\
	   if (diff output/eptorsion2f_1.out eptorsion2f_1.tmp) then true; \
	   else echo "Possible problem with eptorsion2f_1, diffs above"; fi; \
	   ${RM} -f eptorsion2f_1.tmp

runeptorsion2f_2:
	-@${MPIEXEC} -np 2 ./eptorsion2f -tao_smonitor \
	   -tao_method tao_nls > eptorsion2f_2.tmp 2>&1;	\
	   if (diff output/eptorsion2f_2.out eptorsion2f_2.tmp) then true; \
	   else echo "Possible problem with eptorsion2f_2, diffs above"; fi; \
	   ${RM} -f eptorsion2f_2.tmp

runeptorsion2f_3:
	-@${MPIEXEC} -np 3 ./eptorsion2f -mx 16 -my 16 -tao_smonitor \
	   -tao_method tao_lmvm > eptorsion2f_3.tmp 2>&1;	\
	   if (diff output/eptorsion2f_3.out eptorsion2f_3.tmp) then true; \
	   else echo "Possible problem with eptorsion2f_3, diffs above"; fi; \
	   ${RM} -f eptorsion2f_3.tmp

runminsurf1:
	-@${MPIEXEC} -np 1 ./minsurf1 -tao_smonitor -mx 10 -my 8\
           -tao_method tao_nls > minsurf1_1.tmp 2>&1;   \
	   if (diff output/minsurf1_1.out minsurf1_1.tmp) then true; \
	   else echo "Possible problem with minsurf1_1, diffs above"; fi; \
	   ${RM} -f minsurf1_1.tmp

runminsurf2:
	-@${MPIEXEC} -np 1 ./minsurf2 -tao_smonitor -mx 10 -my 8\
           -tao_method tao_lmvm > minsurf2_1.tmp 2>&1;   \
	   if (diff output/minsurf2_1.out minsurf2_1.tmp) then true; \
	   else echo "Possible problem with minsurf2_1, diffs above"; fi; \
	   ${RM} -f minsurf2_1.tmp

runminsurf2_2:
	-@${MPIEXEC} -np 2 ./minsurf2 -tao_smonitor \
           -tao_method tao_nls > minsurf2_2.tmp 2>&1;   \
	   if (diff output/minsurf2_2.out minsurf2_2.tmp) then true; \
	   else echo "Possible problem with minsurf2_2, diffs above"; fi; \
	   ${RM} -f minsurf2_2.tmp

runminsurf2_3:
	-@${MPIEXEC} -np 3 ./minsurf2 -tao_smonitor -mx 10 -my 10\
           -tao_method tao_cg -tao_cg_type fr > minsurf2_3.tmp 2>&1;   \
	   if (diff output/minsurf2_3.out minsurf2_3.tmp) then true; \
	   else echo "Possible problem with minsurf2_3, diffs above"; fi; \
	   ${RM} -f minsurf2_3.tmp

runminsurf2_4:
	-@${MPIEXEC} -np 2 ./minsurf2 -tao_smonitor -mx 6 -my 8\
           -tao_method tao_ntr > minsurf2_4.tmp 2>&1;   \
	   if (diff output/minsurf2_4.out minsurf2_4.tmp) then true; \
	   else echo "Possible problem with minsurf2_4, diffs above"; fi; \
	   ${RM} -f minsurf2_4.tmp

runminsurf2_5:
	-@${MPIEXEC} -np 3 ./minsurf2 -tao_smonitor -mx 23 -my 17\
           -tao_method tao_nls -tao_fddefault > minsurf2_5.tmp 2>&1;   \
	   if (diff output/minsurf2_5.out minsurf2_5.tmp) then true; \
	   else echo "Possible problem with minsurf2_5, diffs above"; fi; \
	   ${RM} -f minsurf2_5.tmp

runminsurf2_6:
	-@${MPIEXEC} -np 1 ./minsurf2 -tao_smonitor -mx 4 -my 20 -random 2 \
           -tao_method tao_ntr > minsurf2_6.tmp 2>&1;   \
	   if (diff output/minsurf2_6.out minsurf2_6.tmp) then true; \
	   else echo "Possible problem with minsurf2_6, diffs above"; fi; \
	   ${RM} -f minsurf2_6.tmp

runrosenbrock1:
	-@${MPIEXEC} -np 1 ./rosenbrock1 -tao_smonitor \
	   -tao_method tao_lmvm > rosenbrock1_1.tmp 2>&1;	\
	   if (diff output/rosenbrock1_1.out rosenbrock1_1.tmp) then true; \
	   else echo "Possible problem with rosenbrock1_1, diffs above"; fi; \
	   ${RM} -f rosenbrock1_1.tmp

runrosenbrock1f:
	-@${MPIEXEC} -np 1 ./rosenbrock1f -tao_smonitor \
	   -tao_method tao_lmvm > rosenbrock1f_1.tmp 2>&1;	\
	   if (diff output/rosenbrock1f_1.out rosenbrock1f_1.tmp) then true; \
	   else echo "Possible problem with rosenbrock1f_1, diffs above"; fi; \
	   ${RM} -f rosenbrock1f_1.tmp



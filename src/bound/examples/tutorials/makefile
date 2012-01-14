
ALL: plate2

MANSEC		 = Bound
CFLAGS		 = 
FFLAGS		 = 
CPPFLAGS         =
FPPFLAGS         =
LOCDIR		 = src/bound/examples/tutorials/
TUTORIALS	 = plate2.c jbearing2.c plate2f.F

# See $TAO_DIR/bmake/tao_common_variables for definitions of
# example categories.

TESTEXAMPLES_C	 = plate2.tao runplate2 runplate2_2 runplate2_3 plate2.rm\
		   jbearing2.tao runjbearing2 runjbearing2_2 \
		   runjbearing2_3 jbearing2.rm clean
TESTEXAMPLES_FORTRAN  = plate2f.tao runplate2f runplate2f_2 plate2f.rm clean
EXAMPLESC	 = plate2.c jbearing2.c
EXAMPLESF	 = plate2f.F
EXAMPLESCH	 = 
EXAMPLESFH       = plate2f.h
DOCS		 = 
DIRS		 = 

include ${TAO_DIR}/bmake/tao_common

plate2: plate2.o  tao_chkopts
	-${CLINKER} -o plate2 plate2.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} plate2.o

plate2f: plate2f.o  tao_chkopts
	-${FLINKER} -o plate2f plate2f.o ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_SNES_LIB}
	${RM} plate2f.o

jbearing2: jbearing2.o  tao_chkopts
	-${CLINKER} -o jbearing2 jbearing2.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} jbearing2.o

runplate2:
	-@${MPIEXEC} -np 1 ./plate2 -tao_smonitor -mx 8 -my 6 -bmx 3 -bmy 3 \
           -bheight 0.2 -tao_method tao_tron > plate2_1.tmp 2>&1; \
	   if (diff output/plate2_1.out plate2_1.tmp) then true; \
	   else echo "Possible problem with plate2_1, diffs above"; fi; \
	   ${RM} -f plate2_1.tmp

runplate2_2:
	-@${MPIEXEC} -np 2 ./plate2 -tao_smonitor -mx 8 -my 8 -bmx 2 -bmy 5 \
           -bheight 0.3 -tao_method tao_blmvm > plate2_2.tmp 2>&1; \
	   if (diff output/plate2_2.out plate2_2.tmp) then true; \
	   else echo "Possible problem with plate2_2, diffs above"; fi; \
	   ${RM} -f plate2_2.tmp

runplate2_3:
	-@${MPIEXEC} -np 3 ./plate2 -tao_smonitor -mx 8 -my 12 -bmx 4 -bmy 10 \
           -bheight 0.1 -tao_method tao_tron > plate2_3.tmp 2>&1; \
	   if (diff output/plate2_3.out plate2_3.tmp) then true; \
	   else echo "Possible problem with plate2_3, diffs above"; fi; \
	   ${RM} -f plate2_3.tmp

runplate2f:
	-@${MPIEXEC} -np 1 ./plate2f -tao_smonitor -mx 8 -my 6 -bmx 3 -bmy 3 \
           -bheight 0.2 -tao_method tao_blmvm > plate2f_1.tmp 2>&1; \
	   if (diff output/plate2f_1.out plate2f_1.tmp) then true; \
	   else echo "Possible problem with plate2f_1, diffs above"; fi; \
	   ${RM} -f plate2f_1.tmp

runplate2f_2:
	-@${MPIEXEC} -np 2 ./plate2f -tao_smonitor -mx 8 -my 6 -bmx 3 -bmy 3 \
           -bheight 0.2 -tao_method tao_blmvm > plate2f_2.tmp 2>&1; \
	   if (diff output/plate2f_2.out plate2f_2.tmp) then true; \
	   else echo "Possible problem with plate2f_2, diffs above"; fi; \
	   ${RM} -f plate2f_2.tmp

runjbearing2:
	-@${MPIEXEC} -np 1 ./jbearing2 -tao_smonitor -mx 8 -my 12\
	   -tao_method tao_tron > jbearing2_1.tmp 2>&1;	\
	   if (diff output/jbearing2_1.out jbearing2_1.tmp) then true; \
	   else echo "Possible problem with jbearing2_1, diffs above"; fi; \
	   ${RM} -f jbearing2_1.tmp

runjbearing2_2:
	-@${MPIEXEC} -np 2 ./jbearing2 -tao_smonitor -mx 50 -my 50 -ecc 0.99\
	   -tao_method tao_gpcg > jbearing2_2.tmp 2>&1;	\
	   if (diff output/jbearing2_2.out jbearing2_2.tmp) then true; \
	   else echo "Possible problem with jbearing2_2, diffs above"; fi; \
	   ${RM} -f jbearing2_2.tmp

runjbearing2_3:
	-@${MPIEXEC} -np 2 ./jbearing2 -tao_smonitor -mx 10 -my 16 -ecc 0.9 \
	   -tao_method tao_bqpip > jbearing2_3.tmp 2>&1;	\
	   if (diff output/jbearing2_3.out jbearing2_3.tmp) then true; \
	   else echo "Possible problem with jbearing2_3, diffs above"; fi; \
	   ${RM} -f jbearing2_3.tmp





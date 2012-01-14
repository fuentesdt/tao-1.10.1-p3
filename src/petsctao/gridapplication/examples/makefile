# $Id: makefile,v 1.1 2002/07/10 17:36:10 benson Exp benson $ 

ALL: minsurf3 jbearing3

#ADIC	    = /home/derivs/adic
#ADIC_ARCH  = linux

#CFLAGS   	= ${ADIC_INCLUDE} -Dad_GRAD_MAX=4  # Needed only if using the ADIC code
CFLAGS   	= ${ADIC_INCLUDE} -Dad_GRAD_MAX=16  # Needed only if using the ADIC code
CPPFLAGS 	= 
MANSEC		 = daapp
FFLAGS		 = 
FPPFLAGS         =
LOCDIR		 = src/petsctao/gridapplication/examples/
TUTORIALS	 = eptorsion3.c jbearing3.c minsurf3.c combustion3.c

TESTEXAMPLES_C	 = eptorsion3.tao \
	           eptorsion3.rm\
		   jbearing3.tao runjbearing3 runjbearing3_2 \
		   jbearing3.rm \
		   minsurf3.tao runminsurf3 runminsurf3_2 runminsurf3_3 \
	           minsurf3.rm \
		   combustion3.tao runcombustion3 runcombustion3_2 combustion3.rm \
		   clean
TESTEXAMPLES_FORTRAN = jbearing3f.tao runjbearing3f runjbearing3f_2 \
                   jbearing3f.rm clean

EXAMPLESC	 = eptorsion3.c jbearing3.c minsurf3.c combustion3.c

# See $PETSC_DIR/bmake/test for definitions of example categories.

include ${TAO_DIR}/bmake/tao_common

.PRECIOUS: %.ad.c

%.ad.c: %.c
	TaoComputeADICgradient $<

minsurf3: minsurf3.o  tao_chkopts
	-${CLINKER} -o minsurf3 minsurf3.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} minsurf3.o

minsurf4: minsurf3.o minsurff.ad.o  tao_chkopts
	-${CLINKER} -o minsurf3 minsurf3.o minsurff.ad.o ${TAO_ADIC_LIB} ${PETSC_SNES_LIB}
	${RM} minsurf3.o minsurff.ad.o

minsurf5: minsurf4.o  tao_chkopts
	-${CLINKER} -o minsurf4 minsurf4.o  ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} minsurf4.o 

jbearing3: jbearing3.o  tao_chkopts
	-${CLINKER} -o jbearing3 jbearing3.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} jbearing3.o

jbearing3f: jbearing3f.o tao_chkopts
	-${FLINKER} -o jbearing3f jbearing3f.o ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	${RM} jbearing3f.o

jbearing4: jbearing3.o jbearingf.ad.o  tao_chkopts
	-${CLINKER} -o jbearing3 jbearing3.o jbearingf.ad.o ${TAO_ADIC_LIB} ${PETSC_SNES_LIB}
	${RM} jbearing3.o jbearingf.ad.o

combustion3: combustion3.o  tao_chkopts
	-${CLINKER} -o combustion3 combustion3.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} combustion3.o

combustion4: combustion3.o combustionf.ad.o tao_chkopts
	-${CLINKER} -o combustion3 combustion3.o combustionf.ad.o ${TAO_ADIC_LIB} ${PETSC_SNES_LIB}
	${RM} combustion3.o combustionf.ad.o

eptorsion3: eptorsion3.o  tao_chkopts
	-${CLINKER} -o eptorsion3 eptorsion3.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} eptorsion3.o

eptorsion4: eptorsion3.o eptorsionf.ad.o  tao_chkopts
	-${CLINKER} -o eptorsion3 eptorsion3.o eptorsionf.ad.o ${TAO_ADIC_LIB} ${PETSC_SNES_LIB}
	${RM} eptorsion3.o eptorsionf.ad.o

ginzland3: ginzland3.o tao_chkopts
	-${CLINKER} -o ginzland3 ginzland3.o ${TAO_LIB} ${PETSC_SNES_LIB}
	${RM} ginzland3.o

ginzland4: ginzland3.o ginzlandf.ad.o  tao_chkopts
	-${CLINKER} -o ginzland3 ginzland3.o ginzlandf.ad.o ${TAO_ADIC_LIB} ${PETSC_SNES_LIB}
	${RM} ginzland3.o ginzlandf.ad.o

runeptorsion3:
	-@${MPIEXEC} -np 3 ./eptorsion3 -tao_smonitor \
	   -tao_method tao_gpcg -byelement > eptorsion3_1.tmp 2>&1;	\
	   if (diff output/eptorsion3_1.out eptorsion3_1.tmp) then true; \
	   else echo "Possible problem with eptorsion3_1, diffs above"; fi; \
	   ${RM} -f eptorsion3_1.tmp

runeptorsion3_2:
	-@${MPIEXEC} -np 2 ./eptorsion3 -tao_smonitor -mx 21 -my 21  -tao_fd_coloring \
	   -tao_method tao_bnls -ksp_type cg -tao_da_multigrid > eptorsion3_2.tmp 2>&1;	\
	   if (diff output/eptorsion3_2.out eptorsion3_2.tmp) then true; \
	   else echo "Possible problem with eptorsion3_2, diffs above"; fi; \
	   ${RM} -f eptorsion3_2.tmp

runjbearing3:
	-@${MPIEXEC} -np 2 ./jbearing3 -tao_smonitor  -tao_fd_coloring  -byelement\
	   -tao_method tao_gpcg -byelement -ksp_type cg > jbearing3_1.tmp 2>&1;	\
	   if (diff output/jbearing3_1.out jbearing3_1.tmp) then true; \
	   else echo "Possible problem with jbearing3_1, diffs above"; fi; \
	   ${RM} -f jbearing3_1.tmp

runjbearing3_2:
	-@${MPIEXEC} -np 3 ./jbearing3 -tao_smonitor  -mx 5 -my 15 \
	   -tao_method tao_blmvm -tao_frtol 1.0e-4 -nlevels 6 > jbearing3_2.tmp 2>&1;	\
	   if (diff output/jbearing3_2.out jbearing3_2.tmp) then true; \
	   else echo "Possible problem with jbearing3_2, diffs above"; fi; \
	   ${RM} -f jbearing3_2.tmp

runjbearing3f:
	-@${MPIEXEC} -np 2 ./jbearing3f -tao_smonitor  -tao_fd_coloring  -byelement\
	   -tao_method tao_gpcg -byelement -ksp_type cg > jbearing3f_1.tmp 2>&1;	\
	   if (diff output/jbearing3f_1.out jbearing3f_1.tmp) then true; \
	   else echo "Possible problem with jbearing3f_1, diffs above"; fi; \
	   ${RM} -f jbearing3f_1.tmp

runjbearing3f_2:
	-@${MPIEXEC} -np 3 ./jbearing3f -tao_smonitor  -mx 5 -my 15 \
	   -tao_method tao_blmvm -tao_frtol 1.0e-4 -nlevels 6 > jbearing3f_2.tmp 2>&1;	\
	   if (diff output/jbearing3f_2.out jbearing3f_2.tmp) then true; \
	   else echo "Possible problem with jbearing3f_2, diffs above"; fi; \
	   ${RM} -f jbearing3f_2.tmp

runcombustion3:
	-@${MPIEXEC} -np 2 ./combustion3 -tao_smonitor -mx 4 -my 21  \
	   -tao_method tao_tron -byelement  -tao_fd_coloring > combustion3_1.tmp 2>&1;	\
	   if (diff output/combustion3_1.out combustion3_1.tmp) then true; \
	   else echo "Possible problem with combustion3_1, diffs above"; fi; \
	   ${RM} -f combustion3_1.tmp

runcombustion3_2:
	-@${MPIEXEC} -np 3 ./combustion3 -tao_smonitor -mx 3 -my 9 \
	   -tao_method tao_bnls -tao_frtol 1.0e-4 -nlevels 6 > combustion3_2.tmp 2>&1;	\
	   if (diff output/combustion3_2.out combustion3_2.tmp) then true; \
	   else echo "Possible problem with combustion3_2, diffs above"; fi; \
	   ${RM} -f combustion3_2.tmp

runminsurf3:
	-@${MPIEXEC} -np 2 ./minsurf3 -tao_smonitor -cops -byelement -mx 14 -my 5 \
	   -tao_method tao_tron -byelement > minsurf3_1.tmp 2>&1;	\
	   if (diff output/minsurf3_1.out minsurf3_1.tmp) then true; \
	   else echo "Possible problem with minsurf3_1, diffs above"; fi; \
	   ${RM} -f minsurf3_1.tmp

runminsurf3_2:
	-@${MPIEXEC} -np 3 ./minsurf3 -tao_smonitor -mx 14 -my 11 -byelement\
	   -tao_method tao_bnls -tao_fd_coloring -nlevels 4 > minsurf3_2.tmp 2>&1;	\
	   if (diff output/minsurf3_2.out minsurf3_2.tmp) then true; \
	   else echo "Possible problem with minsurf3_2, diffs above"; fi; \
	   ${RM} -f minsurf3_2.tmp

runminsurf3_3:
	-@${MPIEXEC} -np 4 ./minsurf3 -tao_smonitor  -mx 7 -my 6\
	   -tao_method tao_blmvm -tao_frtol 1.0e-4 -nlevels 6 > minsurf3_3.tmp 2>&1;	\
	   if (diff output/minsurf3_3.out minsurf3_3.tmp) then true; \
	   else echo "Possible problem with minsurf3_3, diffs above"; fi; \
	   ${RM} -f minsurf3_3.tmp


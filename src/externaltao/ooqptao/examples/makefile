include ${TAO_DIR}/bmake/packages
#CFLAGS   = ${OOQP_INCLUDE} ${OOQP_DEFINES}

FFLAGS   =
SOURCEC  = qpssolver.c
SOURCEF  = 
SOURCEH  = 
OBJSC    = 
OBJSF    = 
LIBBASE  = 
MANSEC   = 
LOCDIR   = src/externaltao/ooqptao/examples
#LIBS = ${OOQP_LIB} 

include ${TAO_DIR}/bmake/tao_common



qpssolver: qpssolver.c tao_chkopts 
	-${CLINKER} ${OOQP_INCLUDE} \
        -o qpssolver  qpssolver.c \
	${TAO_OOQP_LIB} ${TAO_LIB} ${PETSC_SNES_LIB} 


veryclean: clean
	-rm -f qpssolver

#qpssolver.o TaoOoqpGen.o : TaoOoqpGen.h
#TaoOoqpGen.o : QpGenTao.h
#qpssolver.o TaoMpsReader.o : TaoMpsReader.h

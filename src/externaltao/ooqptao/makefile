ALL: lib

include ${TAO_DIR}/bmake/packages

CFLAGS   = ${OOQP_INCLUDE} ${PETSC_INCLUDE}
CPPFLAGS = 
FFLAGS   =
SOURCEC  = QpGenTao.c  TaoOoqpGen.c QpGenTaoLinsys.c TaoMpsReader.c
SOURCEF  = 
SOURCEH  = QpGenTao.h  TaoOoqpGen.h QpGenTaoLinsys.h TaoMpsReader.h
OBJSC    = QpGenTao.o  TaoOoqpGen.o QpGenTaoLinsys.o TaoMpsReader.o
OBJSF    = 
LIBBASE  = libooqptao
MANSEC   = 
LOCDIR   = src/externaltao/ooqptao/
LIBS = 

include ${TAO_DIR}/bmake/tao_common


veryclean: clean
	-rm -f qpssolver



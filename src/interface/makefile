include ${TAO_DIR}/bmake/packages
ALL: lib

CFLAGS   =
FFLAGS   =
CPPFLAGS = -DTAO_LIB_DIR='"${TAO_LIB_DIR}"'
SOURCEC  = tao_init.c tao.c tao_util.c tao_reg.c tao_fghj.c\
	   tao_lsolver.c dlregis.c line.c newsolver.c fdtest.c
SOURCEF  =
SOURCEH  = 
OBJSC    = tao_init.o tao.o tao_util.o tao_reg.o tao_fghj.o\
	   tao_lsolver.o dlregis.o line.o newsolver.o fdtest.o
OBJSF    = 
LIBBASE  = libtao
MANSEC   = solver
DIRS     = 
EDIRS    = 
LOCDIR   = src/interface/

include ${TAO_DIR}/bmake/tao_common


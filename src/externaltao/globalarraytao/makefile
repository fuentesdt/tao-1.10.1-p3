ALL: lib

#GA_ROOT must be defined in bmake/packages.${PETSC_ARCH} file
include ${TAO_DIR}/bmake/packages

CFLAGS   = -I${GA_ROOT}/include/ 
CPPFLAGS = ${GA_INCLUDE}
FFLAGS   =
SOURCEC  = taoapp_ga.c taomat_ga.c taovec_ga.c taolinearsolver_ga.c
SOURCEF  = 
SOURCEH  = taoapp_ga.h taomat_ga.h taovec_ga.h taolinearsolver_ga.h
OBJSC    = taoapp_ga.o taomat_ga.o taovec_ga.o taolinearsolver_ga.o
OBJSF    = 
LIBBASE  = libtaoga
MANSEC   = gaapp
LOCDIR   = src/external/globalarraytao/


include ${TAO_DIR}/bmake/tao_common

#Uncomment to avoid need for PETSc headers
#PETSC_CCPPFLAGS=-DTAO_NO_PETSC ${MPI_INCLUDE}


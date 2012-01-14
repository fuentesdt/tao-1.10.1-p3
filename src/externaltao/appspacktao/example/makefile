include ${TAO_DIR}/bmake/packages


SOURCEC = apps_driver_paraboloid.c
CFLAGS   = ${APPSPACK_INCLUDE} ${APPSPACK_DEFINES} -I${TAO_DIR}/src/externaltao/appspacktao
LOCDIR   = src/externaltao/appspacktao/example

include ${TAO_DIR}/bmake/tao_common


apps_driver_paraboloid: tao_chkopts apps_driver_paraboloid.o 
	-${CLINKER}   -o paraboloid apps_driver_paraboloid.o  \
                 ${TAO_LIB} ${TAO_APPSPACK_LIB} ${PETSC_SNES_LIB} 
	${RM} *.o










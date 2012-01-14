ALL: rosenbrockcc rosenbrockc rosenbrockf

FFLAGS=$(TAO_SIDL_FLAGS_F77)
FPPFLAGS=
EXAMPLESC=rosenbrockc
EXAMPLESCC=rosenbrockcc
EXAMPLESF=rosenbrockf
LOCDIR= src/sidl/examples/
DIRS=
CPPFLAGS =  $(TAO_SIDL_FLAGS_CXX)
CFLAGS =    $(TAO_SIDL_FLAGS_C)

CLEANFILES=rosenbrockc rosenbrockcc rosenbrockf

include ${TAO_DIR}/bmake/tao_common

rosenbrockcc: rosenbrockcc.o tao_chkopts
	-${CLINKER} -o rosenbrockcc rosenbrockcc.o  -L${TAO_LIB_DIR} ${TAO_SIDLCLIENT_CXX_LIB} ${PETSC_LIB}
	${RM} rosenbrockcc.o

rosenbrockc: rosenbrockc.o tao_chkopts
	-${CLINKER} -o rosenbrockc rosenbrockc.o  -L${TAO_LIB_DIR} ${TAO_SIDLCLIENT_C_LIB} ${TAO_SIDLCLIENT_CXX_LIB} ${PETSC_LIB} 
	${RM} rosenbrockc.o

rosenbrockf: rosenbrockf.o tao_chkopts
	-${FLINKER} -o rosenbrockf rosenbrockf.o  -L${TAO_LIB_DIR} ${TAO_SIDLCLIENT_F77_LIB} ${TAO_SIDLCLIENT_CXX_LIB} ${PETSC_LIB}

test: test.o tao_chkopts
	-${CLINKER} -o test test.o  ${SIDL_LIB} \
          ${CLINKER_SLFLAG}${TAO_LIB_DIR} -L${TAO_LIB_DIR} -lOptimize-client-C++ -lsidl-client-C++
	${RM} test.o

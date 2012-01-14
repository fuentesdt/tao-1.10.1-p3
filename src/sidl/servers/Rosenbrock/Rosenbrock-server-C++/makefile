ALL: lib
lib: libRosenbrock-server-C++.so
libfast: libRosenbrock-server-C++.so
#get source and header filenames for this sidl package
-include babel.make

#LIBBASE=libRosenbrock-server-C++
SOURCEC=$(IMPLSRCS) $(STUBSRCS) $(SKELSRCS) $(IORSRCS)
SOURCEH=$(IMPLHDRS) $(IORHDRS) $(STUBHDRS)
CFLAGS=-I$(BABEL_HOME)/include -I$(TAO_DIR)/src/sidl/clients/sidl/sidl-client-C++
OBJSC = $(IMPLSRCS:.cxx=.o) $(IORSRCS:.c=.o) $(SKELSRCS:.cxx=.o) $(STUBSRCS:.cxx=.o)
CLEANFILES=$(STUBSRCS) $(SKELSRCS) $(IORSRCS) $(IORHDRS) $(STUBHDRS) babel.make

CFLAGS = ${TAO_SIDL_FLAGS_CXX} \
	-I$(CCA_SPEC_BABEL_INCLUDE) \
        -I$(CCA_H_DIR)

MYLIBS =  ${TAO_C_SH_LIB_PATH}  ${TAO_SIDLCLIENT_CXX_LIB} \
        ${TAO_LIB} ${PETSC_C_SH_LIB_PATH} ${PETSC_LIB}

libRosenbrock-server-C++.so: ${OBJSC} tao_chkopts
	@${CLINKER} -shared ${OBJSC} -o ${TAO_DIR}/lib/${PETSC_ARCH}/$@ ${MYLIBS} ${MYLIBS}


include $(TAO_DIR)/bmake/tao_common

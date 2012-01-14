ALL: libTaoapi-server-C++.so
libfast: libTaoapi-server-C++.so
lib: libTaoapi-server-C++.so
-include babel.make


#LIBBASE=libTaoapi-server-C++
SOURCEC=$(IMPLSRCS) $(STUBSRCS) $(SKELSRCS) $(IORSRCS)
SOURCEH=$(IMPLHDRS) $(IORHDRS) $(STUBHDRS)
OBJSC = $(IMPLSRCS:.cxx=.o) $(IORSRCS:.c=.o) $(SKELSRCS:.cxx=.o) $(STUBSRCS:.cxx=.o)
CLEANFILES=$(STUBSRCS) $(SKELSRCS) $(IORSRCS) $(IORHDRS) $(STUBHDRS) babel.make




CFLAGS = $(TAO_INCLUDE) $(PETSC_INCLUDE) -I$(BABEL_HOME)/include \
          -I$(TAO_DIR)/src/externaltao/petsctao \
	-I$(TAO_DIR)/src/sidl/interface \
	  $(TAO_SIDL_FLAGS_CXX)



LIBS = ${PETSC_C_SH_LIB_PATH} ${TAO_C_SH_LIB_PATH} ${TAO_LIB} ${SIDL_LIB} \
	${TAO_SIDLCLIENT_CXX_LIB} ${TAO_SIDL_SERVER_LIB} ${PETSC_LIB}


libTaoapi-server-C++.so: ${OBJSC} tao_chkopts
	@${CLINKER} -shared ${OBJSC} -o ${TAO_DIR}/lib/${PETSC_ARCH}/$@ ${LIBS}



include $(TAO_DIR)/bmake/tao_common


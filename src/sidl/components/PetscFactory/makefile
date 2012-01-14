ALL: libfast
lib: libPetscFactory.so
libfast: libPetscFactory.so
-include babel.make


LIBBASE=libPetscFactory
SOURCEC=$(IMPLSRCS) $(STUBSRCS) $(SKELSRCS) $(IORSRCS)
SOURCEH=$(IMPLHDRS) $(IORHDRS) $(STUBHDRS)
OBJSC = $(IMPLSRCS:.cxx=.o) $(IORSRCS:.c=.o) $(SKELSRCS:.cxx=.o) $(STUBSRCS:.cxx=.o)
CLEANFILES=$(STUBSRCS) $(SKELSRCS) $(IORSRCS) $(IORHDRS) $(STUBHDRS) babel.make

CFLAGS = $(TAO_SIDL_FLAGS_CXX) $(PETSC_INCLUDE) $(TAO_INCLUDE)

LIBS = ${TAO_SIDLCLIENT_CXX_LIB} \
	${TAO_C_SH_LIB_PATH} ${TAO_LIB} \
    ${CC_LINKER_SLFLAG}${CCASO_DIR}  -L${CCASO_DIR} -lcca ${SIDL_LIB} ${PETSC_LIB} 


libPetscFactory.so: ${OBJSC} tao_chkopts
	@${CLINKER} -shared ${OBJSC} -o ${TAO_DIR}/lib/${PETSC_ARCH}/$@ ${LIBS}

	@sed s^@TAO_DIR@^$(TAO_DIR)^g ../cca_template | \
           sed s^@PETSC_ARCH@^$(PETSC_ARCH)^g |\
           sed s^@COMPONENT@^PetscLA.Factory^g |\
           sed s^@LIBNAME@^$@^g \
	   > $(TAO_DIR)/lib/$(PETSC_ARCH)/$@.cca



include $(TAO_DIR)/bmake/tao_common


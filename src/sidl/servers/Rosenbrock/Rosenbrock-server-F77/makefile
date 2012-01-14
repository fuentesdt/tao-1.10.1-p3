ALL: lib
-include babel.make

LIBBASE=libRosenbrock-server-F77
SOURCEC=$(STUBSRCS) $(SKELSRCS) $(IORSRCS)
SOURCEF=$(IMPLSRCS)
SOURCEH=$(IORHDRS) $(STUBHDRS)
CFLAGS=-I$(BABEL_HOME)/include
OBJSC=$(IMPLSRCS:.f=.o) $(IORSRCS:.c=.o) $(SKELSRC:.c=.o) $(STUBSRCS:.c=.o)
CLEANFILES=$(STUBSRCS) $(SKELSRCS) $(IORSRCS) $(IORSRCS) $(IORHDRS) $(STUBHDRS) babel.make



include $(TAO_DIR)/bmake/tao_common

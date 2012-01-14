ALL: lib
-include babel.make

LIBBASE=libRosenbrock-server-C
SOURCEC=$(IMPLSRCS) $(STUBSRCS) $(SKELSRCS) $(IORSRCS)
SOURCEH=$(IMPLHDRS) $(IORHDRS) $(STUBHDRS)
CFLAGS=-I$(BABEL_HOME)/include
OBJSC=$(IMPLSRCS:.c=.o) $(IORSRCS:.c=.o) $(SKELSRC:.c=.o) $(STUBSRCS:.c=.o)
CLEANFILES=$(STUBSRCS) $(SKELSRCS) $(IORSRCS) $(IORSRCS) $(IORHDRS) $(STUBHDRS) babel.make


clean:
	$(RM) $(BABEL_OBJS) $(LIBPATH)/libRosenbrock-server-C.so

cleanstubs:
	$(RM) $(IORSRCS) $(IORHDRS) $(SKELSRCS) $(STUBHDRS) $(STUBSRCS) babel.make



include $(TAO_DIR)/bmake/tao_common

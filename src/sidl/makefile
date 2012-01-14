DIRS = clients interface servers ${TAOCCAFFEINEDIR}

LOCDIR=src/sidl
CLEANFILES = xml/* .xml

VPATH = sidlfiles

# Order of sidl files is important!!
SIDLFILES = Optimize.sidl LAFactory.sidl Solver.sidl Taoapi.sidl Rosenbrock.sidl 

ALL: lib
lib: xml
libfast: xml

xml: .xml
.xml: $(SIDLFILES)
	@cd sidlfiles && $(BABEL) --text=xml -o ../xml $(SIDLFILES) > /dev/null
	@touch .xml



include ${TAO_DIR}/bmake/tao_common

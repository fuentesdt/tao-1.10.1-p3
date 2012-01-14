#!/bin/sh
all: lib
lib: xml  ports componentstubs  CCAFERC CCAFERC_gui sclfiles
libfast: xml  ports componentstubs  CCAFERC CCAFERC_gui sclfiles

PORTS = $(foreach l,$(SIDL_CLIENT_LANGUAGES),OptimizationPorts-$l) 
COMPONENTS = TaoDriver TaoSolver Minsurf PetscFactory $(GAFACTORY_SIDL:.sidl=)
DIRS = $(PORTS) 
#$(COMPONENTS)
VPATH = sidlfiles
SIDLFILES = OptimizationPorts.sidl Minsurf.sidl TaoDriver.sidl TaoSolver.sidl PetscFactory.sidl $(GAFACTORY_SIDL)
REPOSITORY = -R $(TAO_DIR)/src/sidl/xml -R $(CCA_XML_DIR)
EXCLUDES=--exclude="^sidl.*" --exclude="^gov.*" --exclude="^Optimize.*"\
         --exclude="^Solver.*" --exclude="^LinearAlgebra.*"

CLEANFILES = .xml  .stubs .ports .templates .componentstubs .sclfiles

xml: .xml
.xml: $(SIDLFILES) ../.xml
	@$(RM) -f ../xml/*.xml ../.xml
	@cd .. && $(MAKE) xml
	@cd sidlfiles && $(BABEL) --text=xml -o ../../xml $(SIDLFILES) $(REPOSITORY) > /dev/null
	@touch .xml



ports: .ports templates
.ports: .xml  
	@$(foreach l, $(SIDL_CLIENT_LANGUAGES), \
            if [ "$lx" = "Pythonx" ]; then \
            $(BABEL)  -c$(l) -oOptimizationPorts-$(l)/ OptimizationPorts.SolverPort OptimizationPorts.ModelPort OptimizationPorts.LAFactoryPort  $(REPOSITORY) > /dev/null; \
	    else  \
            $(BABEL) $(EXCLUDES)  -c$(l) -oOptimizationPorts-$(l)  OptimizationPorts.SolverPort OptimizationPorts.ModelPort OptimizationPorts.LAFactoryPort $(REPOSITORY) > /dev/null; \
            fi;)

	@touch .ports


templates: .templates 
.templates:
	@$(foreach l, $(SIDL_CLIENT_LANGUAGES), \
	 mkdir -p OptimizationPorts-$l; \
         sed "s^@SIDL_PACKAGE_NAME@^OptimizationPorts^g"  \
         $(TAO_DIR)/src/sidl/clients/makefile_templates/makefile-$(l)\
         > $(TAO_DIR)/src/sidl/components/OptimizationPorts-$(l)/makefile;)
	-@sed "s^OptimizationPorts/^^" OptimizationPorts-Python/makefile > tempMK
	-@sed "s^-client^^" tempMK  > OptimizationPorts-Python/makefile
	-@${RM} tempMK

	@touch .templates

componentstubs: .componentstubs
.componentstubs: .xml 
	@$(BABEL) $(EXCLUDES) --exclude="^OptimizationPorts.*" -sC++ -oMinsurf Minsurf.Component $(REPOSITORY)> /dev/null
	@$(BABEL) $(EXCLUDES) --exclude="^OptimizationPorts.*" -sC++ -oTaoDriver  TaoDriver.Component $(REPOSITORY) > /dev/null
	@$(BABEL) $(EXCLUDES) --exclude="^OptimizationPorts.*" -sC++ -oTaoSolver  TaoSolver.Component $(REPOSITORY) > /dev/null
	@$(BABEL) $(EXCLUDES) --exclude="^OptimizationPorts.*" -sC++ -oPetscFactory  PetscLA.Factory $(REPOSITORY) > /dev/null
#	@$(BABEL) $(EXCLUDES) --exclude="^OptimizationPorts.*" -sC++ -oGAFactory  GA.Factory $(REPOSITORY) > /dev/null

	-@${foreach i, ./TaoSolver/TaoSolver_Component_IOR.c ./PetscFactory/PetscLA_Factory_IOR.c ./Minsurf/Minsurf_Component_IOR.c ./TaoDriver/TaoDriver_Component_IOR.c, \
	  if [ "`grep c++-enable $i`foo" = foo ] ; then  \
	    mv $i tmp; \
	    echo '/*c++-enable*/' > $i; \
	    echo '#ifdef __cplusplus' >> $i;\
	    echo 'extern "C" {' >> $i;\
	    echo '#endif' >> $i;\
	    cat tmp >> $i; \
	    echo '#ifdef __cplusplus' >> $i;\
	    echo '}' >> $i;\
	    echo '#endif' >> $i;\
	    rm -f tmp;\
	  fi; \
        }


	@touch .componentstubs

components:
	@$(foreach c,${COMPONENTS}, \
	  cd ${TAO_DIR}/src/sidl/components/$c && ${OMAKE};)



CCAFERC: CCAFERC.in tao_chkopts
	-@cat $(TAO_DIR)/src/sidl/components/CCAFERC.in | \
            sed s^@TAO_DIR@^$(TAO_DIR)^g |\
            sed s^@PETSC_ARCH@^$(PETSC_ARCH)^g \
            > $(TAO_DIR)/src/sidl/components/CCAFERC 


#cca: .cca
#.cca:
#	@$(foreach c,${COMPONENTS}, \
#	  cat $(TAO_DIR)/src/sidl/components/components.cca | \
#            sed s^@TAO_DIR@^$(TAO_DIR)^g |\
#            sed s^@PETSC_ARCH@^$(PETSC_ARCH)^g | \
#            sed s^@COMPONENT@^$(c)^g | \
#            sed s^@LIBNAME@^$(c)^g \
#	 > $(TAO_DIR)/lib/${PETSC_ARCH}/components.cca;)


sclfiles: .sclfiles
.sclfiles:



CCAFERC_gui: CCAFERC_gui.in
	-@cat $(TAO_DIR)/src/sidl/components/CCAFERC_gui.in | \
            sed s^@TAO_DIR@^$(TAO_DIR)^g |\
            sed s^@PETSC_ARCH@^$(PETSC_ARCH)^g \
            > $(TAO_DIR)/src/sidl/components/CCAFERC_gui

include ${TAO_DIR}/bmake/tao_common

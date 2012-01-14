ALL: lib
libfast:  serverstubs sclfiles
lib:  serverstubs sclfiles

DIRS = Rosenbrock Taoapi 
LOCDIR=src/sidl/servers/
EXCLUDES=--exclude="^sidl.*" --exclude="^gov.*" --exclude="^Solver.*" --exclude="^Optimize.*" --exclude="^LinearAlgebra.*"
VPATH=${TAO_DIR}/src/sidl/sidlfiles
CLEANFILES=.serverstubs  .sclfile

serverstubs:  $(STUBS) .serverstubs chkopts
.serverstubs: Optimize.sidl Solver.sidl LAFactory.sidl
	@$(BABEL) -sC++ $(EXCLUDES)\
          -o $(TAO_DIR)/src/sidl/servers/Rosenbrock/Rosenbrock-server-C++ \
          -R $(TAO_DIR)/src/sidl/xml Rosenbrock.RosenbrockModel > /dev/null

	@$(BABEL)  -sC++ $(EXCLUDES)\
         -o$(TAO_DIR)/src/sidl/servers/Taoapi/Taoapi-server-C++ \
          -R $(TAO_DIR)/src/sidl/xml TAO.Solver TAO.Environment > /dev/null

	-@${foreach i, ./Rosenbrock/Rosenbrock-server-C++/Rosenbrock_RosenbrockModel_IOR.c ./Taoapi/Taoapi-server-C++/TAO_Solver_IOR.c ./Taoapi/Taoapi-server-C++/TAO_Environment_IOR.c, \
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


	@touch .serverstubs

sclfiles:
	@cat Taoapi/Taoapi-server-C++/libTaoapi-server-C++.so.scl | \
	 sed s^@TAO_DIR@^$(TAO_DIR)^g | \
         sed s^@PETSC_ARCH@^$(PETSC_ARCH)^g > \
         $(TAO_DIR)/lib/$(PETSC_ARCH)/libTaoapi-server-C++.so.scl

	@cat Rosenbrock/Rosenbrock-server-C++/libRosenbrock-server-C++.so.scl|\
	 sed s^@TAO_DIR@^$(TAO_DIR)^g | \
         sed s^@PETSC_ARCH@^$(PETSC_ARCH)^g > \
         $(TAO_DIR)/lib/$(PETSC_ARCH)/libRosenbrock-server-C++.so.scl


include ${TAO_DIR}/bmake/tao_common


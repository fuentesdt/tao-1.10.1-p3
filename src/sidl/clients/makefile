
ALL: lib
lib: clientstubs templates
libfast: clientstubs templates

LOCDIR=src/sidl/clients/
VPATH = $(TAO_DIR)/src/sidl/sidlfiles
DIRS=Optimize Solver  LAFactory 
CLEANFILES=.clientstubs .templates include 

EXCLUDES=--exclude="^sidl.*" --exclude="^gov.*"


clientstubs: .clientstubs tao_chkopts
.clientstubs: Optimize.sidl Solver.sidl LAFactory.sidl 
	@mkdir -p $(DIRS)
#	@$(foreach i, $(SIDL_CLIENT_LANGUAGES), \
#          $(BABEL) -c$(i) -o ./sidl/sidl-client-$(i)/ \
#          sidl -R $(TAO_DIR)/src/sidl/xml > /dev/null;)


	@$(foreach i,$(SIDL_CLIENT_LANGUAGES), \
           $(BABEL) -c$(i) -o ./Optimize/Optimize-client-$(i)/ $(EXCLUDES) \
           Optimize.OptimizationModel -R $(TAO_DIR)/src/sidl/xml > /dev/null;)

	@$(foreach i,$(SIDL_CLIENT_LANGUAGES), \
           $(BABEL) -c$(i) -o ./Solver/Solver-client-$(i)/ $(EXCLUDES) \
           Solver.OptimizationSolver Solver.ProjectState Solver.TerminateReason \
           -R $(TAO_DIR)/src/sidl/xml  > /dev/null;)

	@$(foreach i,$(SIDL_CLIENT_LANGUAGES), \
	   $(BABEL) -c$(i) -o ./LAFactory/LAFactory-client-$(i)/ $(EXCLUDES) \
           LinearAlgebra.Factory -R $(TAO_DIR)/src/sidl/xml  > /dev/null;)

	@touch .clientstubs

templates: .templates clientstubs tao_chkopts
.templates: 
	@$(foreach i, $(DIRS),  \
          $(foreach l, $(SIDL_CLIENT_LANGUAGES), \
             mkdir -p $(TAO_DIR)/src/sidl/clients/$(i)/$(i)-client-$(l);\
             sed "s^@SIDL_PACKAGE_NAME@^$(i)^g" \
               $(TAO_DIR)/src/sidl/clients/makefile_templates/makefile-$(l)\
               > $(TAO_DIR)/src/sidl/clients/$(i)/$(i)-client-$(l)/makefile;)\
	  rm -f $i/makefile;\
	  echo "ALL: lib" > $i/makefile;\
	  echo "LOCDIR=src/sidl/clients/$i" >> $i/makefile;\
          echo "DIRS=$(foreach l, $(SIDL_CLIENT_LANGUAGES),$i-client-$l)" >> $i/makefile;\
	  echo 'include $${TAO_DIR}/bmake/tao_common' >> $i/makefile;)

	  @touch .templates
# Make client stubs



include ${TAO_DIR}/bmake/tao_common



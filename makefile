# $Id$ 
#
# This is the makefile for installing TAO. See the file
# docs/installation.html for directions on installing TAO.
# See also bmake/common for additional commands.
#
ALL: all

DIRS	   = src include docs 

include ${TAO_DIR}/bmake/packages
include ${TAO_DIR}/bmake/tao_common

#
# Basic targets to build TAO libraries.
# all     : builds the C/C++ and Fortran libraries
all       : info tao_chkcxx chktao_dir tao_chklib_dir tao_deletelibs tao_build_c tao_build_fortran tao_shared 
#
# Prints information about the system and version of TAO being compiled
#
info:
	-@echo "=========================================="
	-@echo " "
	-@echo "See docs/troubleshooting.html and docs/bugreporting.html"
	-@echo "for help with installation problems. Please send EVERYTHING"
	-@echo "printed out below when reporting problems."
	-@echo " "
	-@echo "To subscribe to the TAO users mailing list, send mail to "
	-@echo "majordomo@mcs.anl.gov with the message: "
	-@echo "subscribe tao-news"
	-@echo " "
	-@echo "=========================================="
	-@echo On `date` on `hostname`
	-@echo Machine characteristics: `uname -a`
	-@echo "-----------------------------------------"
	-@echo "Using PETSc directory: ${PETSC_DIR}"
	-@echo "Using PETSc arch: ${PETSC_ARCH}"
	-@echo "Using TAO directory: ${TAO_DIR}"
	-@echo "-----------------------------------------"
	-@grep "define PETSC_VERSION" ${PETSC_DIR}/include/petscversion.h | ${SED} "s/........//"
	-@grep TAO_VERSION_NUMBER include/tao_version.h | sed "s/........//"
	-@echo "-----------------------------------------"
	-@echo "Using include paths: ${TAO_INCLUDE}"
	-@echo "------------------------------------------"
	-@echo "Using C/C++ compiler: ${CC} ${CC_FLAGS} ${COPTFLAGS} ${CFLAGS}"
	-@echo "C/C++ Compiler version: " `${CCV}`
	-@if [ "${FC}" != "" ]; then \
	   echo "Using Fortran compiler: ${FC} ${FC_FLAGS} ${FFLAGS} ${FPP_FLAGS}";\
	   echo "Fortran Compiler version: " `${FCV}`;\
         fi
	-@echo "-----------------------------------------"
	-@echo "Using C/C++ linker: ${CC_LINKER}"
	-@if [ "${FC}" != "" ]; then \
	   echo "Using Fortran linker: ${FC_LINKER}";\
         fi
	-@echo "-----------------------------------------"
	-@echo "Using libraries: ${TAO_LIB}"
	-@echo "------------------------------------------"
	-@echo "Using mpirun: ${MPIEXEC}"
	-@echo "=========================================="


#
# Builds the TAO libraries
# This target also builds fortran77 and f90 interface
# files (except compiling *.F files).
#
tao_build_c:
	-@echo "BEGINNING TO COMPILE TAO LIBRARIES IN ALL DIRECTORIES"
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH} ACTION=libfast  tree 
	${RANLIB} ${TAO_LIB_DIR}/*.${AR_LIB_SUFFIX}
	-@echo "Completed building TAO libraries"
	-@echo "========================================="

#
# Builds TAO Fortran source
# Note: tao_libfast cannot run on .F files on certain machines, so we
# use libf to compile the Fortran source files.
#
tao_build_fortran:
	-@echo "BEGINNING TO COMPILE TAO FORTRAN SOURCE"
	-@echo "========================================="
	-@cd src/fortran/custom; \
	  ${OMAKE} PETSC_ARCH=${PETSC_ARCH} libf clean 
	-${RANLIB} ${TAO_LIB_DIR}/libtaofortran.${AR_LIB_SUFFIX}
	${RANLIB} ${TAO_LIB_DIR}/libtao.${AR_LIB_SUFFIX}
	-@echo "Completed compiling TAO Fortran source"
	-@echo "========================================="

# Builds TAO test examples for a given architecture
tao_testexamples: info chkopts
	-@echo "BEGINNING TO COMPILE AND RUN TAO TEST EXAMPLES"
	-@echo "Due to different numerical round-off on certain"
	-@echo "machines some of the numbers may not match exactly."
	-@echo "========================================="
	-@${OMAKE}  PETSC_ARCH=${PETSC_ARCH} \
	   ACTION=testexamples_C  tree 
	-@echo "Completed compiling and running test examples"
	-@echo "========================================="

# Builds TAO test examples for a given architecture
tao_testfortran: info chkopts
	-@echo "BEGINNING TO COMPILE AND RUN TAO FORTRAN TEST EXAMPLES"
	-@echo "========================================="
	-@echo "Due to different numerical round-off on certain"
	-@echo "machines or the way Fortran formats numbers"
	-@echo "some of the results may not match exactly."
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH} \
	   ACTION=testexamples_Fortran  tree 
	-@echo "Completed compiling and running Fortran test examples"
	-@echo "========================================="

# Builds TAO test examples for a given architecture
tao_testexamples_uni: info chkopts
	-@echo "BEGINNING TO COMPILE AND RUN TEST UNI-PROCESSOR EXAMPLES"
	-@echo "Due to different numerical round-off on certain"
	-@echo "machines some of the numbers may not match exactly."
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH} \
	   ACTION=testexamples_C_X11_MPIUni  tree 
	-@echo "Completed compiling and running uniprocessor test examples"
	-@echo "========================================="
tao_testfortran_uni: info chkopts
	-@echo "BEGINNING TO COMPILE AND RUN TEST UNI-PROCESSOR FORTRAN EXAMPLES"
	-@echo "Due to different numerical round-off on certain"
	-@echo "machines some of the numbers may not match exactly."
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH} \
	   ACTION=testexamples_Fortran_MPIUni tree 
	-@echo "Completed compiling and running uniprocessor fortran test examples"
	-@echo "========================================="

# Ranlib on the libraries
tao_ranlib:
	${RANLIB} ${TAO_LIB_DIR}/*.${AR_LIB_SUFFIX}

# Deletes TAO libraries
tao_deletelibs:
	-${RM} -f ${TAO_LIB_DIR}/*

tao_shared: shared
	@if [ "${TAOSIDLDIR}x" = "sidlx" ]; then \
	  cd ${TAO_DIR}/src/sidl/servers/Taoapi/Taoapi-server-C++ && make; \
         fi;
	@if [ "${TAO_BUILD_COMPONENTS}x" != "x" ]; then \
	  cd ${TAO_DIR}/src/sidl/components && make components; \
         fi;

tao_chklib_dir:chklib_dir

# ------------------------------------------------------------------
#
# All remaining actions are intended for TAO developers only.
# TAO users should not generally need to use these commands.
#
# To access the tags in EMACS, type M-x visit-tags-table and specify
# the file tao/TAGS.	
# 1) To move to where a TAO function is defined, enter M-. and the
#     function name.
# 2) To search for a string and move to the first occurrence,
#     use M-x tags-search and the string.
#     To locate later occurrences, use M-,
# Builds all etags files
tao_alletags:
	-@maint/generateetags.py


BMAKEFILES = bmake/tao_common* 
DOCS	   = bmake/readme bmake/petscconf.defs
SCRIPTS    = maint/addlinks maint/builddist maint/buildlinks maint/wwwman \
	     maint/xclude maint/crontab  \
	     maint/autoftp include/foldinclude/generateincludes

# Builds all the documentation - should be done every night
tao_alldoc: tao_allmanpages
	cd docs/tex/manual; ${OMAKE} manual.dvi manual.pdf manual.html 



# Deletes man pages (HTML version)
tao_deletemanualpages:
	${RM} -f ${TAO_DIR}/docs/manualpages/*/*.html \
                 ${TAO_DIR}/docs/manualpages/manualpages.cit 

# Deletes man pages (LaTeX version)
tao_deletelatexpages:
	${RM} -f ${TAO_DIR}/docs/tex/rsum/*sum*.tex

# Builds all versions of the man pages
#tao_allmanpages: tao_allmanualpages tao_alllatexpages
tao_allmanpages: tao_allmanualpages tao_htmlpages

tao_allmanualpages: tao_deletemanualpages
	-${OMAKE} ACTION=tao_manualpages_buildcite tree
	-${OMAKE} ACTION=tao_manualpages tree
	-${OMAKE} ACTION=tao_manexamples tree LOC=${LOC}
	-maint/wwwindex.py ${TAO_DIR}
	-maint/examplesindex.tcl
	-maint/htmlkeywords.tcl

tao_htmlpages: 
	-${OMAKE} ACTION=tao_html TAO_DIR=${TAO_DIR} PETSC_DIR=${PETSC_DIR} alltree LOC=${TAO_DIR}

tao_alllatexpages: tao_deletelatexpages
	-${OMAKE} ACTION=tao_latexpages tree

# Builds Fortran stub files
tao_allfortranstubs:
	-@maint/generatefortranstubs.py ${BFORT}


# -------------------------------------------------------------------------------
#
# Some macros to check if the fortran interface is up-to-date.
#
tao_countfortranfunctions: 
	-@cd ${TAO_DIR}/src/fortran; egrep '^void' custom/*.c auto/*.c | \
	cut -d'(' -f1 | tr -s  ' ' | cut -d' ' -f3 | uniq | egrep -v "(^$$|Tao)" | \
	sed "s/_$$//" | sort > /tmp/countfortranfunctions

tao_countcfunctions:
	-@ grep extern ${TAO_DIR}/include/*.h  ${TAO_DIR}/src/petsctao/include/*.h | grep "(" | tr -s ' ' |  sed s/static// |\
	cut -d'(' -f1 | cut -d' ' -f3 | grep -v "\*" | tr -s '\012' |  \
	tr 'A-Z' 'a-z' |  sort > /tmp/countcfunctions

tao_difffortranfunctions: tao_countfortranfunctions tao_countcfunctions
	-@echo -------------- Functions missing in the Fortran interface ---------------------
	-@diff /tmp/countcfunctions /tmp/countfortranfunctions | grep "^<" | cut -d' ' -f2
	-@echo ----------------- Functions missing in the C interface ------------------------
	-@diff /tmp/countcfunctions /tmp/countfortranfunctions | grep "^>" | cut -d' ' -f2
	-@${RM}  /tmp/countcfunctions /tmp/countfortranfunctions

tao_checkbadfortranstubs:
	-@echo "========================================="
	-@echo "Functions with MPI_Comm as an Argument"
	-@echo "========================================="
	-@cd ${TAO_DIR}/src/fortran/auto; grep '^void' *.c | grep 'MPI_Comm' | \
	tr -s ' ' | tr -s ':' ' ' |cut -d'(' -f1 | cut -d' ' -f1,3
	-@echo "========================================="
	-@echo "Functions with a String as an Argument"
	-@echo "========================================="
	-@cd ${TAO_DIR}/src/fortran/auto; grep '^void' *.c | grep 'char \*' | \
	tr -s ' ' | tr -s ':' ' ' |cut -d'(' -f1 | cut -d' ' -f1,3
	-@echo "========================================="
	-@echo "Functions with Pointers to TAO Objects as Argument"
	-@echo "========================================="
	-@cd ${TAO_DIR}/src/fortran/auto; \
	_p_OBJ=`grep _p_ ${TAO_DIR}/include/*.h | tr -s ' ' | \
	cut -d' ' -f 3 | tr -s '\012' | grep -v '{' | cut -d'*' -f1 | \
	sed "s/_p_//g" | tr -s '\012 ' ' *|' ` ; \
	for OBJ in $$_p_OBJ; do \
	grep "$$OBJ \*" *.c | tr -s ' ' | tr -s ':' ' ' | \
	cut -d'(' -f1 | cut -d' ' -f1,3; \
	done 


#
# Automatically generates TAO exercises in html from the tutorial examples.
#
# The introduction for each section is obtained from docs/manualpages/bop.${MANSEC} is under RCS and may be edited
#  (used also in introductions to the manual pages)
# The overall introduction is in docs/exercises/introduction.html and is under RCS and may be edited
# The list of exercises is from TUTORIALS in each directory's makefile
#
# DO NOT EDIT the pageform.txt or *.htm files generated since they will be automatically replaced.
# The pagemaker rule is in the file bmake/common (at the bottom)
#
# Eventually the line below will replace the two cd in the rule below, it is just this way now for speed
#	-@${OMAKE} TAO_DIR=${TAO_DIR} pagemaker
#
exercises:
	-@echo "========================================="
	-@echo "Generating HTML tutorial exercises"
	-@rm -f docs/pageform.txt
	-@echo "title=\"TAO Exercises\""                >  docs/pageform.txt 
	-@echo "access_title=Exercise Sections"              >>  docs/pageform.txt 
	-@echo "access_format=short"                                            >> pageform.txt
	-@echo "startpage=../exercises/introduction.htm"  >> docs/pageform.txt
	-@echo "NONE title=\"Introduction\" command=link src=../exercises/introduction.htm" >> docs/pageform.txt
	-@echo "Generating HTML for individual directories"
	-@echo "========================================="
	cd src/unconstrained/examples/tutorials; ${OMAKE} TAO_DIR=${TAO_DIR} taopagemaker
	cd src/bound/examples/tutorials; ${OMAKE} TAO_DIR=${TAO_DIR} taopagemaker
	cd src/least_squares/examples/tutorials; ${OMAKE} TAO_DIR=${TAO_DIR} taopagemaker
	cd src/complementarity/examples/tutorials; ${OMAKE} TAO_DIR=${TAO_DIR} taopagemaker
	cd src/petsctao/gridapplication/examples; ${OMAKE} TAO_DIR=${TAO_DIR} taopagemaker
	cd src/externaltao/globalarraytao/examples; ${OMAKE} TAO_DIR=${TAO_DIR} taopagemaker
	-@echo "Completed HTML for individual directories"
	-@echo "NONE title=\"<HR>\" " >> docs/pageform.txt; 
	-@echo "NONE title=\"TAO Documentation\" command=link src=http://www-fp.mcs.anl.gov/tao/docs/index.htm" >> docs/pageform.txt
	/home/sarich/software/makepage/makepage.new -pageform=docs/pageform.txt -access_extra=/dev/null -outdir=docs/exercises 
	-@echo "========================================="


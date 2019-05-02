###############################################################################
#									      #
#       Main Makefile for CNSsolve 					      #
#       Paul Adams 21-10-98						      #
#       copyright Yale University					      #
#									      #
###############################################################################

SHELL=/bin/sh

CSHELL=/bin/csh

default:
	@ echo "options: "; \
	  echo "   install          [install]"; \
	  echo "   update           [update the installation]"; \
	  echo "   utils            [make utility programs]"; \
	  echo "   patch            [install patch from file]"; \
	  echo "   run_inputs       [run all input files]"; \
	  echo "   win32_inputs     [make run_inputs.bat for Windows-NT]"; \
	  echo "   g77install       [install using g77]"; \
	  echo "   g77update        [update the g77 installation]"; \
	  echo "   cns_solve        [recompile]"; \
	  echo "   makefiles        [remake makefiles]"; \
	  echo "   machinefile      [link machine specific files to source]"; \
	  echo "   fftfile          [link fft file to source]"; \
	  echo "   fftfile-generic  [link generic fft file to source]"; \
	  echo "   relink           [relink source files]"; \
	  echo "   syntax-manual    [make online syntax manual and index]"; \
	  echo "   compiler-test    [test C and Fortran compilers]"; \
	  echo "   space            [remove object files]"; \
	  echo "   clean            [remove object and executable]"; \
	  echo "   static-fastm     [use the static fastm Irix libraries]"; \
	  echo "   no-fastm         [do not use fastm Irix libraries]"; \
	  echo "   tidy             [remove stray files]"

install:
	@ $(CSHELL) -cf 'unsetenv CNS_ARCH; \
	                 source ./cns_solve_env; ./bin/install'

update:
	@ $(CSHELL) -cf 'unsetenv CNS_ARCH; \
	                 source ./cns_solve_env; ./bin/install -update'

g77install:
	@ $(CSHELL) -cf 'unsetenv CNS_ARCH; \
	                 source ./cns_solve_env; ./bin/install -g77'

g77update:
	@ $(CSHELL) -cf 'unsetenv CNS_ARCH; \
	                 source ./cns_solve_env; ./bin/install -update -g77'

patch:
	@ $(CSHELL) -cf 'source ./cns_solve_env; \
	                 cp ./bin/apply_patch /usr/tmp/apply_patch; \
	                                      /usr/tmp/apply_patch; \
	                 /bin/rm -f /usr/tmp/apply_patch'

utils: valid_solve valid_comp
	@ $(CSHELL) -cf 'source $(CNS_INST)/arch_env; \
	                 cd $(CNS_INST); make -f Makefile.header compile-utils'

cns_solve: valid_solve valid_comp
	@ $(CSHELL) -cf 'source $(CNS_INST)/arch_env; \
	                 cd $(CNS_SOURCE); make -k cns_solve'

compiler-test: valid_solve valid_comp
	@ $(CSHELL) -cf 'source $(CNS_INST)/arch_env; \
	                 cd $(CNS_SOURCE); make compiler-test'

fftfile: valid_solve valid_comp
	@ $(CSHELL) -cf 'source $(CNS_INST)/arch_env; cd $(CNS_SOURCE); \
	                 make -f Makefile.proto -k fftfile OBJS=none'

fftfile-generic: valid_solve valid_comp
	@ $(CSHELL) -cf 'source $(CNS_INST)/arch_env; cd $(CNS_SOURCE); \
	             make -f Makefile.proto -k fftfile CNS_FFTLIB="" OBJS=none'

no-fastm: valid_solve valid_comp
	@ $(CSHELL) -cf 'setenv CNS_NOFASTM ""; \
	                 source $(CNS_INST)/arch_env; cd $(CNS_SOURCE); \
	                 make -f Makefile.proto -k remove-exe OBJS=none; \
	                 make -k cns_solve'

static-fastm: valid_solve valid_comp
	@ $(CSHELL) -cf 'setenv CNS_FASTM_STATIC ""; \
	                 source $(CNS_INST)/arch_env; cd $(CNS_SOURCE); \
	                 make -f Makefile.proto -k remove-exe OBJS=none; \
	                 make -k cns_solve'

run_inputs:
	@ for dir in `/bin/ls -d inputs/*/run`; \
	  do echo "running input files in $$dir"; \
	    cd $(CNS_SOLVE)/$$dir; \
	    make link; \
	    make; \
	    make action=diff; \
	  done

win32_inputs:
	@ if [ -f windows-nt/run_inputs.bat ]; \
	    then /bin/rm -f windows-nt/run_inputs.bat; \
	  fi
	@ if [ -d $(CNS_SOLVE)/inputs ]; \
	  then \
	    for dir in `/bin/ls -d inputs/*/run`; \
	    do echo 'cd %CNS_SOLVE%/'"$$dir" | tr '/' '\\' \
	                 >> $(CNS_SOLVE)/windows-nt/run_inputs.bat; \
	      cd $(CNS_SOLVE)/$$dir; \
	      make link; \
	      make action=create; \
	      make action=echo >> $(CNS_SOLVE)/windows-nt/run_inputs.bat; \
	    done; \
	  fi

run_tests: valid_solve
	@ if [ ! -d $(CNS_INST)/test ]; \
	    then echo "making test directory $(CNS_INST)/test"; \
	         mkdir $(CNS_INST)/test; \
	  fi
	@ cd $(CNS_INST)/test; \
	  $(CNS_SOLVE)/bin/run_tests -link; \
	  $(CNS_SOLVE)/bin/run_tests *.inp; \
	  $(CNS_SOLVE)/bin/run_tests -tidy; \
	  $(CNS_SOLVE)/bin/run_diffs *.inp

relink: valid_solve
	@ cd $(CNS_SOURCE); make -k relink

makefiles: valid_solve
	@ cd $(CNS_SOURCE); make -f Makefile.proto -k makefiles OBJS=none

machinefile: valid_solve
	@ cd $(CNS_SOURCE); make -f Makefile.proto -k machinefile OBJS=none

tidy: valid_solve
	-@ cd $(CNS_SOURCE); make -k tidy

space: valid_solve
	-@ cd $(CNS_SOURCE); make -k space

clean: valid_solve
	-@ cd $(CNS_SOURCE); make -k clean

syntax-manual:
	@ if [ -n "$(CNS_SOLVE)" ]; \
	  then if [ -d $(CNS_SOLVE)/helplib ]; \
	      then echo "generating online syntax manual and index"; \
	        $(CNS_SOLVE)/bin/syntax_convert \
	          -source=$(CNS_SOLVE)/helplib \
	          -dest=$(CNS_SOLVE)/doc/html/syntax_manual/html; \
	      else echo "run make helplib first"; \
	      fi; \
	  else \
	    echo '$$CNS_SOLVE undefined - source cns_solve_env'; \
	    exit 1; \
	  fi

valid_solve:
	@ if [ -z "$(CNS_SOLVE)" ]; \
	    then echo '$$CNS_SOLVE undefined - source cns_solve_env'; \
	    exit 1; \
	  fi
	@ if [ "$(CNS_SOLVE)" = "_CNSsolve_location_" ]; \
	    then echo '$$CNS_SOLVE needs to be set'; \
	         echo "modify CNS_SOLVE in cns_solve_env"; \
	    exit 1; \
	  fi
	@ if [ -z "$(CNS_INST)" ]; \
	    then echo "machine directory undefined"; \
	         echo "source cns_solve_env or make install required"; \
	    exit 1; \
	  fi
	@ if [ "$(CNS_SOLVE)/$(CNS_ARCH)" != "$(CNS_INST)" ]; \
	    then \
	      echo "CNS installation does not match current machine type"; \
	      echo "CNS_ARCH is: $(CNS_ARCH)"; \
	      echo "CNS is installed: $(CNS_INST)"; \
	    exit 1; \
	  fi

valid_comp:
	@ if [ ! -r "$(CNS_INST)/arch_env" ]; \
	    then echo '$$CNS_INST/arch_env missing - check installation'; \
	    exit 1; \
	  fi

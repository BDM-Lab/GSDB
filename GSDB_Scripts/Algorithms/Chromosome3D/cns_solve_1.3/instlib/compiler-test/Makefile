###############################################################################
# 									      #
#	Makefile for CNSsolve						      #
#	Paul Adams  22-10-98						      #
#       copyright Yale University					      #
# 									      #
###############################################################################

.SUFFIXES:
.SUFFIXES: .f .c .o

RM=/bin/rm

SHELL=/bin/sh

# Fortran compiler test
fortran-test:
	@ cd $(CNS_INST)/compiler-test; \
	  echo "compiling: $(F77) $(F77FLAGS)"; \
	  $(F77) $(F77FLAGS) -c test_f.f; \
	  echo "linking: $(LD) $(LDFLAGS)"; \
	  $(LD) -o test_f test_f.o $(LDFLAGS); \
	  test=`./test_f | tr -d " "`; \
	  $(RM) -f test_f.o test_f; \
	  if [ X$$test = XHello_FORTRAN_world ]; \
	  then \
	    echo 'Fortran compiler passes test'; \
	  else \
	    echo '***** ERROR: problem with Fortran compiler *****'; \
	    exit 2; \
	  fi; \

# C compiler test
c-test:
	@ cd $(CNS_INST)/compiler-test; \
	  echo "compiling: $(CC) $(CCFLAGS)"; \
	  $(CC) $(CCFLAGS) -o test_c test_c.c; \
	  test=`./test_c | tr -d " "`; \
	  $(RM) -f test_c; \
	  if [ X$$test = XHello_C_world ]; \
	  then \
	    echo 'C compiler passes test'; \
	  else \
	    echo '***** ERROR: problem with C compiler *****'; \
	    exit 3; \
	  fi

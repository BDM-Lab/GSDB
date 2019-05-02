###############################################################################
#									      #
#       Makefile for CNSsolve utility programs				      #
#       Paul Adams 22-10-98						      #
#	copyright Yale University					      #
#									      #
###############################################################################

SHELL=/bin/sh

FTARGETS= PSmapx PSmapy PSmapz

CTARGETS= to_cns cluster_struc

LTARGETS= refloat

default: $(FTARGETS) $(CTARGETS) $(LTARGETS)

utils:
	make relink
	make default
	make clean

.f:
	$(F77) -o $@ $(F77FLAGS) $(@).f $(F77LINK)

.c:
	$(CC) -o $@ $(CCFLAGS) $(@).c $(CCLINK)

.cpp:
	$(CPP) -o $@ $(CCFLAGS) $(@).cpp $(CCLINK)

.l:
	lex $(@).l
	$(CC) $(CCFLAGS) -o $@ lex.yy.c $(CCLINK) -l$(LEXLIB)

relink:
	@ cd $(CNS_INST)/utils; touch x; rm -f x `ls * | grep -v Makefile`
	@ cd $(CNS_SOLVE)/utils; \
	  for file in `/bin/ls *`; \
	    do if [ ! -f $(CNS_INST)/utils/$$file ]; \
	      then cd $(CNS_INST)/utils; \
	           ln -s ../../utils/$$file .; \
	      fi; \
	  done

clean:
	@ rm -f *.o *.f *.c *.l

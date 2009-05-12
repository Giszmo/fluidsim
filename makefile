#==============================================================================
# G R A P H I K  -  P R O G R A M M I E R U N G   U N D   A N W E N D U N G E N
#------------------------------------------------------------------------------
#                         Diplomarbeit: Wassersimulator
#------------------------------------------------------------------------------
# $Id: Makefile,v 1.1 2004/04/16 12:58:41 DOMAIN-I15+prkipfer Exp $	
#------------------------------------------------------------------------------
# Bearbeitet von: Leo Wandersleb, e-mail: Leo.Wandersleb@gmx.de
#==============================================================================

# check which system we're running on an setup variables
SHELL :=/bin/sh
HOSTNAME  := $(shell hostname)
ifeq "$(VARIANT)" ""
override VARIANT := std
endif
OS := $(shell uname -s)
ifeq "$(OS)" "IRIX64"
OS := IRIX
endif
ifeq "$(OSVER)" ""
OSVER := $(shell uname -r)
endif
OSMVER := $(firstword $(subst ., ,$(OSVER)))
MK := $(OS).$(OSMVER)_$(COMP_VER)$(VARIANT)

# choose compiler
CC := gcc -DPCCTS_USE_NAMESPACE_STD
CFLAGS := -g -DVERBOSE 
ifeq "$(VARIANT)" "dbg"
CFLAGS += -DTRACE -Wall -pg
endif
ifeq "$(VARIANT)" "opt"
CFLAGS := -O3 -fomit-frame-pointer -ffast-math -fexpensive-optimizations -funroll-loops -fprefetch-loop-arrays -mmmx -march=x86-64 -mno-sse2 -msse -mfpmath=sse -mfpmath=sse,387
#CFLAGS := -O3 -fomit-frame-pointer -ffast-math -fexpensive-optimizations -funroll-loops -fprefetch-loop-arrays -mmmx -march=pentium4 -mno-sse2 -msse -mfpmath=sse -mfpmath=sse,387
#CFLAGS := -O3 -march=pentium4 -pipe -ffast-math -fPIC -mno-sse2 -mmmx -msse -mfpmath=sse,387 -falign-functions=4 -fomit-frame-pointer
#CFLAGS := -O3 -fomit-frame-pointer -ffast-math -fexpensive-optimizations -funroll-loops -fprefetch-loop-arrays -mmmx -march=pentium4 -msse2 -mfpmath=sse
#CFLAGS := -O3 -fomit-frame-pointer -ffast-math -fexpensive-optimizations -funroll-loops -fprefetch-loop-arrays -mmmx -m3dnow -march=athlon-xp -msse -mfpmath=sse
endif
PCCTSDIR = /usr/lib/pccts

ifeq "$(OS)" "IRIX"
CC := CC -n32
CFLAGS := -g -DVERBOSE -MDupdate Makefile.d
ifeq "$(VARIANT)" "dbg"
CFLAGS += -DTRACE -fullwarn
endif
ifeq "$(VARIANT)" "opt"
CFLAGS := -O3
endif
PCCTSDIR = /share/pckg/pccts
endif

# setup misc
INCS = -I/share/GL/include/ -I/share/GL/include/Inventor  -I/usr/X11R6/include
MAKE = make

LIBS     =  -lboost_thread -lm -L/share/GL/lib -L/usr/X11R6/lib -lglut -L/usr/X11/lib -lGLU -lGL 

SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))
EXEC := main.$(MK)



.C.o:
	$(CC) -c $(CFLAGS) $(INCS) $<

$(EXEC):  $(PCCTSOBJ) $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(PCCTSOBJ) -o $(EXEC) $(INCS) $(LIBS) 

# ----- PCCTS RULES -----
# ANTLR: generate <G>P.cpp, <G>.cpp <G>.dlg <G>.h from <G>.g
# DLG: generate   <G>L.cpp and <G>L.h from <G>.dlg
%L.cpp %P.cpp %.cpp %.dlg %.h: %.g
	antlr -fl $(patsubst %.g,%.dlg,$<) -ft $(patsubst %.g,%.h,$<) -CC -ga -gl $<
	dlg -C2 -CC -cl $(patsubst %.g,%L,$<) -ga $(patsubst %.g,%.dlg,$<)

# Compile PCCTS Base classes
%.o: $(PCCTSDIR)/h/%.cpp
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

# Compile PCCTS generated cpp files
%.o: %.cpp
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

depend:	$(SRCS)
	makedepend $(INCS) $(SRCS) -fMakefile.d

#
# target variants
#

dbg:
	$(MAKE) VARIANT=dbg

opt:
	$(MAKE) VARIANT=opt

#
# clearing all object files and the executable
#
clean:
	\rm -f $(OBJS) y.tab.c ;

realclean:
	$(MAKE) clean
	\rm -f $(EXEC) *~ \#* lint.output main.pure core Makefile.bak Makefile~ Makefile.d* $(GENERATED)
	\rm -rf ii_files doc
	touch Makefile.d

#
# some special things
#

doc:
	doxygen doxygen.conf

lint:
	lint $(INCS) $(SRCS) > lint.output

tags:	
	etags -S -C $(SRCS) *.hh parser.g Makefile

.PRECIOUS: %.cpp

# DO NOT DELETE THIS LINE -- make depend depends on it.
